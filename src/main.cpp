#include "motorsim/ingest.hpp"
#include "motorsim/io_csv.hpp"
#include "motorsim/io_vtk.hpp"
#include "motorsim/mechanical.hpp"
#include "motorsim/probes.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"
#include "motorsim/circuit.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cctype>
#include <filesystem>
#include <deque>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
#include <string_view>
#include <mutex>
#include <system_error>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <thread>

namespace {

void printUsage() {
    std::cout << "Usage: motor_sim [--scenario PATH] [--solve] [--write-midline]"
                 " [--list-outputs] [--outputs IDs] [--max-iters N] [--tol VALUE]"
                 " [--omega VALUE] [--parallel-frames] [--vtk-series PATH]"
                " [--solver {sor|cg|harmonic}] [--pc {none|jacobi|ssor}]"
                " [--harmonic-frequency HZ] [--harmonic-omega RAD_S]"
                " [--warm-start] [--no-warm-start] [--use-prolongation]"
                 " [--no-prolongation] [--coarse-nx N] [--coarse-ny N]"
                 " [--progress-every SEC] [--snapshot-every N]"
                 " [--progress-history PATH] [--quiet] [--no-quiet]\n";
}

std::vector<std::string> splitCommaSeparated(const std::string& text) {
    std::vector<std::string> items;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) {
            items.push_back(item);
        }
    }
    return items;
}

void ensureParentDirectory(const std::filesystem::path& path) {
    const std::filesystem::path parent = path.parent_path();
    if (!parent.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(parent, ec);
    }
}

struct OutputFilter {
    bool restrict{false};
    bool skip{false};
    const std::unordered_set<std::string>* requested{nullptr};
};

struct FrameRunResult {
    bool solverSuccess{false};
    bool outputError{false};
    std::string stdoutLog;
    std::string stderrLog;
    std::unordered_map<std::string, double> backEmfFluxes;
    std::unordered_map<std::string, motorsim::StressTensorResult> torqueSamples;
    std::vector<double> AzSolution;
    std::vector<motorsim::ProgressSample> progressHistory;
    bool hasMagneticCoenergy{false};
    double magneticCoenergy{0.0};
    struct VtkSeriesFrame {
        std::string id;
        std::string path;
        double time{0.0};
    };
    struct BoreSample {
        std::string id;
        double time{0.0};
        double bx{0.0};
        double by{0.0};
        double magnitude{0.0};
    };
    std::vector<VtkSeriesFrame> vtkSeriesFrames;
    std::vector<BoreSample> boreSamples;
};

struct FrameProgressSnapshot {
    std::size_t frameIndex{0};
    double frameTime{0.0};
    bool hasTime{false};
    bool active{false};
    bool completed{false};
    bool success{false};
    std::size_t lastIteration{0};
    double lastResidual{0.0};
    std::chrono::steady_clock::time_point lastUpdate{};
    double targetTol{1e-6};
    double etaSeconds{-1.0};
};

class FrameProgressState {
public:
    FrameProgressState(std::size_t index, double timeValue, bool hasTimeValue, double tol)
        : frameIndex(index), frameTime(timeValue), hasTime(hasTimeValue), targetTol(tol),
          lastUpdate(Clock::now()) {}

    void markStarted() {
        std::lock_guard<std::mutex> lock(mutex);
        active = true;
        completed = false;
        success = false;
        lastIteration = 0;
        lastResidual = 0.0;
        lastUpdate = Clock::now();
        history.clear();
    }

    void update(std::size_t iteration, double residual) {
        std::lock_guard<std::mutex> lock(mutex);
        active = true;
        lastIteration = iteration;
        lastResidual = residual;
        const auto now = Clock::now();
        lastUpdate = now;
        history.push_back({iteration, residual, now});
        if (history.size() > maxHistory) {
            history.pop_front();
        }
    }

    void markFinished(std::size_t iteration, double residual, bool didSucceed) {
        std::lock_guard<std::mutex> lock(mutex);
        active = false;
        completed = true;
        success = didSucceed;
        lastIteration = iteration;
        lastResidual = residual;
        lastUpdate = Clock::now();
        history.push_back({iteration, residual, lastUpdate});
        if (history.size() > maxHistory) {
            history.pop_front();
        }
    }

    FrameProgressSnapshot snapshot() const {
        std::lock_guard<std::mutex> lock(mutex);
        FrameProgressSnapshot snap{};
        snap.frameIndex = frameIndex;
        snap.frameTime = frameTime;
        snap.hasTime = hasTime;
        snap.active = active;
        snap.completed = completed;
        snap.success = success;
        snap.lastIteration = lastIteration;
        snap.lastResidual = lastResidual;
        snap.lastUpdate = lastUpdate;
        snap.targetTol = targetTol;
        snap.etaSeconds = estimateEtaSeconds();
        return snap;
    }

    const std::size_t frameIndex;
    const double frameTime;
    const bool hasTime;
    const double targetTol;

private:
    using Clock = std::chrono::steady_clock;

    struct HistorySample {
        std::size_t iteration{0};
        double residual{0.0};
        Clock::time_point time{};
    };

    double estimateEtaSeconds() const {
        if (history.size() < 2 || !(targetTol > 0.0)) {
            return -1.0;
        }
        const HistorySample& latest = history.back();
        if (latest.residual <= targetTol) {
            return 0.0;
        }
        for (auto it = history.rbegin() + 1; it != history.rend(); ++it) {
            if (it->iteration == latest.iteration) {
                continue;
            }
            if (!(it->residual > 0.0) || !(latest.residual > 0.0)) {
                continue;
            }
            const std::size_t deltaIter = latest.iteration - it->iteration;
            if (deltaIter == 0) {
                continue;
            }
            const double deltaLog = std::log(it->residual) - std::log(latest.residual);
            if (!(deltaLog > 0.0)) {
                continue;
            }
            const double decayPerIter = deltaLog / static_cast<double>(deltaIter);
            if (!(decayPerIter > 0.0)) {
                continue;
            }
            const double remainingLog = std::log(latest.residual) - std::log(targetTol);
            if (!(remainingLog > 0.0)) {
                return 0.0;
            }
            const double remainingIter = remainingLog / decayPerIter;
            if (!(remainingIter > 0.0)) {
                continue;
            }
            const double deltaTime = std::chrono::duration<double>(latest.time - it->time).count();
            if (!(deltaTime > 0.0)) {
                continue;
            }
            const double secondsPerIter = deltaTime / static_cast<double>(deltaIter);
            return remainingIter * secondsPerIter;
        }
        return -1.0;
    }

    mutable std::mutex mutex;
    bool active{false};
    bool completed{false};
    bool success{false};
    std::size_t lastIteration{0};
    double lastResidual{0.0};
    Clock::time_point lastUpdate;
    std::deque<HistorySample> history;
    static constexpr std::size_t maxHistory = 6;
};

class ProgressPrinter {
public:
    ProgressPrinter(std::vector<std::shared_ptr<FrameProgressState>> states,
                    double intervalSec,
                    bool quiet)
        : states_(std::move(states)), intervalSec_(intervalSec), quiet_(quiet), start_(Clock::now()) {
        if (!quiet_ && !states_.empty() && intervalSec_ > 0.0) {
            worker_ = std::thread([this]() { run(); });
        }
    }

    ~ProgressPrinter() { stop(); }

    void stop() {
        const bool wasRunning = !stopped_.exchange(true);
        if (wasRunning && worker_.joinable()) {
            worker_.join();
        }
        if (!quiet_ && !states_.empty() && !finalPrinted_.exchange(true)) {
            printStatus(true);
        }
    }

private:
    using Clock = std::chrono::steady_clock;

    void run() {
        while (!stopped_.load()) {
            std::this_thread::sleep_for(std::chrono::duration<double>(intervalSec_));
            if (stopped_.load()) {
                break;
            }
            const bool finished = printStatus(false);
            if (finished) {
                finalPrinted_.store(true);
                break;
            }
        }
    }

    static std::string formatElapsed(double seconds) {
        const int totalSeconds = static_cast<int>(seconds);
        const int minutes = totalSeconds / 60;
        const int remSeconds = totalSeconds % 60;
        std::ostringstream oss;
        oss << std::setfill('0') << std::setw(2) << minutes << ':' << std::setw(2) << remSeconds;
        return oss.str();
    }

    bool printStatus(bool force) {
        if (states_.empty() || quiet_) {
            return true;
        }

        std::vector<FrameProgressSnapshot> snapshots;
        snapshots.reserve(states_.size());
        for (const auto& state : states_) {
            snapshots.push_back(state->snapshot());
        }

        std::size_t completed = 0;
        std::size_t active = 0;
        std::size_t failed = 0;
        for (const auto& snap : snapshots) {
            if (snap.completed) {
                ++completed;
                if (!snap.success) {
                    ++failed;
                }
            }
            if (snap.active) {
                ++active;
            }
        }

        if (!force && completed == 0 && active == 0) {
            return false;
        }

        const std::size_t total = snapshots.size();
        const std::size_t queued = total > (completed + active) ? total - completed - active : 0;

        const double elapsed = std::chrono::duration<double>(Clock::now() - start_).count();

        std::ostringstream oss;
        oss << '[' << formatElapsed(elapsed) << "] frames: " << completed << '/' << total
            << " done | active: " << active << " | queued: " << queued;
        if (failed > 0) {
            oss << " | failed: " << failed;
        }
        oss << '\n';

        if (active > 0) {
            oss << "        ";
            bool first = true;
            for (const auto& snap : snapshots) {
                if (!snap.active) {
                    continue;
                }
                if (!first) {
                    oss << " | ";
                }
                first = false;
                oss << 'F' << snap.frameIndex << " iter=" << std::setw(6) << snap.lastIteration
                    << " relRes=";
                oss << std::scientific << std::setprecision(3) << snap.lastResidual;
                oss << std::defaultfloat;
                double progressRatio = 0.0;
                if (snap.lastIteration > 0 && snap.lastResidual > 0.0 && snap.targetTol > 0.0) {
                    progressRatio = std::clamp(snap.targetTol / snap.lastResidual, 0.0, 1.0);
                }
                oss << " prog=" << std::fixed << std::setprecision(0) << progressRatio * 100.0 << '%';
                oss << std::defaultfloat;
                if (snap.etaSeconds >= 0.0) {
                    oss << " eta≈" << std::fixed << std::setprecision(1) << snap.etaSeconds << 's';
                    oss << std::defaultfloat;
                }
            }
            oss << '\n';
        }

        std::cout << oss.str();
        std::cout.flush();

        return completed == total;
    }

    std::vector<std::shared_ptr<FrameProgressState>> states_;
    double intervalSec_{0.0};
    bool quiet_{false};
    Clock::time_point start_;
    std::atomic<bool> stopped_{false};
    std::atomic<bool> finalPrinted_{false};
    std::thread worker_{};
};

class FrameProgressGuard {
public:
    explicit FrameProgressGuard(FrameProgressState* state) : state_(state) {
        if (state_) {
            state_->markStarted();
        }
    }

    FrameProgressGuard(const FrameProgressGuard&) = delete;
    FrameProgressGuard& operator=(const FrameProgressGuard&) = delete;

    ~FrameProgressGuard() {
        if (state_) {
            state_->markFinished(finalIterations_, finalResidual_, success_);
        }
    }

    void markFinal(std::size_t iterations, double residual, bool success) {
        finalIterations_ = iterations;
        finalResidual_ = residual;
        success_ = success;
    }

private:
    FrameProgressState* state_{nullptr};
    std::size_t finalIterations_{0};
    double finalResidual_{0.0};
    bool success_{false};
};

std::filesystem::path makeFramePath(const std::string& basePath, bool timelineActive,
                                    std::size_t frameIndex, std::size_t digits) {
    if (!timelineActive) {
        return std::filesystem::path(basePath);
    }
    std::filesystem::path base(basePath);
    const std::filesystem::path dir = base.parent_path();
    std::ostringstream stem;
    const std::string stemName = base.stem().string();
    constexpr std::string_view frameTag{"_frame"};
    if (stemName.size() >= frameTag.size() &&
        stemName.compare(stemName.size() - frameTag.size(), frameTag.size(), frameTag) == 0) {
        stem << stemName << '_';
    } else {
        stem << stemName << "_frame_";
    }
    stem << std::setfill('0') << std::setw(static_cast<int>(digits)) << frameIndex;
    const std::string extension = base.extension().string();
    return dir / (stem.str() + extension);
}

bool pointInPolygon(double x, double y, const std::vector<double>& xs, const std::vector<double>& ys) {
    const std::size_t count = xs.size();
    if (count == 0 || count != ys.size()) {
        return false;
    }
    bool inside = false;
    for (std::size_t i = 0, j = count - 1; i < count; j = i++) {
        const double xi = xs[i];
        const double yi = ys[i];
        const double xj = xs[j];
        const double yj = ys[j];
        const bool intersect = ((yi > y) != (yj > y)) &&
                               (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (intersect) {
            inside = !inside;
        }
    }
    return inside;
}

motorsim::VtkOutlineLoop makeRectangleLoop(motorsim::VtkOutlineLoop::Kind kind,
                                           std::string label,
                                           double minX,
                                           double maxX,
                                           double minY,
                                           double maxY) {
    motorsim::VtkOutlineLoop loop;
    loop.kind = kind;
    loop.label = std::move(label);
    loop.xs = {minX, maxX, maxX, minX};
    loop.ys = {minY, minY, maxY, maxY};
    return loop;
}

std::vector<motorsim::VtkOutlineLoop> buildOutlineLoops(const motorsim::ScenarioSpec& spec) {
    using motorsim::VtkOutlineLoop;
    std::vector<VtkOutlineLoop> loops;

    std::vector<std::string> polygonGroups(spec.polygons.size());
    std::vector<std::string> polygonNames(spec.polygons.size());
    std::vector<std::string> magnetGroups(spec.magnetRegions.size());
    std::vector<std::string> magnetNames(spec.magnetRegions.size());
    std::vector<std::string> wireGroups(spec.wires.size());
    std::vector<std::string> wireNames(spec.wires.size());

    for (std::size_t rotorIdx = 0; rotorIdx < spec.rotors.size(); ++rotorIdx) {
        const auto& rotor = spec.rotors[rotorIdx];
        std::string rotorLabel = rotor.name;
        if (rotorLabel.empty()) {
            rotorLabel = "rotor_" + std::to_string(rotorIdx);
        }
        for (std::size_t local = 0; local < rotor.polygonIndices.size(); ++local) {
            const std::size_t idx = rotor.polygonIndices[local];
            if (idx >= polygonGroups.size()) {
                continue;
            }
            polygonGroups[idx] = rotorLabel;
            polygonNames[idx] = rotorLabel + ":polygon_" + std::to_string(local);
        }
        for (std::size_t local = 0; local < rotor.magnetIndices.size(); ++local) {
            const std::size_t idx = rotor.magnetIndices[local];
            if (idx >= magnetGroups.size()) {
                continue;
            }
            magnetGroups[idx] = rotorLabel;
            magnetNames[idx] = rotorLabel + ":magnet_" + std::to_string(local);
        }
        for (std::size_t local = 0; local < rotor.wireIndices.size(); ++local) {
            const std::size_t idx = rotor.wireIndices[local];
            if (idx >= wireGroups.size()) {
                continue;
            }
            wireGroups[idx] = rotorLabel;
            wireNames[idx] = rotorLabel + ":wire_" + std::to_string(local);
        }
    }

    const double minX = spec.originX;
    const double maxX = spec.originX + static_cast<double>(spec.nx - 1) * spec.dx;
    const double minY = spec.originY;
    const double maxY = spec.originY + static_cast<double>(spec.ny - 1) * spec.dy;

    loops.push_back(makeRectangleLoop(VtkOutlineLoop::Kind::Domain, "domain", minX, maxX, minY, maxY));

    for (std::size_t i = 0; i < spec.polygons.size(); ++i) {
        const auto& poly = spec.polygons[i];
        if (poly.xs.size() != poly.ys.size() || poly.xs.size() < 3) {
            continue;
        }
        VtkOutlineLoop loop;
        loop.kind = VtkOutlineLoop::Kind::Material;
        if (!polygonNames[i].empty()) {
            loop.label = polygonNames[i];
        } else {
            loop.label = "material_polygon_" + std::to_string(i);
        }
        loop.groupLabel = polygonGroups[i];
        loop.xs = poly.xs;
        loop.ys = poly.ys;
        loops.push_back(std::move(loop));
    }

    for (std::size_t i = 0; i < spec.magnetRegions.size(); ++i) {
        const auto& magnet = spec.magnetRegions[i];
        VtkOutlineLoop loop;
        loop.kind = VtkOutlineLoop::Kind::Magnet;
        if (!magnetNames[i].empty()) {
            loop.label = magnetNames[i];
        } else {
            loop.label = "magnet_" + std::to_string(i);
        }
        loop.groupLabel = magnetGroups[i];
        if (magnet.xs.size() == magnet.ys.size() && magnet.xs.size() >= 3) {
            loop.xs = magnet.xs;
            loop.ys = magnet.ys;
            loops.push_back(std::move(loop));
        } else {
            loops.push_back(makeRectangleLoop(loop.kind, loop.label, magnet.min_x, magnet.max_x, magnet.min_y,
                                              magnet.max_y));
        }
    }

    const double twoPi = 6.28318530717958647692;

    for (std::size_t i = 0; i < spec.wires.size(); ++i) {
        const auto& wire = spec.wires[i];
        if (!(wire.radius > 0.0)) {
            continue;
        }
        VtkOutlineLoop loop;
        loop.kind = VtkOutlineLoop::Kind::Wire;
        if (!wireNames[i].empty()) {
            loop.label = wireNames[i];
        } else {
            loop.label = "wire_" + std::to_string(i);
        }
        loop.groupLabel = wireGroups[i];
        constexpr std::size_t segments = 64;
        loop.xs.reserve(segments);
        loop.ys.reserve(segments);
        for (std::size_t s = 0; s < segments; ++s) {
            const double theta = (twoPi * static_cast<double>(s)) / static_cast<double>(segments);
            loop.xs.push_back(wire.x + wire.radius * std::cos(theta));
            loop.ys.push_back(wire.y + wire.radius * std::sin(theta));
        }
        loops.push_back(std::move(loop));
    }

    for (std::size_t i = 0; i < spec.currentRegions.size(); ++i) {
        const auto& region = spec.currentRegions[i];
        if (region.xs.size() != region.ys.size() || region.xs.size() < 3) {
            continue;
        }
        VtkOutlineLoop loop;
        loop.kind = VtkOutlineLoop::Kind::CurrentRegion;
        if (!region.id.empty()) {
            loop.label = region.id;
        } else {
            loop.label = "current_region_" + std::to_string(i);
        }
        loop.groupLabel = region.phase;
        loop.xs = region.xs;
        loop.ys = region.ys;
        loops.push_back(std::move(loop));
    }

    return loops;
}

bool tryParseSolverKind(const std::string& text, motorsim::SolverKind& out) {
    std::string lower;
    lower.reserve(text.size());
    for (unsigned char ch : text) {
        lower.push_back(static_cast<char>(std::tolower(ch)));
    }
    if (lower == "sor") {
        out = motorsim::SolverKind::SOR;
        return true;
    }
    if (lower == "cg") {
        out = motorsim::SolverKind::CG;
        return true;
    }
    if (lower == "harmonic" || lower == "harmonic_cg") {
        out = motorsim::SolverKind::Harmonic;
        return true;
    }
    return false;
}

bool tryParsePreconditioner(const std::string& text, motorsim::PreconditionerKind& out) {
    std::string lower;
    lower.reserve(text.size());
    for (unsigned char ch : text) {
        lower.push_back(static_cast<char>(std::tolower(ch)));
    }
    if (lower == "none") {
        out = motorsim::PreconditionerKind::None;
        return true;
    }
    if (lower == "jacobi") {
        out = motorsim::PreconditionerKind::Jacobi;
        return true;
    }
    if (lower == "ssor") {
        out = motorsim::PreconditionerKind::SSOR;
        return true;
    }
    return false;
}

class FrameProgressSinkImpl : public motorsim::ProgressSink {
public:
    using SnapshotProvider = std::function<std::optional<std::string>(const motorsim::ProgressSample&)>;

    FrameProgressSinkImpl(FrameProgressState* state,
                          SnapshotProvider provider,
                          std::vector<motorsim::ProgressSample>* history)
        : state_(state), provider_(std::move(provider)), history_(history) {}

    bool onProgress(const motorsim::ProgressSample& sample) override {
        if (state_) {
            state_->update(sample.iter, sample.relResidual);
        }
        if (history_) {
            history_->push_back(sample);
        }
        return true;
    }

    std::optional<std::string> requestFieldDump(const motorsim::ProgressSample& sample) override {
        if (!provider_) {
            return std::nullopt;
        }
        return provider_(sample);
    }

private:
    FrameProgressState* state_{nullptr};
    SnapshotProvider provider_{};
    std::vector<motorsim::ProgressSample>* history_{nullptr};
};

FrameRunResult solveFrame(motorsim::ScenarioFrame& frame, const motorsim::SolveOptions& options,
                          const OutputFilter& filter, bool timelineActive, std::size_t frameDigits,
                          bool writeMidline, FrameProgressState* progressState,
                          const std::vector<double>* warmStartGuess,
                          motorsim::CircuitSimulator* circuitSim,
                          const std::vector<double>* transientPrevAz,
                          double transientDt) {
    FrameRunResult result{};
    std::ostringstream out;
    std::ostringstream err;
    FrameProgressGuard progressGuard(progressState);

    if (circuitSim) {
        circuitSim->apply_currents(frame);
    }

    motorsim::Grid2D grid(frame.spec.nx, frame.spec.ny, frame.spec.dx, frame.spec.dy);
    try {
        motorsim::rasterizeScenarioToGrid(frame.spec, grid);
    } catch (const std::exception& ex) {
        err << "Frame " << frame.index << ": rasterisation error: " << ex.what() << '\n';
        progressGuard.markFinal(0, 0.0, false);
        result.stderrLog = err.str();
        return result;
    }

    motorsim::SolveOptions solveOptions = options;

    const bool transientActive = frame.spec.transient.enabled && transientPrevAz != nullptr && (transientDt > 0.0);

    FrameProgressSinkImpl::SnapshotProvider snapshotProvider;
    if (progressState && options.snapshotEveryIters > 0) {
        snapshotProvider = [timelineActive, frameDigits, frame](const motorsim::ProgressSample& sample) {
            std::ostringstream name;
            name << "outputs/snapshots/frame_";
            if (timelineActive) {
                name << std::setfill('0') << std::setw(static_cast<int>(frameDigits)) << frame.index;
            } else {
                name << std::setw(3) << std::setfill('0') << frame.index;
            }
            name << "_iter_" << std::setw(6) << std::setfill('0') << sample.iter << ".csv";
            return std::optional<std::string>(name.str());
        };
    }
    FrameProgressSinkImpl sink(progressState, snapshotProvider, &result.progressHistory);

    motorsim::InitialGuess guess{};
    const motorsim::InitialGuess* guessPtr = nullptr;
    if (!transientActive && warmStartGuess && warmStartGuess->size() == grid.Az.size()) {
        guess.Az0 = warmStartGuess;
        guessPtr = &guess;
    }

    std::vector<double> emptyPrev;
    const std::vector<double>& transientPrev = (transientActive && transientPrevAz && transientPrevAz->size() == grid.Az.size())
                                                   ? *transientPrevAz
                                                   : emptyPrev;

    if (transientActive && transientPrev.empty() && transientPrevAz && transientPrevAz->size() != grid.Az.size()) {
        emptyPrev.assign(grid.Az.size(), 0.0);
    }

    if (transientActive && transientPrevAz && transientPrevAz->size() == grid.Az.size()) {
        grid.Az = *transientPrevAz;
    }

    motorsim::SolveResult report{};
    if (transientActive) {
        report = motorsim::solveTransientStep(grid, solveOptions, transientDt, transientPrev,
                                              progressState ? &sink : nullptr);
    } else {
        report = motorsim::solveAz(grid, solveOptions, guessPtr, progressState ? &sink : nullptr);
    }
    if (!report.converged) {
        err << "Frame " << frame.index
            << ": solver did not converge. relResidual=" << report.relResidual << '\n';
        progressGuard.markFinal(report.iters, report.relResidual, false);
        result.stderrLog = err.str();
        return result;
    }

    motorsim::computeB(grid);
    if (circuitSim) {
        circuitSim->record_solved_frame(frame, grid);
    }

    const bool scenarioHasOutputs = !frame.spec.outputs.fieldMaps.empty() ||
                                    !frame.spec.outputs.lineProbes.empty() ||
                                    !frame.spec.outputs.probes.empty() ||
                                    !frame.spec.outputs.backEmfProbes.empty() ||
                                    !frame.spec.outputs.vtkSeries.empty() ||
                                    !frame.spec.outputs.boreProbes.empty();

    const auto shouldEmit = [&](const std::string& id) {
        if (filter.skip) {
            return false;
        }
        if (!filter.restrict) {
            return true;
        }
        if (!filter.requested) {
            return false;
        }
        return filter.requested->find(id) != filter.requested->end();
    };

    bool fieldMapNeedsH = false;
    bool fieldMapNeedsEnergy = false;
    bool lineProbeNeedsH = false;
    bool lineProbeNeedsEnergy = false;
    bool vtkSeriesNeedsH = false;
    bool vtkSeriesNeedsEnergy = false;
    bool torqueProbeNeedsEnergy = false;
    bool outputError = false;

    std::vector<motorsim::VtkOutlineLoop> outlineLoops;
    bool outlinesPrepared = false;
    bool outlinesWritten = false;

    const auto prepareOutlines = [&]() {
        if (!outlinesPrepared) {
            outlineLoops = buildOutlineLoops(frame.spec);
            outlinesPrepared = true;
        }
    };

    if (scenarioHasOutputs && !filter.skip) {
        for (const auto& request : frame.spec.outputs.fieldMaps) {
            if (!shouldEmit(request.id)) {
                continue;
            }
            if (request.format == "vti") {
                fieldMapNeedsH = true;
                fieldMapNeedsEnergy = true;
            }
            if (request.quantity == "H" || request.quantity == "BH") {
                fieldMapNeedsH = true;
            } else if (request.quantity == "energy_density") {
                fieldMapNeedsH = true;
                fieldMapNeedsEnergy = true;
            }
        }
        for (const auto& request : frame.spec.outputs.lineProbes) {
            if (!shouldEmit(request.id)) {
                continue;
            }
            if (request.quantity == "Hx" || request.quantity == "Hy" || request.quantity == "Hmag") {
                lineProbeNeedsH = true;
            } else if (request.quantity == "energy_density") {
                lineProbeNeedsH = true;
                lineProbeNeedsEnergy = true;
            }
        }

        for (const auto& request : frame.spec.outputs.probes) {
            if (!shouldEmit(request.id)) {
                continue;
            }
            if (request.quantity == motorsim::ScenarioSpec::Outputs::Probe::Quantity::Torque ||
                request.quantity == motorsim::ScenarioSpec::Outputs::Probe::Quantity::ForceAndTorque) {
                torqueProbeNeedsEnergy = true;
            }
        }

        for (const auto& request : frame.spec.outputs.vtkSeries) {
            if (!shouldEmit(request.id)) {
                continue;
            }
            if (request.includeEnergy) {
                vtkSeriesNeedsEnergy = true;
                vtkSeriesNeedsH = true;
            } else if (request.includeH) {
                vtkSeriesNeedsH = true;
            }
        }

    }

    const bool needHField = fieldMapNeedsH || fieldMapNeedsEnergy || lineProbeNeedsH || lineProbeNeedsEnergy ||
                            vtkSeriesNeedsH || vtkSeriesNeedsEnergy || torqueProbeNeedsEnergy;
    if (needHField) {
        motorsim::computeH(grid);
        if (torqueProbeNeedsEnergy) {
            try {
                result.magneticCoenergy =
                    motorsim::compute_magnetic_coenergy(grid, frame.spec.dx, frame.spec.dy);
                result.hasMagneticCoenergy = true;
            } catch (const std::exception& ex) {
                err << "Frame " << frame.index << ": failed to compute magnetic co-energy: " << ex.what() << '\n';
                outputError = true;
            }
        }
    }

    out << "Frame " << frame.index;
    if (timelineActive) {
        out << " (t=" << frame.time << " s)";
    }
    out << ": solved in " << report.iters << " iterations, relResidual=" << report.relResidual
        << '\n';

    if (scenarioHasOutputs && !filter.skip) {
        std::vector<double> gridXs;
        std::vector<double> gridYs;
        std::vector<double> gridBx;
        std::vector<double> gridBy;
        std::vector<double> gridBmag;
        std::vector<double> gridHx;
        std::vector<double> gridHy;
        std::vector<double> gridHmag;
        std::vector<double> gridEnergyDensity;
        const bool fieldMapRequiresHData = fieldMapNeedsH || fieldMapNeedsEnergy;
        bool fieldPrepared = false;
        const auto prepareFieldVectors = [&]() {
            if (fieldPrepared) {
                return;
            }
            const std::size_t total = frame.spec.nx * frame.spec.ny;
            gridXs.resize(total);
            gridYs.resize(total);
            gridBx.resize(total);
            gridBy.resize(total);
            gridBmag.resize(total);
            if (fieldMapRequiresHData) {
                gridHx.resize(total);
                gridHy.resize(total);
                gridHmag.resize(total);
                if (fieldMapNeedsEnergy) {
                    gridEnergyDensity.resize(total);
                }
            }
            std::size_t idx = 0;
            for (std::size_t j = 0; j < frame.spec.ny; ++j) {
                const double y = frame.spec.originY + static_cast<double>(j) * frame.spec.dy;
                for (std::size_t i = 0; i < frame.spec.nx; ++i) {
                    const double x = frame.spec.originX + static_cast<double>(i) * frame.spec.dx;
                    const std::size_t gridIdx = grid.idx(i, j);
                    gridXs[idx] = x;
                    gridYs[idx] = y;
                    gridBx[idx] = grid.Bx[gridIdx];
                    gridBy[idx] = grid.By[gridIdx];
                    gridBmag[idx] = std::hypot(grid.Bx[gridIdx], grid.By[gridIdx]);
                    if (fieldMapRequiresHData) {
                        const double hx = grid.Hx[gridIdx];
                        const double hy = grid.Hy[gridIdx];
                        gridHx[idx] = hx;
                        gridHy[idx] = hy;
                        gridHmag[idx] = std::hypot(hx, hy);
                        if (fieldMapNeedsEnergy) {
                            gridEnergyDensity[idx] =
                                0.5 * (grid.Bx[gridIdx] * hx + grid.By[gridIdx] * hy);
                        }
                    }
                    ++idx;
                }
            }
            fieldPrepared = true;
        };

        for (const auto& request : frame.spec.outputs.fieldMaps) {
            if (!shouldEmit(request.id)) {
                continue;
            }
            const std::filesystem::path outPath =
                makeFramePath(request.path, timelineActive, frame.index, frameDigits);
            ensureParentDirectory(outPath);
            try {
                if (request.format == "csv") {
                    prepareFieldVectors();
                    std::vector<motorsim::FieldColumnView> columns;
                    if (request.quantity == "B") {
                        columns = {{"Bx", &gridBx}, {"By", &gridBy}, {"Bmag", &gridBmag}};
                    } else if (request.quantity == "H") {
                        columns = {{"Hx", &gridHx}, {"Hy", &gridHy}, {"Hmag", &gridHmag}};
                    } else if (request.quantity == "BH") {
                        columns = {{"Bx", &gridBx}, {"By", &gridBy}, {"Bmag", &gridBmag},
                                   {"Hx", &gridHx}, {"Hy", &gridHy}, {"Hmag", &gridHmag}};
                    } else {
                        columns = {{"Bx", &gridBx},        {"By", &gridBy},
                                   {"Bmag", &gridBmag},   {"Hx", &gridHx},
                                   {"Hy", &gridHy},       {"Hmag", &gridHmag},
                                   {"EnergyDensity", &gridEnergyDensity}};
                    }
                    motorsim::write_csv_field_map(outPath.string(), gridXs, gridYs, columns);
                } else {
                    const bool includeH = true;
                    const bool includeEnergy = true;
                    motorsim::write_vti_field_map(outPath.string(), frame.spec.nx, frame.spec.ny,
                                                  frame.spec.originX, frame.spec.originY, frame.spec.dx,
                                                  frame.spec.dy, grid.Bx, grid.By, &grid.Hx, &grid.Hy,
                                                  includeH, includeEnergy);
                    if (!outlinesWritten) {
                        prepareOutlines();
                        if (!outlineLoops.empty()) {
                            std::filesystem::path outlinePath = outPath;
                            outlinePath.replace_extension();
                            outlinePath += "_outlines.vtp";
                            ensureParentDirectory(outlinePath);
                            motorsim::write_vtp_outlines(outlinePath.string(), outlineLoops);
                            out << "Frame " << frame.index << ": wrote geometry_outlines to " << outlinePath
                                << '\n';
                        }
                        outlinesWritten = true;
                    }
                }
                out << "Frame " << frame.index << ": wrote field_map '" << request.id
                    << "' to " << outPath << '\n';
            } catch (const std::exception& ex) {
                err << "Frame " << frame.index << ": failed to write field_map '" << request.id
                    << "': " << ex.what() << '\n';
                outputError = true;
            }
        }

        for (const auto& request : frame.spec.outputs.vtkSeries) {
            if (!shouldEmit(request.id)) {
                continue;
            }
            if (!request.includeB) {
                err << "Frame " << frame.index << ": vtk_series '" << request.id
                    << "' must include B field" << '\n';
                outputError = true;
                continue;
            }

            std::filesystem::path basePath;
            if (request.directory.empty()) {
                basePath = std::filesystem::path(request.basename + ".vti");
            } else {
                basePath = std::filesystem::path(request.directory) / (request.basename + ".vti");
            }
            const std::filesystem::path outPath =
                makeFramePath(basePath.string(), timelineActive, frame.index, frameDigits);
            ensureParentDirectory(outPath);
            try {
                const bool includeH = request.includeH || request.includeEnergy;
                const bool includeEnergy = request.includeEnergy;
                const std::vector<double>* hxPtr = includeH ? &grid.Hx : nullptr;
                const std::vector<double>* hyPtr = includeH ? &grid.Hy : nullptr;
                motorsim::write_vti_field_map(outPath.string(), frame.spec.nx, frame.spec.ny,
                                              frame.spec.originX, frame.spec.originY, frame.spec.dx,
                                              frame.spec.dy, grid.Bx, grid.By, hxPtr, hyPtr, includeH,
                                              includeEnergy);
                result.vtkSeriesFrames.push_back(
                    FrameRunResult::VtkSeriesFrame{request.id, outPath.string(), frame.time});
                out << "Frame " << frame.index << ": wrote vtk_series '" << request.id << "' to " << outPath
                    << '\n';
            } catch (const std::exception& ex) {
                err << "Frame " << frame.index << ": failed to write vtk_series '" << request.id
                    << "': " << ex.what() << '\n';
                outputError = true;
            }
        }

        for (const auto& request : frame.spec.outputs.boreProbes) {
            if (!shouldEmit(request.id)) {
                continue;
            }
            if (request.xs.size() != request.ys.size() || request.xs.size() < 3) {
                err << "Frame " << frame.index << ": bore probe '" << request.id
                    << "' polygon is invalid" << '\n';
                outputError = true;
                continue;
            }
            double sumBx = 0.0;
            double sumBy = 0.0;
            double count = 0.0;
            for (std::size_t j = 0; j < frame.spec.ny; ++j) {
                const double y = frame.spec.originY + static_cast<double>(j) * frame.spec.dy;
                for (std::size_t i = 0; i < frame.spec.nx; ++i) {
                    const double x = frame.spec.originX + static_cast<double>(i) * frame.spec.dx;
                    if (!pointInPolygon(x, y, request.xs, request.ys)) {
                        continue;
                    }
                    const std::size_t idx = grid.idx(i, j);
                    sumBx += grid.Bx[idx];
                    sumBy += grid.By[idx];
                    count += 1.0;
                }
            }
            if (!(count > 0.0)) {
                err << "Frame " << frame.index << ": bore probe '" << request.id
                    << "' contains no grid points" << '\n';
                outputError = true;
                continue;
            }
            const double avgBx = sumBx / count;
            const double avgBy = sumBy / count;
            const double magnitude = std::hypot(avgBx, avgBy);
            result.boreSamples.push_back(
                FrameRunResult::BoreSample{request.id, frame.time, avgBx, avgBy, magnitude});
            out << "Frame " << frame.index << ": bore probe '" << request.id << "' avg|B|=" << magnitude
                << '\n';
        }

        for (const auto& request : frame.spec.outputs.lineProbes) {
            if (!shouldEmit(request.id)) {
                continue;
            }

            std::vector<double> xs;
            std::vector<double> ys;
            std::vector<double> values;

            if (request.axis == "x") {
                const double coord = (request.value - frame.spec.originX) / frame.spec.dx;
                const auto idx = static_cast<long long>(std::llround(coord));
                if (idx < 0 || idx >= static_cast<long long>(frame.spec.nx)) {
                    err << "Frame " << frame.index << ": line_probe '" << request.id
                        << "' is outside the domain (x=" << request.value << ")\n";
                    outputError = true;
                    continue;
                }
                const double actual =
                    frame.spec.originX + static_cast<double>(idx) * frame.spec.dx;
                if (std::abs(actual - request.value) > 0.5 * frame.spec.dx + 1e-9) {
                    err << "Frame " << frame.index << ": line_probe '" << request.id
                        << "' does not align with grid column (requested x=" << request.value
                        << ", snapped to " << actual << ")\n";
                    outputError = true;
                    continue;
                }

                xs.resize(frame.spec.ny, actual);
                ys.resize(frame.spec.ny);
                values.resize(frame.spec.ny);
                for (std::size_t j = 0; j < frame.spec.ny; ++j) {
                    const double y = frame.spec.originY + static_cast<double>(j) * frame.spec.dy;
                    const std::size_t gridIdx = grid.idx(static_cast<std::size_t>(idx), j);
                    ys[j] = y;
                    if (request.quantity == "Bx") {
                        values[j] = grid.Bx[gridIdx];
                    } else if (request.quantity == "By") {
                        values[j] = grid.By[gridIdx];
                    } else if (request.quantity == "Bmag") {
                        values[j] = std::hypot(grid.Bx[gridIdx], grid.By[gridIdx]);
                    } else if (request.quantity == "Hx") {
                        values[j] = grid.Hx[gridIdx];
                    } else if (request.quantity == "Hy") {
                        values[j] = grid.Hy[gridIdx];
                    } else if (request.quantity == "Hmag") {
                        values[j] = std::hypot(grid.Hx[gridIdx], grid.Hy[gridIdx]);
                    } else {
                        values[j] =
                            0.5 * (grid.Bx[gridIdx] * grid.Hx[gridIdx] +
                                   grid.By[gridIdx] * grid.Hy[gridIdx]);
                    }
                }
            } else if (request.axis == "y") {
                const double coord = (request.value - frame.spec.originY) / frame.spec.dy;
                const auto idx = static_cast<long long>(std::llround(coord));
                if (idx < 0 || idx >= static_cast<long long>(frame.spec.ny)) {
                    err << "Frame " << frame.index << ": line_probe '" << request.id
                        << "' is outside the domain (y=" << request.value << ")\n";
                    outputError = true;
                    continue;
                }
                const double actual =
                    frame.spec.originY + static_cast<double>(idx) * frame.spec.dy;
                if (std::abs(actual - request.value) > 0.5 * frame.spec.dy + 1e-9) {
                    err << "Frame " << frame.index << ": line_probe '" << request.id
                        << "' does not align with grid row (requested y=" << request.value
                        << ", snapped to " << actual << ")\n";
                    outputError = true;
                    continue;
                }

                xs.resize(frame.spec.nx);
                ys.resize(frame.spec.nx, actual);
                values.resize(frame.spec.nx);
                for (std::size_t i = 0; i < frame.spec.nx; ++i) {
                    const double x = frame.spec.originX + static_cast<double>(i) * frame.spec.dx;
                    const std::size_t gridIdx = grid.idx(i, static_cast<std::size_t>(idx));
                    xs[i] = x;
                    if (request.quantity == "Bx") {
                        values[i] = grid.Bx[gridIdx];
                    } else if (request.quantity == "By") {
                        values[i] = grid.By[gridIdx];
                    } else if (request.quantity == "Bmag") {
                        values[i] = std::hypot(grid.Bx[gridIdx], grid.By[gridIdx]);
                    } else if (request.quantity == "Hx") {
                        values[i] = grid.Hx[gridIdx];
                    } else if (request.quantity == "Hy") {
                        values[i] = grid.Hy[gridIdx];
                    } else if (request.quantity == "Hmag") {
                        values[i] = std::hypot(grid.Hx[gridIdx], grid.Hy[gridIdx]);
                    } else {
                        values[i] =
                            0.5 * (grid.Bx[gridIdx] * grid.Hx[gridIdx] +
                                   grid.By[gridIdx] * grid.Hy[gridIdx]);
                    }
                }
            }

            if (values.empty()) {
                continue;
            }

            const std::filesystem::path outPath =
                makeFramePath(request.path, timelineActive, frame.index, frameDigits);
            ensureParentDirectory(outPath);
            try {
                motorsim::write_csv_line_profile(outPath.string(), xs, ys, values);
                out << "Frame " << frame.index << ": wrote line_probe '" << request.id
                    << "' to " << outPath << '\n';
            } catch (const std::exception& ex) {
                err << "Frame " << frame.index << ": failed to write line_probe '" << request.id
                    << "': " << ex.what() << '\n';
                outputError = true;
            }
        }

        for (const auto& request : frame.spec.outputs.probes) {
            if (!shouldEmit(request.id)) {
                continue;
            }

            try {
                const motorsim::StressTensorResult probeResult = motorsim::evaluate_maxwell_stress_probe(
                    grid, frame.spec.originX, frame.spec.originY, frame.spec.dx, frame.spec.dy,
                    request.loopXs, request.loopYs);

                result.torqueSamples[request.id] = probeResult;

                const std::filesystem::path outPath =
                    makeFramePath(request.path, timelineActive, frame.index, frameDigits);
                ensureParentDirectory(outPath);

                std::ofstream ofs(outPath);
                if (!ofs.is_open()) {
                    throw std::runtime_error("Failed to open output file");
                }
                ofs << "Fx,Fy,Tz";
                if (result.hasMagneticCoenergy) {
                    ofs << ",CoEnergy";
                }
                ofs << "\n";
                ofs << std::scientific << std::setprecision(12) << probeResult.forceX << "," << probeResult.forceY
                    << "," << probeResult.torqueZ;
                if (result.hasMagneticCoenergy) {
                    ofs << "," << result.magneticCoenergy;
                }
                ofs << "\n";

                const char* quantityLabel = "torque";
                switch (request.quantity) {
                    case motorsim::ScenarioSpec::Outputs::Probe::Quantity::Force:
                        quantityLabel = "force";
                        break;
                    case motorsim::ScenarioSpec::Outputs::Probe::Quantity::ForceAndTorque:
                        quantityLabel = "force_and_torque";
                        break;
                    case motorsim::ScenarioSpec::Outputs::Probe::Quantity::Torque:
                    default:
                        quantityLabel = "torque";
                        break;
                }

                out << "Frame " << frame.index << ": wrote probe '" << request.id << "' (" << quantityLabel
                    << ", stress_tensor) to " << outPath << " [Fx=" << probeResult.forceX
                    << ", Fy=" << probeResult.forceY << ", Tz=" << probeResult.torqueZ;
                if (result.hasMagneticCoenergy) {
                    out << ", Wm=" << result.magneticCoenergy;
                }
                out << "]\n";
            } catch (const std::exception& ex) {
                err << "Frame " << frame.index << ": failed to evaluate probe '" << request.id
                    << "': " << ex.what() << '\n';
                outputError = true;
            }
        }

        if (!frame.spec.outputs.backEmfProbes.empty()) {
            for (const auto& request : frame.spec.outputs.backEmfProbes) {
                if (!shouldEmit(request.id)) {
                    continue;
                }
                if (!request.frameIndices.empty() &&
                    std::find(request.frameIndices.begin(), request.frameIndices.end(), frame.index) ==
                        request.frameIndices.end()) {
                    continue;
                }

                try {
                    double flux = 0.0;
                    if (request.shape == motorsim::ScenarioSpec::Outputs::BackEmfProbe::RegionShape::Polygon) {
                        flux = motorsim::integrate_polygon_flux_component(
                            grid, frame.spec.originX, frame.spec.originY, frame.spec.dx, frame.spec.dy, request.xs,
                            request.ys, request.minX, request.maxX, request.minY, request.maxY, request.component);
                    } else {
                        flux = motorsim::integrate_rect_flux_component(
                            grid, frame.spec.originX, frame.spec.originY, frame.spec.dx, frame.spec.dy, request.minX,
                            request.maxX, request.minY, request.maxY, request.component);
                    }

                    result.backEmfFluxes[request.id] = flux;
                    out << "Frame " << frame.index << ": sampled back_emf '" << request.id << "' flux=" << flux
                        << '\n';
                } catch (const std::exception& ex) {
                    err << "Frame " << frame.index << ": failed to evaluate back_emf '" << request.id
                        << "': " << ex.what() << '\n';
                    outputError = true;
                }
            }
        }
    }

    if (writeMidline) {
        const std::filesystem::path basePath{"outputs/twowire_midline.csv"};
        const std::filesystem::path csvPath =
            makeFramePath(basePath.string(), timelineActive, frame.index, frameDigits);
        ensureParentDirectory(csvPath);

        const std::size_t iMid = frame.spec.nx / 2;
        std::vector<double> xs;
        std::vector<double> ys;
        std::vector<double> bmag;
        xs.reserve(frame.spec.ny);
        ys.reserve(frame.spec.ny);
        bmag.reserve(frame.spec.ny);
        for (std::size_t j = 0; j < frame.spec.ny; ++j) {
            const double x = frame.spec.originX + static_cast<double>(iMid) * frame.spec.dx;
            const double y = frame.spec.originY + static_cast<double>(j) * frame.spec.dy;
            const std::size_t idx = grid.idx(iMid, j);
            xs.push_back(x);
            ys.push_back(y);
            bmag.push_back(std::hypot(grid.Bx[idx], grid.By[idx]));
        }

        try {
            motorsim::write_csv_line_profile(csvPath.string(), xs, ys, bmag);
            out << "Frame " << frame.index << ": wrote midline profile to " << csvPath << '\n';
        } catch (const std::exception& ex) {
            err << "Frame " << frame.index << ": failed to write midline CSV: " << ex.what()
                << '\n';
            outputError = true;
        }
    }

    result.solverSuccess = report.converged;
    result.outputError = outputError;
    result.stdoutLog = out.str();
    result.stderrLog = err.str();
    result.AzSolution = grid.Az;
    progressGuard.markFinal(report.iters, report.relResidual, report.converged);
    return result;
}

}  // namespace

int main(int argc, char** argv) {
    using namespace motorsim;

    std::optional<std::string> scenarioPath;
    bool runSolve = false;
    bool writeMidline = false;
    bool listOutputs = false;
    bool parallelFramesFlag = false;
    std::optional<std::string> outputsFilterArg;
    std::optional<int> maxItersOverride;
    std::optional<double> tolOverride;
    std::optional<double> omegaOverride;
    std::optional<double> harmonicFreqOverride;
    std::optional<double> harmonicOmegaOverride;
    std::optional<std::string> vtkSeriesPath;
    std::optional<std::string> solverOverride;
    std::optional<std::string> preconditionerOverride;
    std::optional<bool> warmStartOverride;
    std::optional<bool> prolongationOverride;
    std::optional<std::size_t> coarseNxOverride;
    std::optional<std::size_t> coarseNyOverride;
    std::optional<double> progressEveryOverride;
    std::optional<std::size_t> snapshotEveryOverride;
    std::optional<std::string> progressHistoryPath;
    std::optional<bool> quietOverride;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--scenario") {
            if (i + 1 >= argc) {
                std::cerr << "--scenario requires a path argument\n";
                printUsage();
                return 1;
            }
            scenarioPath = std::string(argv[++i]);
        } else if (arg == "--solve") {
            runSolve = true;
        } else if (arg == "--write-midline") {
            writeMidline = true;
        } else if (arg == "--list-outputs") {
            listOutputs = true;
        } else if (arg == "--outputs") {
            if (i + 1 >= argc) {
                std::cerr << "--outputs requires a comma-separated list or 'all'/'none'\n";
                printUsage();
                return 1;
            }
            outputsFilterArg = std::string(argv[++i]);
        } else if (arg == "--parallel-frames") {
            parallelFramesFlag = true;
        } else if (arg == "--vtk-series") {
            if (i + 1 >= argc) {
                std::cerr << "--vtk-series requires a path argument\n";
                printUsage();
                return 1;
            }
            vtkSeriesPath = std::string(argv[++i]);
        } else if (arg == "--max-iters") {
            if (i + 1 >= argc) {
                std::cerr << "--max-iters requires an integer argument\n";
                printUsage();
                return 1;
            }
            int value = 0;
            try {
                value = std::stoi(argv[++i]);
            } catch (const std::exception&) {
                std::cerr << "--max-iters requires a valid integer argument\n";
                return 1;
            }
            if (value <= 0) {
                std::cerr << "--max-iters must be positive\n";
                return 1;
            }
            maxItersOverride = value;
        } else if (arg == "--tol") {
            if (i + 1 >= argc) {
                std::cerr << "--tol requires a positive floating-point argument\n";
                printUsage();
                return 1;
            }
            double value = 0.0;
            try {
                value = std::stod(argv[++i]);
            } catch (const std::exception&) {
                std::cerr << "--tol requires a valid floating-point argument\n";
                return 1;
            }
            if (!(value > 0.0)) {
                std::cerr << "--tol must be positive\n";
                return 1;
            }
            tolOverride = value;
        } else if (arg == "--omega") {
            if (i + 1 >= argc) {
                std::cerr << "--omega requires a floating-point argument\n";
                printUsage();
                return 1;
            }
            double value = 0.0;
            try {
                value = std::stod(argv[++i]);
            } catch (const std::exception&) {
                std::cerr << "--omega requires a valid floating-point argument\n";
                return 1;
            }
            if (!(value > 0.0)) {
                std::cerr << "--omega must be positive\n";
                return 1;
            }
            omegaOverride = value;
        } else if (arg == "--harmonic-frequency") {
            if (i + 1 >= argc) {
                std::cerr << "--harmonic-frequency requires a floating-point argument\n";
                printUsage();
                return 1;
            }
            double value = 0.0;
            try {
                value = std::stod(argv[++i]);
            } catch (const std::exception&) {
                std::cerr << "--harmonic-frequency requires a valid floating-point argument\n";
                return 1;
            }
            if (!(value > 0.0)) {
                std::cerr << "--harmonic-frequency must be positive\n";
                return 1;
            }
            harmonicFreqOverride = value;
        } else if (arg == "--harmonic-omega") {
            if (i + 1 >= argc) {
                std::cerr << "--harmonic-omega requires a floating-point argument\n";
                printUsage();
                return 1;
            }
            double value = 0.0;
            try {
                value = std::stod(argv[++i]);
            } catch (const std::exception&) {
                std::cerr << "--harmonic-omega requires a valid floating-point argument\n";
                return 1;
            }
            if (!(value > 0.0)) {
                std::cerr << "--harmonic-omega must be positive\n";
                return 1;
            }
            harmonicOmegaOverride = value;
        } else if (arg == "--solver") {
            if (i + 1 >= argc) {
                std::cerr << "--solver requires an argument (sor, cg, or harmonic)\n";
                printUsage();
                return 1;
            }
            const std::string value = argv[++i];
            if (value != "sor" && value != "cg" && value != "harmonic" && value != "harmonic_cg") {
                std::cerr << "--solver must be 'sor', 'cg', or 'harmonic'\n";
                return 1;
            }
            solverOverride = value;
        } else if (arg == "--pc") {
            if (i + 1 >= argc) {
                std::cerr << "--pc requires an argument (none, jacobi, or ssor)\n";
                printUsage();
                return 1;
            }
            const std::string value = argv[++i];
            if (value != "none" && value != "jacobi" && value != "ssor") {
                std::cerr << "--pc must be 'none', 'jacobi', or 'ssor'\n";
                return 1;
            }
            preconditionerOverride = value;
        } else if (arg == "--warm-start") {
            warmStartOverride = true;
        } else if (arg == "--no-warm-start") {
            warmStartOverride = false;
        } else if (arg == "--use-prolongation") {
            prolongationOverride = true;
        } else if (arg == "--no-prolongation") {
            prolongationOverride = false;
        } else if (arg == "--coarse-nx") {
            if (i + 1 >= argc) {
                std::cerr << "--coarse-nx requires an integer argument\n";
                return 1;
            }
            std::size_t value = 0;
            try {
                value = static_cast<std::size_t>(std::stoul(argv[++i]));
            } catch (const std::exception&) {
                std::cerr << "--coarse-nx requires a valid integer argument\n";
                return 1;
            }
            if (value < 3) {
                std::cerr << "--coarse-nx must be at least 3\n";
                return 1;
            }
            coarseNxOverride = value;
        } else if (arg == "--coarse-ny") {
            if (i + 1 >= argc) {
                std::cerr << "--coarse-ny requires an integer argument\n";
                return 1;
            }
            std::size_t value = 0;
            try {
                value = static_cast<std::size_t>(std::stoul(argv[++i]));
            } catch (const std::exception&) {
                std::cerr << "--coarse-ny requires a valid integer argument\n";
                return 1;
            }
            if (value < 3) {
                std::cerr << "--coarse-ny must be at least 3\n";
                return 1;
            }
            coarseNyOverride = value;
        } else if (arg == "--progress-every") {
            if (i + 1 >= argc) {
                std::cerr << "--progress-every requires a floating-point argument\n";
                return 1;
            }
            double value = 0.0;
            try {
                value = std::stod(argv[++i]);
            } catch (const std::exception&) {
                std::cerr << "--progress-every requires a valid floating-point argument\n";
                return 1;
            }
            if (value < 0.0) {
                std::cerr << "--progress-every must be non-negative\n";
                return 1;
            }
            progressEveryOverride = value;
        } else if (arg == "--snapshot-every") {
            if (i + 1 >= argc) {
                std::cerr << "--snapshot-every requires an integer argument\n";
                return 1;
            }
            std::size_t value = 0;
            try {
                value = static_cast<std::size_t>(std::stoul(argv[++i]));
            } catch (const std::exception&) {
                std::cerr << "--snapshot-every requires a valid integer argument\n";
                return 1;
            }
            snapshotEveryOverride = value;
        } else if (arg == "--progress-history") {
            if (i + 1 >= argc) {
                std::cerr << "--progress-history requires a path argument\n";
                return 1;
            }
            const std::string pathArg = std::string(argv[++i]);
            if (pathArg.empty()) {
                std::cerr << "--progress-history requires a non-empty path\n";
                return 1;
            }
            progressHistoryPath = pathArg;
        } else if (arg == "--quiet") {
            quietOverride = true;
        } else if (arg == "--no-quiet") {
            quietOverride = false;
        } else if (arg == "--help" || arg == "-h") {
            printUsage();
            return 0;
        } else {
            std::cerr << "Unrecognised argument: " << arg << "\n";
            printUsage();
            return 1;
        }
    }

    if (!scenarioPath) {
        printUsage();
        return 0;
    }

    ScenarioSpec spec;
    try {
        spec = loadScenarioFromJson(*scenarioPath);
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load scenario: " << ex.what() << "\n";
        return 1;
    }

    if (vtkSeriesPath && spec.outputs.vtkSeries.empty()) {
        std::cerr << "--vtk-series specified but scenario defines no vtk_field_series outputs\n";
        return 1;
    }

    std::unordered_set<std::string> availableOutputIds;
    availableOutputIds.reserve(spec.outputs.fieldMaps.size() + spec.outputs.lineProbes.size() +
                              spec.outputs.probes.size() + spec.outputs.backEmfProbes.size() +
                              spec.outputs.mechanicalTraces.size());
    for (const auto& request : spec.outputs.fieldMaps) {
        availableOutputIds.insert(request.id);
    }
    for (const auto& request : spec.outputs.vtkSeries) {
        availableOutputIds.insert(request.id);
    }
    for (const auto& request : spec.outputs.lineProbes) {
        availableOutputIds.insert(request.id);
    }
    for (const auto& request : spec.outputs.probes) {
        availableOutputIds.insert(request.id);
    }
    for (const auto& request : spec.outputs.backEmfProbes) {
        availableOutputIds.insert(request.id);
    }
    for (const auto& request : spec.outputs.polylineOutlines) {
        availableOutputIds.insert(request.id);
    }
    for (const auto& request : spec.outputs.boreProbes) {
        availableOutputIds.insert(request.id);
    }
    for (const auto& request : spec.outputs.mechanicalTraces) {
        availableOutputIds.insert(request.id);
    }
    for (const auto& request : spec.outputs.circuitTraces) {
        availableOutputIds.insert(request.id);
    }

    if (listOutputs) {
        std::cout << "Scenario outputs (" << availableOutputIds.size() << "):\n";
        for (const auto& request : spec.outputs.fieldMaps) {
            std::cout << "  - [field_map] id=" << request.id << ", quantity=" << request.quantity
                      << ", format=" << request.format << ", path=" << request.path << "\n";
        }
        for (const auto& request : spec.outputs.lineProbes) {
            std::cout << "  - [line_probe] id=" << request.id << ", axis=" << request.axis
                      << ", value=" << request.value << ", quantity=" << request.quantity
                      << ", format=" << request.format << ", path=" << request.path << "\n";
        }
        for (const auto& request : spec.outputs.probes) {
            const char* quantityLabel = "torque";
            switch (request.quantity) {
                case motorsim::ScenarioSpec::Outputs::Probe::Quantity::Force:
                    quantityLabel = "force";
                    break;
                case motorsim::ScenarioSpec::Outputs::Probe::Quantity::ForceAndTorque:
                    quantityLabel = "force_and_torque";
                    break;
                case motorsim::ScenarioSpec::Outputs::Probe::Quantity::Torque:
                default:
                    quantityLabel = "torque";
                    break;
            }
            std::cout << "  - [probe] id=" << request.id << ", method=stress_tensor, quantity="
                      << quantityLabel << ", vertices=" << request.loopXs.size() << ", path="
                      << request.path << "\n";
        }
        for (const auto& request : spec.outputs.backEmfProbes) {
            std::string componentLabel = "Bmag";
            switch (request.component) {
                case FluxComponent::Bx:
                    componentLabel = "Bx";
                    break;
                case FluxComponent::By:
                    componentLabel = "By";
                    break;
                case FluxComponent::Bmag:
                default:
                    componentLabel = "Bmag";
                    break;
            }
            std::cout << "  - [back_emf_probe] id=" << request.id << ", component=" << componentLabel
                      << ", path=" << request.path;
            if (!request.frameIndices.empty()) {
                std::cout << ", frames=[";
                for (std::size_t idx = 0; idx < request.frameIndices.size(); ++idx) {
                    if (idx > 0) {
                        std::cout << ',';
                    }
                    std::cout << request.frameIndices[idx];
                }
                std::cout << ']';
            }
            std::cout << "\n";
        }
        for (const auto& request : spec.outputs.mechanicalTraces) {
            std::cout << "  - [mechanical_trace] id=" << request.id << ", rotors=";
            for (std::size_t idx = 0; idx < request.rotors.size(); ++idx) {
                if (idx > 0) {
                    std::cout << ',';
                }
                std::cout << request.rotors[idx];
            }
            std::cout << ", path=" << request.path << "\n";
        }
        for (const auto& request : spec.outputs.circuitTraces) {
            std::cout << "  - [circuit_trace] id=" << request.id << ", path=" << request.path;
            if (!request.circuits.empty()) {
                std::cout << ", circuits=[";
                for (std::size_t idx = 0; idx < request.circuits.size(); ++idx) {
                    if (idx > 0) {
                        std::cout << ',';
                    }
                    std::cout << request.circuits[idx];
                }
                std::cout << ']';
            }
            std::cout << "\n";
        }
    }

    std::unordered_set<std::string> requestedOutputIds;
    bool restrictOutputs = false;
    bool skipOutputs = false;
    if (outputsFilterArg) {
        if (outputsFilterArg->empty()) {
            std::cerr << "--outputs argument must not be empty\n";
            return 1;
        }
        if (*outputsFilterArg == "all") {
            restrictOutputs = false;
        } else if (*outputsFilterArg == "none") {
            skipOutputs = true;
        } else {
            if (availableOutputIds.empty()) {
                std::cerr << "Scenario defines no outputs but --outputs was specified\n";
                return 1;
            }
            const auto ids = splitCommaSeparated(*outputsFilterArg);
            if (ids.empty()) {
                std::cerr << "--outputs requires a comma-separated list of ids or 'all'/'none'\n";
                return 1;
            }
            restrictOutputs = true;
            for (const auto& id : ids) {
                if (availableOutputIds.find(id) == availableOutputIds.end()) {
                    std::cerr << "Requested output id not found in scenario: " << id << "\n";
                    return 1;
                }
                requestedOutputIds.insert(id);
            }
        }
    }

    std::vector<ScenarioFrame> frames;
    try {
        frames = expandScenarioTimeline(spec);
    } catch (const std::exception& ex) {
        std::cerr << "Failed to instantiate scenario timeline: " << ex.what() << "\n";
        return 1;
    }

    if (frames.empty()) {
        std::cerr << "Scenario expansion produced no frames\n";
        return 1;
    }

    Grid2D preview(frames.front().spec.nx, frames.front().spec.ny, frames.front().spec.dx,
                   frames.front().spec.dy);
    try {
        rasterizeScenarioToGrid(frames.front().spec, preview);
    } catch (const std::exception& ex) {
        std::cerr << "Rasterisation error: " << ex.what() << "\n";
        return 1;
    }

    if (!runSolve) {
        if (spec.timeline.empty()) {
            std::cout << "Scenario loaded. Use --solve to run the magnetostatic solver."
                      << std::endl;
        } else {
            std::cout << "Scenario loaded with " << frames.size()
                      << " timeline frame(s). Use --solve to run the magnetostatic solver."
                      << std::endl;
        }

        return 0;
    }

    SolveOptions options{};
    options.maxIters = maxItersOverride ? static_cast<std::size_t>(*maxItersOverride) : 20000;
    options.tol = tolOverride.value_or(1e-6);
    options.omega = omegaOverride.value_or(1.7);

    const std::string solverId = solverOverride.value_or(spec.solverSettings.solverId);
    if (!tryParseSolverKind(solverId, options.kind)) {
        std::cerr << "Unknown solver identifier: " << solverId << "\n";
        return 1;
    }

    std::string preconditionerId = "none";
    if (spec.solverSettings.preconditionerSpecified) {
        preconditionerId = spec.solverSettings.preconditionerId;
    }
    if (preconditionerOverride) {
        preconditionerId = *preconditionerOverride;
    }
    if (!tryParsePreconditioner(preconditionerId, options.preconditioner)) {
        std::cerr << "Unknown preconditioner identifier: " << preconditionerId << "\n";
        return 1;
    }

    if (options.kind == motorsim::SolverKind::Harmonic) {
        double omega = 0.0;
        if (harmonicOmegaOverride) {
            omega = *harmonicOmegaOverride;
        } else if (harmonicFreqOverride) {
            constexpr double twoPi = 2.0 * 3.14159265358979323846;
            omega = *harmonicFreqOverride * twoPi;
        } else if (spec.solverSettings.harmonicOmegaSpecified) {
            omega = spec.solverSettings.harmonicOmega;
        } else if (spec.solverSettings.harmonicFrequencySpecified) {
            constexpr double twoPi = 2.0 * 3.14159265358979323846;
            omega = spec.solverSettings.harmonicFrequencyHz * twoPi;
        }
        if (!(omega > 0.0)) {
            std::cerr << "Harmonic solver requires a positive drive frequency."
                      << " Provide --harmonic-frequency/--harmonic-omega or set solver.frequency_hz in the scenario.\n";
            return 1;
        }
        options.harmonicOmega = omega;
    }

    bool warmStart = spec.solverSettings.warmStartSpecified ? spec.solverSettings.warmStart : false;
    if (warmStartOverride) {
        warmStart = *warmStartOverride;
    }
    options.warmStart = warmStart;

    bool useProlongation = spec.solverSettings.prolongationSpecified ? spec.solverSettings.prolongationEnabled : false;
    if (prolongationOverride) {
        useProlongation = *prolongationOverride;
    }
    options.useProlongation = useProlongation;
    if (options.useProlongation) {
        if (coarseNxOverride) {
            options.coarseNx = *coarseNxOverride;
        } else if (spec.solverSettings.prolongationNx) {
            options.coarseNx = *spec.solverSettings.prolongationNx;
        }
        if (coarseNyOverride) {
            options.coarseNy = *coarseNyOverride;
        } else if (spec.solverSettings.prolongationNy) {
            options.coarseNy = *spec.solverSettings.prolongationNy;
        }
    }

    double progressEverySec = spec.solverSettings.progressEverySecSpecified ?
                                 spec.solverSettings.progressEverySec :
                                 2.0;
    if (progressEveryOverride) {
        progressEverySec = *progressEveryOverride;
    }
    options.progressEverySec = progressEverySec;

    std::size_t snapshotEvery = spec.solverSettings.snapshotEveryItersSpecified ?
                                    spec.solverSettings.snapshotEveryIters :
                                    static_cast<std::size_t>(0);
    if (snapshotEveryOverride) {
        snapshotEvery = *snapshotEveryOverride;
    }
    options.snapshotEveryIters = snapshotEvery;

    bool quiet = spec.solverSettings.quietSpecified ? spec.solverSettings.quiet : false;
    if (quietOverride) {
        quiet = *quietOverride;
    }

    OutputFilter filter{};
    filter.restrict = restrictOutputs;
    filter.skip = skipOutputs;
    if (restrictOutputs) {
        filter.requested = &requestedOutputIds;
    }

    const bool timelineActive = !spec.timeline.empty();
    std::size_t frameDigits = 0;
    if (timelineActive) {
        const std::size_t maxIndex = frames.back().index;
        frameDigits = std::max<std::size_t>(3, std::to_string(maxIndex).size());
    }

    motorsim::CircuitSimulator circuitSim;
    circuitSim.initialize(frames);
    const bool circuitsActive = circuitSim.is_active();

    motorsim::MechanicalSimulator mechanicalSim;
    mechanicalSim.initialize(spec, frames);
    const bool mechanicalActive = mechanicalSim.is_active();

    bool runParallel = false;
    std::size_t threadCount = 1;
    if (parallelFramesFlag && frames.size() > 1) {
        const unsigned int hw = std::thread::hardware_concurrency();
        if (hw > 1) {
            threadCount = std::min<std::size_t>(frames.size(), static_cast<std::size_t>(hw - 1));
            if (threadCount == 0) {
                threadCount = 1;
            }
            if (threadCount > 1) {
                runParallel = true;
            }
        }
    }

    if (options.warmStart && timelineActive && runParallel) {
        if (!quiet) {
            std::cout << "Warm start enabled; solving frames sequentially to preserve dependency order.\n";
        }
        runParallel = false;
        threadCount = 1;
    }

    if (circuitsActive && runParallel) {
        if (!quiet) {
            std::cout << "Circuit coupling requires sequential frame solves. Disabling parallel execution.\n";
        }
        runParallel = false;
        threadCount = 1;
    }

    if (mechanicalActive && runParallel) {
        if (!quiet) {
            std::cout << "Mechanical coupling requires sequential frame solves. Disabling parallel execution.\n";
        }
        runParallel = false;
        threadCount = 1;
    }

    std::vector<std::shared_ptr<FrameProgressState>> progressStates;
    progressStates.reserve(frames.size());
    for (const auto& frame : frames) {
        const double frameTime = timelineActive ? frame.time : 0.0;
        progressStates.push_back(
            std::make_shared<FrameProgressState>(frame.index, frameTime, timelineActive, options.tol));
    }
    double printerInterval = options.progressEverySec > 0.0 ? options.progressEverySec : 2.0;
    ProgressPrinter progressPrinter(progressStates, printerInterval, quiet);

    if (runParallel && !quiet) {
        std::cout << "Solving " << frames.size() << " frames with " << threadCount
                  << " worker threads\n";
    }

    std::vector<FrameRunResult> results(frames.size());
    if (runParallel) {
        std::atomic<std::size_t> nextIndex{0};
        std::vector<std::thread> workers;
        workers.reserve(threadCount);
        for (std::size_t t = 0; t < threadCount; ++t) {
            workers.emplace_back([&]() {
                while (true) {
                    const std::size_t idx = nextIndex.fetch_add(1);
                    if (idx >= frames.size()) {
                        break;
                    }
                    results[idx] = solveFrame(frames[idx], options, filter, timelineActive,
                                              frameDigits, writeMidline, progressStates[idx].get(),
                                              nullptr, nullptr, nullptr, 0.0);
                }
            });
        }
        for (auto& worker : workers) {
            worker.join();
        }
    } else {
        std::vector<double> warmAz;
        bool haveWarmAz = false;
        std::vector<double> transientAz;
        bool haveTransientAz = false;
        for (std::size_t idx = 0; idx < frames.size(); ++idx) {
            const std::vector<double>* guessPtr =
                (options.warmStart && timelineActive && haveWarmAz) ? &warmAz : nullptr;
            motorsim::CircuitSimulator* circuitPtr = circuitsActive ? &circuitSim : nullptr;
            if (mechanicalActive) {
                mechanicalSim.apply_state(frames[idx]);
            }
            if (circuitsActive) {
                circuitSim.update_for_frame(frames[idx], mechanicalActive ? &mechanicalSim : nullptr);
            }
            const std::vector<double>* transientPrev = nullptr;
            double frameTransientDt = 0.0;
            if (frames[idx].spec.transient.enabled) {
                const std::size_t cellCount = frames[idx].spec.nx * frames[idx].spec.ny;
                if (!haveTransientAz || transientAz.size() != cellCount) {
                    transientAz.assign(cellCount, 0.0);
                    haveTransientAz = true;
                }
                transientPrev = &transientAz;
                frameTransientDt = frames[idx].spec.transient.dt;
            }
            results[idx] = solveFrame(frames[idx], options, filter, timelineActive, frameDigits,
                                      writeMidline, progressStates[idx].get(), guessPtr, circuitPtr,
                                      transientPrev, frameTransientDt);
            if (mechanicalActive && results[idx].solverSuccess) {
                ScenarioFrame* nextFramePtr = (idx + 1 < frames.size()) ? &frames[idx + 1] : nullptr;
                mechanicalSim.handle_solved_frame(frames[idx], results[idx].torqueSamples, nextFramePtr);
            }
            if (circuitsActive && results[idx].solverSuccess && idx + 1 < frames.size()) {
                circuitSim.prepare_next_frame(idx, frames[idx], &frames[idx + 1]);
            }
            if (results[idx].solverSuccess) {
                if (frames[idx].spec.transient.enabled) {
                    transientAz = results[idx].AzSolution;
                    haveTransientAz = true;
                }
                if (options.warmStart && timelineActive) {
                    warmAz = results[idx].AzSolution;
                    haveWarmAz = true;
                } else if (frames[idx].spec.transient.enabled && !options.warmStart) {
                    warmAz = results[idx].AzSolution;
                    haveWarmAz = true;
                }
            }
        }
    }

    progressPrinter.stop();

    for (const auto& result : results) {
        if (!result.stdoutLog.empty()) {
            std::cout << result.stdoutLog;
        }
    }
    for (const auto& result : results) {
        if (!result.stderrLog.empty()) {
            std::cerr << result.stderrLog;
        }
    }

    bool solverFailure = false;
    bool outputFailure = false;
    for (const auto& result : results) {
        if (!result.solverSuccess) {
            solverFailure = true;
        }
        if (result.outputError) {
            outputFailure = true;
        }
    }

    if (progressHistoryPath) {
        for (std::size_t idx = 0; idx < results.size(); ++idx) {
            const auto& history = results[idx].progressHistory;
            const std::filesystem::path outPath =
                makeFramePath(*progressHistoryPath, timelineActive, frames[idx].index, frameDigits);
            ensureParentDirectory(outPath);
            std::ofstream ofs(outPath);
            if (!ofs.is_open()) {
                std::cerr << "Failed to open progress history path " << outPath << '\n';
                outputFailure = true;
                continue;
            }
            ofs << "iter,rel_residual,elapsed_seconds\n";
            ofs << std::scientific << std::setprecision(12);
            for (const auto& sample : history) {
                ofs << sample.iter << ',' << sample.relResidual << ',' << sample.elapsedSeconds << '\n';
            }
            ofs.close();
            if (!ofs) {
                std::cerr << "Failed while writing progress history to " << outPath << '\n';
                outputFailure = true;
                continue;
            }
            if (!quiet) {
                std::cout << "Progress history frame " << frames[idx].index << " wrote " << history.size()
                          << " sample(s) to " << outPath << '\n';
            }
        }
    }

    const auto shouldEmitId = [&](const std::string& id) {
        if (filter.skip) {
            return false;
        }
        if (!filter.restrict) {
            return true;
        }
        if (!filter.requested) {
            return false;
        }
        return filter.requested->find(id) != filter.requested->end();
    };

    if (!solverFailure && !filter.skip && !frames.front().spec.outputs.backEmfProbes.empty()) {
        const auto& backEmfRequests = frames.front().spec.outputs.backEmfProbes;
        for (const auto& request : backEmfRequests) {
            if (!shouldEmitId(request.id)) {
                continue;
            }

            std::vector<std::size_t> indices = request.frameIndices;
            if (indices.empty()) {
                indices.resize(frames.size());
                std::iota(indices.begin(), indices.end(), 0);
            }
            std::sort(indices.begin(), indices.end());
            indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
            if (indices.size() < 2) {
                std::cerr << "Back-EMF output '" << request.id
                          << "' requires at least two frame indices\n";
                outputFailure = true;
                continue;
            }

            std::vector<double> times;
            std::vector<double> fluxes;
            times.reserve(indices.size());
            fluxes.reserve(indices.size());
            bool missingData = false;
            for (std::size_t frameIdx : indices) {
                if (frameIdx >= frames.size()) {
                    std::cerr << "Back-EMF output '" << request.id << "' references invalid frame " << frameIdx
                              << "\n";
                    outputFailure = true;
                    missingData = true;
                    break;
                }
                const auto fluxIt = results[frameIdx].backEmfFluxes.find(request.id);
                if (fluxIt == results[frameIdx].backEmfFluxes.end()) {
                    std::cerr << "Back-EMF flux missing for frame " << frameIdx << " (output '" << request.id
                              << "')\n";
                    outputFailure = true;
                    missingData = true;
                    break;
                }
                times.push_back(frames[frameIdx].time);
                fluxes.push_back(fluxIt->second);
            }
            if (missingData) {
                continue;
            }

            std::vector<motorsim::BackEmfSample> samples;
            try {
                samples = motorsim::compute_back_emf_series(indices, times, fluxes);
            } catch (const std::exception& ex) {
                std::cerr << "Failed to evaluate back-EMF '" << request.id << "': " << ex.what() << "\n";
                outputFailure = true;
                continue;
            }

            const std::filesystem::path outPath(request.path);
            ensureParentDirectory(outPath);
            std::ofstream ofs(outPath);
            if (!ofs.is_open()) {
                std::cerr << "Failed to open back-EMF output path " << outPath << '\n';
                outputFailure = true;
                continue;
            }
            ofs << "frame_start,frame_end,t_start,t_end,delta_t,flux_start,flux_end,emf\n";
            ofs << std::scientific << std::setprecision(12);
            for (const auto& sample : samples) {
                const double deltaT = sample.timeEnd - sample.timeStart;
                ofs << sample.frameStart << ',' << sample.frameEnd << ',' << sample.timeStart << ','
                    << sample.timeEnd << ',' << deltaT << ',' << sample.fluxStart << ',' << sample.fluxEnd << ','
                    << sample.emf << '\n';
            }
            ofs.close();
            std::cout << "Back-EMF '" << request.id << "' wrote " << samples.size()
                      << " interval(s) to " << outPath << '\n';
        }
    }

    if (!solverFailure && !filter.skip && !frames.front().spec.outputs.boreProbes.empty()) {
        std::unordered_map<std::string, std::vector<FrameRunResult::BoreSample>> boreSeries;
        for (const auto& result : results) {
            for (const auto& sample : result.boreSamples) {
                if (!shouldEmitId(sample.id)) {
                    continue;
                }
                boreSeries[sample.id].push_back(sample);
            }
        }
        for (const auto& request : frames.front().spec.outputs.boreProbes) {
            if (!shouldEmitId(request.id)) {
                continue;
            }
            auto seriesIt = boreSeries.find(request.id);
            if (seriesIt == boreSeries.end() || seriesIt->second.empty()) {
                std::cerr << "Bore probe '" << request.id << "' produced no samples\n";
                outputFailure = true;
                continue;
            }
            auto& samples = seriesIt->second;
            std::sort(samples.begin(), samples.end(),
                      [](const FrameRunResult::BoreSample& a, const FrameRunResult::BoreSample& b) {
                          return a.time < b.time;
                      });
            const std::filesystem::path outPath(request.path);
            ensureParentDirectory(outPath);
            std::ofstream ofs(outPath);
            if (!ofs.is_open()) {
                std::cerr << "Failed to open bore probe output path " << outPath << '\n';
                outputFailure = true;
                continue;
            }
            ofs << "time,bx,by,bmag,angle_deg\n";
            ofs << std::scientific << std::setprecision(12);
            for (const auto& sample : samples) {
                const double angle = std::atan2(sample.by, sample.bx) * 180.0 / 3.14159265358979323846;
                ofs << sample.time << ',' << sample.bx << ',' << sample.by << ',' << sample.magnitude << ','
                    << angle << '\n';
            }
            ofs.close();
            if (!ofs) {
                std::cerr << "Failed while writing bore probe output " << outPath << '\n';
                outputFailure = true;
                continue;
            }
            std::cout << "Bore probe '" << request.id << "' wrote " << samples.size() << " samples to "
                      << outPath << '\n';
        }
    }

    if (!solverFailure && !filter.skip && timelineActive && !frames.front().spec.outputs.probes.empty()) {
        for (const auto& request : frames.front().spec.outputs.probes) {
            if (!shouldEmitId(request.id)) {
                continue;
            }

            struct AggregatedTorqueRow {
                double time;
                std::size_t frameIndex;
                motorsim::StressTensorResult sample;
                bool hasCoEnergy;
                double coEnergy;
            };

            std::vector<AggregatedTorqueRow> rows;
            rows.reserve(frames.size());
            std::vector<std::size_t> missingFrames;
            bool includeCoEnergy = false;

            for (std::size_t idx = 0; idx < frames.size(); ++idx) {
                const auto sampleIt = results[idx].torqueSamples.find(request.id);
                if (sampleIt == results[idx].torqueSamples.end()) {
                    if (results[idx].solverSuccess) {
                        missingFrames.push_back(frames[idx].index);
                    }
                    continue;
                }

                AggregatedTorqueRow row{frames[idx].time,
                                        frames[idx].index,
                                        sampleIt->second,
                                        results[idx].hasMagneticCoenergy,
                                        results[idx].magneticCoenergy};
                includeCoEnergy = includeCoEnergy || row.hasCoEnergy;
                rows.push_back(row);
            }

            if (rows.empty()) {
                std::cerr << "Torque probe '" << request.id
                          << "' produced no samples for aggregated timeline output.\n";
                outputFailure = true;
                continue;
            }

            if (!missingFrames.empty()) {
                std::cerr << "Torque probe '" << request.id << "' missing samples for frame indices";
                for (std::size_t idx = 0; idx < missingFrames.size(); ++idx) {
                    std::cerr << (idx == 0 ? ' ' : ',') << missingFrames[idx];
                }
                std::cerr << "\n";
                outputFailure = true;
            }

            const std::filesystem::path outPath(request.path);
            ensureParentDirectory(outPath);
            std::ofstream ofs(outPath);
            if (!ofs.is_open()) {
                std::cerr << "Failed to open torque probe output path " << outPath << '\n';
                outputFailure = true;
                continue;
            }

            ofs << "time_s,frame_index,Fx,Fy,Tz";
            if (includeCoEnergy) {
                ofs << ",CoEnergy";
            }
            ofs << '\n';
            ofs << std::scientific << std::setprecision(12);
            for (const auto& row : rows) {
                ofs << row.time << ',' << row.frameIndex << ',' << row.sample.forceX << ',' << row.sample.forceY << ','
                    << row.sample.torqueZ;
                if (includeCoEnergy) {
                    ofs << ',';
                    if (row.hasCoEnergy) {
                        ofs << row.coEnergy;
                    }
                }
                ofs << '\n';
            }
            ofs.close();
            if (!ofs) {
                std::cerr << "Failed while writing aggregated torque probe output to " << outPath << '\n';
                outputFailure = true;
                continue;
            }
            if (!quiet) {
                std::cout << "Torque probe '" << request.id << "' wrote " << rows.size()
                          << " aggregated sample(s) to " << outPath << '\n';
            }
        }
    }

    if (!solverFailure && !filter.skip && !frames.front().spec.outputs.mechanicalTraces.empty()) {
        if (!mechanicalActive) {
            for (const auto& request : frames.front().spec.outputs.mechanicalTraces) {
                if (!shouldEmitId(request.id)) {
                    continue;
                }
                std::cerr << "Mechanical trace '" << request.id
                          << "' requested but mechanical simulator is inactive (timeline overrides rotor angles).\n";
            }
            outputFailure = true;
        } else {
            const auto& histories = mechanicalSim.history();
            constexpr double kPi = 3.14159265358979323846;
            constexpr double kTwoPi = 2.0 * kPi;
            for (const auto& request : frames.front().spec.outputs.mechanicalTraces) {
                if (!shouldEmitId(request.id)) {
                    continue;
                }

                std::vector<std::string> rotorNames = request.rotors;
                if (rotorNames.empty()) {
                    rotorNames.reserve(histories.size());
                    for (const auto& entry : histories) {
                        rotorNames.push_back(entry.first);
                    }
                    std::sort(rotorNames.begin(), rotorNames.end());
                    rotorNames.erase(std::unique(rotorNames.begin(), rotorNames.end()), rotorNames.end());
                }

                if (rotorNames.empty()) {
                    std::cerr << "Mechanical trace '" << request.id
                              << "' could not resolve any rotor histories to record.\n";
                    outputFailure = true;
                    continue;
                }

                const std::filesystem::path outPath(request.path);
                ensureParentDirectory(outPath);
                std::ofstream ofs(outPath);
                if (!ofs.is_open()) {
                    std::cerr << "Failed to open mechanical trace output path " << outPath << '\n';
                    outputFailure = true;
                    continue;
                }

                ofs << "time_s,rotor,angle_deg,omega_rad_s,omega_rpm,torque_Nm\n";
                ofs << std::scientific << std::setprecision(12);

                std::size_t totalSamples = 0;
                bool missingRotor = false;

                for (const auto& rotorName : rotorNames) {
                    const auto histIt = histories.find(rotorName);
                    if (histIt == histories.end()) {
                        std::cerr << "Mechanical trace '" << request.id << "' references unknown rotor '"
                                  << rotorName << "'\n";
                        outputFailure = true;
                        missingRotor = true;
                        continue;
                    }
                    for (const auto& sample : histIt->second) {
                        const double angleDeg = sample.angleRad * 180.0 / kPi;
                        const double rpm = sample.omega * (60.0 / kTwoPi);
                        ofs << sample.time << ',' << rotorName << ',' << angleDeg << ',' << sample.omega << ',' << rpm
                            << ',' << sample.torque << '\n';
                        ++totalSamples;
                    }
                }

                ofs.close();
                if (!ofs) {
                    std::cerr << "Failed while writing mechanical trace to " << outPath << '\n';
                    outputFailure = true;
                    continue;
                }

                if (!missingRotor && !quiet) {
                    std::cout << "Mechanical trace '" << request.id << "' wrote " << totalSamples
                              << " sample(s) to " << outPath << '\n';
                }
            }
        }
    }

    if (!filter.skip && circuitsActive && !frames.front().spec.outputs.circuitTraces.empty()) {
        const auto& histories = circuitSim.history();
        const auto& metadata = circuitSim.trace_metadata();
        for (const auto& request : frames.front().spec.outputs.circuitTraces) {
            if (!shouldEmitId(request.id)) {
                continue;
            }
            const std::filesystem::path outPath(request.path);
            ensureParentDirectory(outPath);
            std::ofstream ofs(outPath);
            if (!ofs.is_open()) {
                std::cerr << "Failed to open circuit trace output path " << outPath << '\n';
                outputFailure = true;
                continue;
            }

            ofs << "time_s,circuit,coil,region,turns,orientation,current_A,effective_current_A,ampere_turns,back_emf_V\n";
            ofs << std::scientific << std::setprecision(12);

            std::unordered_set<std::string> allowedCircuits(request.circuits.begin(), request.circuits.end());
            const bool restrictCircuits = !allowedCircuits.empty();

            std::size_t totalSamples = 0;
            for (const auto& pair : histories) {
                const auto metaIt = metadata.find(pair.first);
                if (metaIt == metadata.end()) {
                    continue;
                }
                const auto& meta = metaIt->second;
                if (restrictCircuits && allowedCircuits.find(meta.circuitId) == allowedCircuits.end()) {
                    continue;
                }
                for (const auto& sample : pair.second) {
                    ofs << sample.time << ',' << meta.circuitId << ',' << meta.coilId << ',' << meta.regionId << ','
                        << meta.turns << ',' << sample.orientation << ',' << sample.current << ','
                        << sample.effectiveCurrent << ',' << sample.ampereTurns << ',' << sample.backEmf << '\n';
                    ++totalSamples;
                }
            }

            ofs.close();
            if (!ofs) {
                std::cerr << "Failed while writing circuit trace to " << outPath << '\n';
                outputFailure = true;
                continue;
            }
            if (!quiet) {
                std::cout << "Circuit trace '" << request.id << "' wrote " << totalSamples
                          << " sample(s) to " << outPath << '\n';
            }
        }
    }

    if (!filter.skip && !frames.front().spec.outputs.polylineOutlines.empty()) {
        const auto loops = buildOutlineLoops(frames.front().spec);
        for (const auto& request : frames.front().spec.outputs.polylineOutlines) {
            if (!shouldEmitId(request.id)) {
                continue;
            }
            try {
                const std::filesystem::path outPath(request.path);
                ensureParentDirectory(outPath);
                motorsim::write_vtp_outlines(outPath.string(), loops);
                std::cout << "Polyline outlines '" << request.id << "' wrote to " << outPath << '\n';
            } catch (const std::exception& ex) {
                std::cerr << "Failed to write polyline outlines '" << request.id << "': " << ex.what()
                          << '\n';
                outputFailure = true;
            }
        }
    }

    if (!solverFailure && vtkSeriesPath) {
        std::unordered_map<std::string, std::vector<motorsim::PvdDataSet>> seriesData;
        for (std::size_t idx = 0; idx < results.size(); ++idx) {
            const auto& frameResult = results[idx];
            for (const auto& frameEntry : frameResult.vtkSeriesFrames) {
                if (!shouldEmitId(frameEntry.id)) {
                    continue;
                }
                seriesData[frameEntry.id].push_back(
                    motorsim::PvdDataSet{frameEntry.time, frameEntry.path});
            }
        }
        if (seriesData.empty()) {
            std::cerr << "--vtk-series specified but no VTK series frames were recorded\n";
            outputFailure = true;
        } else if (seriesData.size() > 1) {
            std::cerr << "--vtk-series currently supports a single VTK series output; scenario provided "
                      << seriesData.size() << '\n';
            outputFailure = true;
        } else {
            const auto& pair = *seriesData.begin();
            auto datasets = pair.second;
            std::sort(datasets.begin(), datasets.end(), [](const motorsim::PvdDataSet& a,
                                                           const motorsim::PvdDataSet& b) {
                return a.time < b.time;
            });
            std::filesystem::path pvdPath(*vtkSeriesPath);
            ensureParentDirectory(pvdPath);
            const std::filesystem::path baseDir = pvdPath.parent_path();
            for (auto& dataset : datasets) {
                std::filesystem::path filePath(dataset.file);
                try {
                    dataset.file = std::filesystem::relative(filePath, baseDir).string();
                } catch (const std::exception&) {
                    dataset.file = filePath.string();
                }
            }
            try {
                motorsim::write_pvd_series(pvdPath.string(), datasets);
                std::cout << "VTK series '" << pair.first << "' wrote " << datasets.size()
                          << " timesteps to " << pvdPath << '\n';
            } catch (const std::exception& ex) {
                std::cerr << "Failed to write VTK series index: " << ex.what() << '\n';
                outputFailure = true;
            }
        }
    }

    if (solverFailure) {
        return 1;
    }
    if (outputFailure) {
        return 1;
    }
    return 0;
}
