#include "motorsim/ingest.hpp"
#include "motorsim/io_csv.hpp"
#include "motorsim/io_vtk.hpp"
#include "motorsim/probes.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
#include <sstream>
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
                 " [--omega VALUE] [--parallel-frames] [--vtk-series PATH]\n";
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
};

class FrameProgressState {
public:
    FrameProgressState(std::size_t index, double timeValue, bool hasTimeValue)
        : frameIndex(index), frameTime(timeValue), hasTime(hasTimeValue), lastUpdate(Clock::now()) {}

    void markStarted() {
        std::lock_guard<std::mutex> lock(mutex);
        active = true;
        completed = false;
        success = false;
        lastIteration = 0;
        lastResidual = 0.0;
        lastUpdate = Clock::now();
    }

    void update(std::size_t iteration, double residual) {
        std::lock_guard<std::mutex> lock(mutex);
        active = true;
        lastIteration = iteration;
        lastResidual = residual;
        lastUpdate = Clock::now();
    }

    void markFinished(std::size_t iteration, double residual, bool didSucceed) {
        std::lock_guard<std::mutex> lock(mutex);
        active = false;
        completed = true;
        success = didSucceed;
        lastIteration = iteration;
        lastResidual = residual;
        lastUpdate = Clock::now();
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
        return snap;
    }

    const std::size_t frameIndex;
    const double frameTime;
    const bool hasTime;

private:
    using Clock = std::chrono::steady_clock;

    mutable std::mutex mutex;
    bool active{false};
    bool completed{false};
    bool success{false};
    std::size_t lastIteration{0};
    double lastResidual{0.0};
    Clock::time_point lastUpdate;
};

class ProgressPrinter {
public:
    ProgressPrinter(std::vector<std::shared_ptr<FrameProgressState>> states,
                    std::chrono::steady_clock::duration interval)
        : states_(std::move(states)), interval_(interval) {
        if (!states_.empty() && interval_ > std::chrono::steady_clock::duration::zero()) {
            worker_ = std::thread([this]() { run(); });
        }
    }

    ~ProgressPrinter() { stop(); }

    void stop() {
        const bool wasRunning = !stopped_.exchange(true);
        if (wasRunning && worker_.joinable()) {
            worker_.join();
        }
        if (!states_.empty() && !finalPrinted_.exchange(true)) {
            printStatus(true);
        }
    }

private:
    void run() {
        while (!stopped_.load()) {
            std::this_thread::sleep_for(interval_);
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

    bool printStatus(bool force) {
        if (states_.empty()) {
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

        std::ostringstream oss;
        oss << "[progress] completed " << completed << "/" << total;
        if (failed > 0) {
            oss << " (failed: " << failed << ")";
        }
        oss << " frames; active: " << active << "; queued: " << queued << '\n';

        const auto now = std::chrono::steady_clock::now();
        for (const auto& snap : snapshots) {
            if (!snap.active) {
                continue;
            }
            oss << "    Frame " << snap.frameIndex;
            if (snap.hasTime) {
                oss << " (t=" << snap.frameTime << " s)";
            }
            oss << " iter=" << snap.lastIteration;
            if (snap.lastIteration > 0) {
                std::ostringstream residualStream;
                residualStream.setf(std::ios::scientific, std::ios::floatfield);
                residualStream << std::setprecision(3) << snap.lastResidual;
                oss << " relResidual=" << residualStream.str();
            }
            const auto ageSeconds =
                std::chrono::duration_cast<std::chrono::seconds>(now - snap.lastUpdate).count();
            if (ageSeconds > 0) {
                oss << " (updated " << ageSeconds << "s ago)";
            }
            oss << '\n';
        }

        std::cout << oss.str();
        std::cout.flush();

        return completed == total;
    }

    std::vector<std::shared_ptr<FrameProgressState>> states_;
    std::chrono::steady_clock::duration interval_;
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
    stem << base.stem().string() << "_frame_";
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

FrameRunResult solveFrame(const motorsim::ScenarioFrame& frame, const motorsim::SolveOptions& options,
                          const OutputFilter& filter, bool timelineActive, std::size_t frameDigits,
                          bool writeMidline, FrameProgressState* progressState) {
    FrameRunResult result{};
    std::ostringstream out;
    std::ostringstream err;
    FrameProgressGuard progressGuard(progressState);

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
    if (progressState) {
        solveOptions.progressCallback = [progressState](const motorsim::SolverProgress& progress) {
            progressState->update(progress.iteration, progress.relResidual);
        };
    }

    const motorsim::SolveReport report = motorsim::solveAz_GS_SOR(grid, solveOptions);
    if (!report.converged) {
        err << "Frame " << frame.index
            << ": solver did not converge. relResidual=" << report.relResidual << '\n';
        progressGuard.markFinal(report.iters, report.relResidual, false);
        result.stderrLog = err.str();
        return result;
    }

    motorsim::computeB(grid);

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

    const bool needHField = fieldMapNeedsH || fieldMapNeedsEnergy || lineProbeNeedsH ||
                            lineProbeNeedsEnergy || vtkSeriesNeedsH || vtkSeriesNeedsEnergy;
    if (needHField) {
        motorsim::computeH(grid);
    }

    out << "Frame " << frame.index;
    if (timelineActive) {
        out << " (t=" << frame.time << " s)";
    }
    out << ": solved in " << report.iters << " iterations, relResidual=" << report.relResidual
        << '\n';

    bool outputError = false;

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
                const motorsim::StressTensorResult result = motorsim::evaluate_maxwell_stress_probe(
                    grid, frame.spec.originX, frame.spec.originY, frame.spec.dx, frame.spec.dy,
                    request.loopXs, request.loopYs);

                const std::filesystem::path outPath =
                    makeFramePath(request.path, timelineActive, frame.index, frameDigits);
                ensureParentDirectory(outPath);

                std::ofstream ofs(outPath);
                if (!ofs.is_open()) {
                    throw std::runtime_error("Failed to open output file");
                }
                ofs << "Fx,Fy,Tz\n";
                ofs << std::scientific << std::setprecision(12) << result.forceX << "," << result.forceY
                    << "," << result.torqueZ << "\n";

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
                    << ", stress_tensor) to " << outPath << " [Fx=" << result.forceX
                    << ", Fy=" << result.forceY << ", Tz=" << result.torqueZ << "]\n";
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

    result.solverSuccess = true;
    result.outputError = outputError;
    result.stdoutLog = out.str();
    result.stderrLog = err.str();
    progressGuard.markFinal(report.iters, report.relResidual, true);
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
    std::optional<std::string> vtkSeriesPath;

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
                              spec.outputs.probes.size() + spec.outputs.backEmfProbes.size());
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
    options.maxIters = maxItersOverride.value_or(20000);
    options.tol = tolOverride.value_or(1e-6);
    options.omega = omegaOverride.value_or(1.7);
    options.progressInterval = std::chrono::seconds(2);

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

    std::vector<std::shared_ptr<FrameProgressState>> progressStates;
    progressStates.reserve(frames.size());
    for (const auto& frame : frames) {
        const double frameTime = timelineActive ? frame.time : 0.0;
        progressStates.push_back(
            std::make_shared<FrameProgressState>(frame.index, frameTime, timelineActive));
    }
    auto printerInterval = options.progressInterval;
    if (printerInterval <= std::chrono::steady_clock::duration::zero()) {
        printerInterval = std::chrono::seconds(2);
    }
    ProgressPrinter progressPrinter(progressStates, printerInterval);

    if (runParallel) {
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
                                              frameDigits, writeMidline, progressStates[idx].get());
                }
            });
        }
        for (auto& worker : workers) {
            worker.join();
        }
    } else {
        for (std::size_t idx = 0; idx < frames.size(); ++idx) {
            results[idx] = solveFrame(frames[idx], options, filter, timelineActive, frameDigits,
                                      writeMidline, progressStates[idx].get());
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
