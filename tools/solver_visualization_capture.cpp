// filename: solver_visualization_capture.cpp
// part of 2D Electromagnetic Motor Simulator
// MIT License

#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"
#include "motorsim/io_csv.hpp"
#include "motorsim/solver.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace
{

namespace fs = std::filesystem;

struct CsvProgressSink : motorsim::ProgressSink
{
    explicit CsvProgressSink(fs::path historyPath, fs::path snapshotDir = {},
                             std::size_t maxSnapshots = 0)
        : historyPath(std::move(historyPath)), snapshotDir(std::move(snapshotDir)),
          maxSnapshots(maxSnapshots)
    {
    }

    bool onProgress(const motorsim::ProgressSample &sample) override
    {
        samples.push_back(sample);
        return true;
    }

    std::optional<std::string> requestFieldDump(const motorsim::ProgressSample &sample) override
    {
        if (snapshotDir.empty())
        {
            return std::nullopt;
        }
        if (maxSnapshots > 0 && snapshotsWritten >= maxSnapshots)
        {
            return std::nullopt;
        }

        // Variable snapshot interval: every 5 iterations before 150, every 20 after
        const std::size_t desiredInterval = (sample.iter < 150) ? 5 : 20;
        if (sample.iter > 0 && sample.iter % desiredInterval != 0)
        {
            return std::nullopt;
        }

        if (!fs::exists(snapshotDir))
        {
            std::error_code ec;
            fs::create_directories(snapshotDir, ec);
        }
        ++snapshotsWritten;
        std::ostringstream name;
        name << "snapshot_iter_" << std::setw(6) << std::setfill('0') << sample.iter << ".csv";
        snapshotIterations.push_back(sample.iter);
        return (snapshotDir / name.str()).string();
    }

    bool flush() const
    {
        if (samples.empty())
        {
            return false;
        }
        std::error_code ec;
        fs::create_directories(historyPath.parent_path(), ec);
        std::ofstream ofs(historyPath);
        if (!ofs.is_open())
        {
            std::cerr << "Failed to open history CSV: " << historyPath << '\n';
            return false;
        }
        ofs << "iter,rel_residual,elapsed_seconds\n";
        for (const auto &sample : samples)
        {
            ofs << sample.iter << ',' << sample.relResidual << ',' << sample.elapsedSeconds << '\n';
        }
        return true;
    }

    fs::path historyPath{};
    fs::path snapshotDir{};
    std::size_t maxSnapshots{0};
    mutable std::vector<motorsim::ProgressSample> samples;
    mutable std::size_t snapshotsWritten{0};
    mutable std::vector<std::size_t> snapshotIterations;
};

struct FrameSolveSummary
{
    motorsim::SolveResult result{};
    std::vector<double> az;
};

FrameSolveSummary solveFrame(const motorsim::ScenarioFrame &frame,
                             const motorsim::SolveOptions &options,
                             const motorsim::InitialGuess *guess, CsvProgressSink &sink)
{
    motorsim::Grid2D grid(frame.spec.nx, frame.spec.ny, frame.spec.dx, frame.spec.dy);
    motorsim::rasterizeScenarioToGrid(frame.spec, grid);
    FrameSolveSummary summary{};
    summary.result = motorsim::solveAz(grid, options, guess, &sink);
    summary.az = grid.Az;
    return summary;
}

bool writeWarmStartRuns(const fs::path &outputRoot)
{
    const fs::path scenarioPath =
        (fs::path("inputs") / "tests" / "magnet_strip_test.json").lexically_normal();
    motorsim::ScenarioSpec baseSpec;
    try
    {
        baseSpec = motorsim::loadScenarioFromJson(scenarioPath.string());
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Failed to load scenario for warm-start capture: " << ex.what() << '\n';
        return false;
    }

    motorsim::ScenarioFrame baseFrame{};
    baseFrame.index = 0;
    baseFrame.time = 0.0;
    baseFrame.spec = baseSpec;

    auto makeMagnetFrame = [&](std::size_t index, double time, double scaleMy, double deltaMx)
    {
        motorsim::ScenarioFrame frame = baseFrame;
        frame.index = index;
        frame.time = time;
        for (auto &magnet : frame.spec.magnetRegions)
        {
            magnet.My *= scaleMy;
            magnet.Mx += deltaMx;
        }
        return frame;
    };

    std::vector<motorsim::ScenarioFrame> frames;
    frames.push_back(baseFrame);
    frames.push_back(makeMagnetFrame(1, 1e-3, 0.97, 4.0e4));
    frames.push_back(makeMagnetFrame(2, 2e-3, 0.94, 8.0e4));

    motorsim::SolveOptions options{};
    options.kind = motorsim::SolverKind::CG;
    options.maxIters = 20000;
    options.tol = 1e-6;
    options.progressEverySec = 0.0;

    const fs::path warmDir = outputRoot / "warm_start";
    const fs::path coldDir = warmDir / "cold";
    const fs::path seededDir = warmDir / "warm";

    std::error_code ec;
    fs::create_directories(coldDir, ec);
    fs::create_directories(seededDir, ec);

    for (const auto &frame : frames)
    {
        const fs::path historyPath = coldDir / ("frame_" + std::to_string(frame.index) + ".csv");
        CsvProgressSink sink(historyPath);
        FrameSolveSummary summary = solveFrame(frame, options, nullptr, sink);
        if (!summary.result.converged)
        {
            std::cerr << "Cold solve failed for frame " << frame.index << '\n';
            return false;
        }
        sink.flush();
    }

    options.warmStart = true;
    std::vector<double> previousAz;

    for (std::size_t idx = 0; idx < frames.size(); ++idx)
    {
        const auto &frame = frames[idx];
        const fs::path historyPath = seededDir / ("frame_" + std::to_string(frame.index) + ".csv");
        CsvProgressSink sink(historyPath);
        motorsim::InitialGuess guess{};
        const motorsim::InitialGuess *guessPtr = nullptr;
        if (idx > 0)
        {
            guess.Az0 = &previousAz;
            guessPtr = &guess;
        }
        FrameSolveSummary summary = solveFrame(frame, options, guessPtr, sink);
        if (!summary.result.converged)
        {
            std::cerr << "Warm solve failed for frame " << frame.index << '\n';
            return false;
        }
        sink.flush();
        previousAz = summary.az;
    }

    return true;
}

bool writeSnapshotRun(const fs::path &outputRoot, motorsim::SolverKind solverKind,
                      std::string_view label)
{
    const fs::path scenarioPath =
        (fs::path("inputs") / "iron_ring_magnet_viz.json").lexically_normal();
    motorsim::ScenarioSpec spec;
    try
    {
        spec = motorsim::loadScenarioFromJson(scenarioPath.string());
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Failed to load snapshot scenario: " << ex.what() << '\n';
        return false;
    }

    motorsim::Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    motorsim::rasterizeScenarioToGrid(spec, grid);

    motorsim::SolveOptions options{};
    options.kind = solverKind;
    options.maxIters = 8000;
    options.tol = 1e-8;
    options.progressEverySec = 0.0;
    options.snapshotEveryIters = 1; // Check every iteration (filtering done in sink)
    if (solverKind == motorsim::SolverKind::SOR)
    {
        options.omega = 1.8;
        options.maxIters = 60000;
        options.tol = 5e-6;
    }

    fs::path snapshotDir = outputRoot / "progress_snapshots" / std::string(label) / "snapshots";
    fs::path historyPath = outputRoot / "progress_snapshots" / std::string(label) / "history.csv";

    CsvProgressSink sink(historyPath, snapshotDir, 150); // Increased max snapshots
    motorsim::SolveResult result = motorsim::solveAz(grid, options, nullptr, &sink);
    if (!result.converged)
    {
        std::cerr << "Snapshot solve (" << label
                  << ") did not converge (residual=" << result.relResidual << ")\n";
        return false;
    }
    sink.flush();
    return true;
}

} // namespace

int main()
{
    try
    {
        const fs::path outputRoot = fs::path("outputs") / "visualization";
        std::error_code ec;
        fs::remove_all(outputRoot, ec);
        fs::create_directories(outputRoot, ec);

        bool ok = writeWarmStartRuns(outputRoot);
        ok = writeSnapshotRun(outputRoot, motorsim::SolverKind::CG, "cg") && ok;
        ok = writeSnapshotRun(outputRoot, motorsim::SolverKind::SOR, "sor") && ok;

        if (!ok)
        {
            std::cerr << "Failed to generate solver visualization datasets." << '\n';
            return 1;
        }

        std::cout << "Solver visualization datasets written to " << outputRoot << '\n';
        return 0;
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Unhandled exception: " << ex.what() << '\n';
        return 1;
    }
}
