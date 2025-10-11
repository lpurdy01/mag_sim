#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"
#include "motorsim/solver.hpp"

#include <filesystem>
#include <iostream>
#include <vector>

namespace {

struct RecordingSink : motorsim::ProgressSink {
    explicit RecordingSink(std::filesystem::path dumpDirectory)
        : directory(std::move(dumpDirectory)) {}

    bool onProgress(const motorsim::ProgressSample& sample) override {
        samples.push_back(sample);
        return true;
    }

    std::optional<std::string> requestFieldDump(const motorsim::ProgressSample& sample) override {
        if (sample.iter == 0 || emittedDump) {
            return std::nullopt;
        }
        emittedDump = true;
        std::filesystem::create_directories(directory);
        const std::filesystem::path outPath = directory / "snapshot.csv";
        return outPath.string();
    }

    std::filesystem::path directory;
    std::vector<motorsim::ProgressSample> samples;
    bool emittedDump{false};
};

}  // namespace

int main() {
    namespace fs = std::filesystem;

    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/two_wire_cancel_test.json").lexically_normal();
    motorsim::ScenarioSpec spec = motorsim::loadScenarioFromJson(scenarioPath.string());

    motorsim::Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    motorsim::rasterizeScenarioToGrid(spec, grid);

    motorsim::SolveOptions options{};
    options.kind = motorsim::SolverKind::SOR;
    options.maxIters = 2000;
    options.tol = 1e-12;
    options.progressEverySec = 0.0;
    options.snapshotEveryIters = 250;

    const fs::path dumpDir = fs::path("outputs") / "progress_test";
    fs::remove_all(dumpDir);

    RecordingSink sink(dumpDir);
    motorsim::SolveResult result = motorsim::solveAz(grid, options, nullptr, &sink);

    if (sink.samples.size() < 2) {
        std::cerr << "Expected multiple progress samples, received " << sink.samples.size() << '\n';
        return 1;
    }

    if (!sink.emittedDump) {
        std::cerr << "Progress sink did not request snapshot" << '\n';
        return 1;
    }

    const fs::path snapshotPath = dumpDir / "snapshot.csv";
    if (!fs::exists(snapshotPath)) {
        std::cerr << "Snapshot file was not created" << '\n';
        return 1;
    }

    if (!result.converged && result.iters < options.maxIters) {
        std::cerr << "Solver terminated early without convergence" << '\n';
        return 1;
    }

    fs::remove(snapshotPath);
    fs::remove(dumpDir);

    return 0;
}

