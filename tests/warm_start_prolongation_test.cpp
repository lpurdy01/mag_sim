#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"
#include "motorsim/solver.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <vector>

namespace {

using motorsim::InitialGuess;
using motorsim::ScenarioFrame;
using motorsim::ScenarioSpec;
using motorsim::SolveOptions;
using motorsim::SolveResult;

struct FrameSolveSummary {
    SolveResult result{};
    std::vector<double> az;
};

FrameSolveSummary solve_frame(const ScenarioFrame& frame,
                              SolveOptions options,
                              const InitialGuess* guess = nullptr) {
    motorsim::Grid2D grid(frame.spec.nx, frame.spec.ny, frame.spec.dx, frame.spec.dy);
    motorsim::rasterizeScenarioToGrid(frame.spec, grid);
    FrameSolveSummary summary{};
    summary.result = motorsim::solveAz(grid, options, guess);
    summary.az = grid.Az;
    return summary;
}

SolveResult solve_single_spec(const ScenarioSpec& spec, SolveOptions options) {
    motorsim::Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    motorsim::rasterizeScenarioToGrid(spec, grid);
    return motorsim::solveAz(grid, options);
}

double max_abs_diff(const std::vector<double>& a, const std::vector<double>& b) {
    const std::size_t n = std::min(a.size(), b.size());
    double maxDiff = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        maxDiff = std::max(maxDiff, std::abs(a[i] - b[i]));
    }
    return maxDiff;
}

}  // namespace

int main() {
    namespace fs = std::filesystem;

    const fs::path magnetPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/magnet_strip_test.json").lexically_normal();
    ScenarioSpec warmStartSpec = motorsim::loadScenarioFromJson(magnetPath.string());

    ScenarioFrame baseFrame{};
    baseFrame.index = 0;
    baseFrame.time = 0.0;
    baseFrame.spec = warmStartSpec;

    auto makeMagnetFrame = [&](std::size_t index, double time, double scaleMy, double deltaMx) {
        ScenarioFrame frame = baseFrame;
        frame.index = index;
        frame.time = time;
        for (auto& magnet : frame.spec.magnetRegions) {
            magnet.My *= scaleMy;
            magnet.Mx += deltaMx;
        }
        return frame;
    };

    std::vector<ScenarioFrame> frames;
    frames.push_back(baseFrame);
    frames.push_back(makeMagnetFrame(1, 1e-3, 0.97, 4.0e4));
    frames.push_back(makeMagnetFrame(2, 2e-3, 0.94, 8.0e4));

    SolveOptions cgOptions{};
    cgOptions.kind = motorsim::SolverKind::CG;
    cgOptions.maxIters = 20000;
    cgOptions.tol = 1e-6;
    cgOptions.warmStart = false;

    std::vector<FrameSolveSummary> coldSummaries;
    coldSummaries.reserve(frames.size());
    for (const ScenarioFrame& frame : frames) {
        FrameSolveSummary summary = solve_frame(frame, cgOptions);
        if (!summary.result.converged || summary.result.relResidual >= cgOptions.tol) {
            std::cerr << "Cold solve failed on frame " << frame.index << '\n';
            return 1;
        }
        coldSummaries.push_back(std::move(summary));
    }

    cgOptions.warmStart = true;
    std::vector<FrameSolveSummary> warmSummaries;
    warmSummaries.reserve(frames.size());
    std::vector<double> previousAz;
    for (std::size_t idx = 0; idx < frames.size(); ++idx) {
        const ScenarioFrame& frame = frames[idx];
        InitialGuess guess{};
        if (idx > 0) {
            guess.Az0 = &previousAz;
        }
        FrameSolveSummary summary = solve_frame(frame, cgOptions, (idx > 0) ? &guess : nullptr);
        if (!summary.result.converged || summary.result.relResidual >= cgOptions.tol) {
            std::cerr << "Warm solve failed on frame " << frame.index << '\n';
            return 1;
        }
        if (summary.az.size() != frame.spec.nx * frame.spec.ny) {
            std::cerr << "Warm solve returned unexpected grid size" << '\n';
            return 1;
        }
        previousAz = summary.az;
        warmSummaries.push_back(std::move(summary));
    }

    std::size_t coldIterSum = 0;
    std::size_t warmIterSum = 0;
    for (std::size_t idx = 1; idx < frames.size(); ++idx) {
        coldIterSum += coldSummaries[idx].result.iters;
        warmIterSum += warmSummaries[idx].result.iters;
        const double diff = max_abs_diff(coldSummaries[idx].az, warmSummaries[idx].az);
        if (diff > 1e-4) {
            std::cerr << "Warm start solution deviates from cold solve on frame " << frames[idx].index
                      << " (maxDiff=" << diff << ")" << '\n';
            return 1;
        }
    }

    if (!(warmIterSum < coldIterSum && warmIterSum * 100 < coldIterSum * 85)) {
        std::cerr << "Warm start did not reduce iteration count sufficiently (cold=" << coldIterSum
                  << ", warm=" << warmIterSum << ")" << '\n';
        return 1;
    }

    const fs::path prolongationPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/magnet_strip_test.json").lexically_normal();
    ScenarioSpec prolongationSpec = motorsim::loadScenarioFromJson(prolongationPath.string());

    SolveOptions fineOptions{};
    fineOptions.kind = motorsim::SolverKind::CG;
    fineOptions.maxIters = 30000;
    fineOptions.tol = 5e-7;
    fineOptions.useProlongation = false;

    SolveResult coldResult = solve_single_spec(prolongationSpec, fineOptions);
    if (!coldResult.converged || coldResult.relResidual >= fineOptions.tol) {
        std::cerr << "Fine-grid baseline solve failed" << '\n';
        return 1;
    }

    fineOptions.useProlongation = true;
    fineOptions.coarseNx = prolongationSpec.nx / 2;
    fineOptions.coarseNy = prolongationSpec.ny / 2;
    SolveResult prolongatedResult = solve_single_spec(prolongationSpec, fineOptions);
    if (!prolongatedResult.converged || prolongatedResult.relResidual >= fineOptions.tol) {
        std::cerr << "Fine-grid prolongated solve failed" << '\n';
        return 1;
    }

    if (!(prolongatedResult.iters < coldResult.iters && prolongatedResult.iters * 100 < coldResult.iters * 95)) {
        std::cerr << "Prolongation did not sufficiently reduce iterations (baseline=" << coldResult.iters
                  << ", prolong=" << prolongatedResult.iters << ")" << '\n';
        return 1;
    }

    return 0;
}

