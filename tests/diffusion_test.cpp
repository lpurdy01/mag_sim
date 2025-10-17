#include "motorsim/ingest.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

int main() {
    using namespace motorsim;
    namespace fs = std::filesystem;

    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/diffusion_test.json").lexically_normal();

    ScenarioSpec spec;
    try {
        spec = loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load diffusion scenario: " << ex.what() << '\n';
        return 1;
    }

    std::vector<ScenarioFrame> frames = expandScenarioTimeline(spec);
    if (frames.empty()) {
        std::cerr << "Diffusion scenario produced no timeline frames\n";
        return 1;
    }

    const double dt = spec.transient.dt > 0.0 ? spec.transient.dt : (frames.size() > 1 ? frames[1].time - frames[0].time : 0.0);
    if (!(dt > 0.0)) {
        std::cerr << "Transient timestep must be positive for diffusion test\n";
        return 1;
    }

    SolveOptions options{};
    options.kind = SolverKind::CG;
    options.maxIters = 25000;
    options.tol = 5e-6;
    options.preconditioner = PreconditionerKind::Jacobi;

    std::vector<double> transientState;
    bool haveTransientState = false;
    Grid2D grid(frames.front().spec.nx, frames.front().spec.ny, frames.front().spec.dx, frames.front().spec.dy);

    for (std::size_t idx = 0; idx < frames.size(); ++idx) {
        ScenarioFrame& frame = frames[idx];
        if (grid.nx != frame.spec.nx || grid.ny != frame.spec.ny) {
            grid.resize(frame.spec.nx, frame.spec.ny, frame.spec.dx, frame.spec.dy);
        }
        try {
            rasterizeScenarioToGrid(frame.spec, grid);
        } catch (const std::exception& ex) {
            std::cerr << "Frame " << idx << ": rasterisation error: " << ex.what() << '\n';
            return 1;
        }

        if (!haveTransientState || transientState.size() != grid.Az.size()) {
            transientState.assign(grid.Az.size(), 0.0);
            haveTransientState = true;
        }

        const SolveResult report = solveTransientStep(grid, options, dt, transientState, nullptr);
        if (!report.converged) {
            std::cerr << "Frame " << idx << ": transient solve did not converge (relResidual=" << report.relResidual
                      << ")\n";
            return 1;
        }

        computeB(grid);

        transientState = grid.Az;
    }

    const double finalTime = frames.back().time;
    if (!(finalTime > 0.0)) {
        std::cerr << "Final frame time must be positive\n";
        return 1;
    }

    const double mu = motorsim::MU0;
    const double sigma = 5.96e7;
    const double alpha = 1.0 / (mu * sigma);

    std::vector<double> depths;
    std::vector<double> measuredBx;

    const double sampleMinY = spec.originY;
    const double sampleMaxY = spec.originY + spec.Ly;

    double bSurface = 0.0;
    std::size_t surfaceCount = 0;

    for (std::size_t i = 0; i < grid.nx; ++i) {
        const double x = spec.originX + static_cast<double>(i) * spec.dx;
        if (x < 0.0) {
            continue;
        }
        double accum = 0.0;
        std::size_t count = 0;
        for (std::size_t j = 0; j < grid.ny; ++j) {
            const double y = spec.originY + static_cast<double>(j) * spec.dy;
            if (y < sampleMinY || y > sampleMaxY) {
                continue;
            }
            const std::size_t idx = grid.idx(i, j);
            accum += grid.Bx[idx];
            ++count;
        }
        if (count == 0) {
            continue;
        }
        const double avg = accum / static_cast<double>(count);
        const double depth = x;
        if (depth < grid.dx * 1.5) {
            bSurface += avg;
            ++surfaceCount;
        }
        depths.push_back(depth);
        measuredBx.push_back(avg);
    }

    if (surfaceCount == 0) {
        std::cerr << "Failed to recover surface field sample\n";
        return 1;
    }

    const double b0 = bSurface / static_cast<double>(surfaceCount);
    if (!(std::abs(b0) > 0.0)) {
        std::cerr << "Surface field magnitude is zero\n";
        return 1;
    }

    double maxRelError = 0.0;
    const double sqrtAlphaT = std::sqrt(alpha * finalTime);
    for (std::size_t idx = 0; idx < depths.size(); ++idx) {
        const double depth = depths[idx];
        const double expected = b0 * std::erfc(depth / (2.0 * sqrtAlphaT));
        const double measured = measuredBx[idx];
        const double denom = std::max(std::abs(expected), 1e-12);
        const double relError = std::abs(measured - expected) / denom;
        if (relError > maxRelError) {
            maxRelError = relError;
        }
    }

    if (!(maxRelError < 0.2)) {
        std::cerr << "Diffusion profile relative error too high: " << maxRelError << '\n';
        return 1;
    }

    return 0;
}
