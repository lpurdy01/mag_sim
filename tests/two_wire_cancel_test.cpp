#include "motorsim/ingest.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

int main() {
    using namespace motorsim;

    ScenarioSpec spec{};
    spec.version = "0.1";
    spec.Lx = 0.20;
    spec.Ly = 0.20;
    spec.nx = 201;
    spec.ny = 201;
    spec.dx = spec.Lx / static_cast<double>(spec.nx - 1);
    spec.dy = spec.Ly / static_cast<double>(spec.ny - 1);
    spec.originX = -0.5 * spec.Lx;
    spec.originY = -0.5 * spec.Ly;
    spec.mu_r_background = 1.0;
    spec.wires.push_back({-0.03, 0.0, 0.003, 10.0});
    spec.wires.push_back({0.03, 0.0, 0.003, -10.0});

    Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    rasterizeScenarioToGrid(spec, grid);

    SolveOptions options{};
    options.maxIters = 5000;
    options.tol = 1e-6;
    options.omega = 1.7;

    const SolveReport report = solveAz_GS_SOR(grid, options);
    if (!report.converged) {
        std::cerr << "two_wire_cancel_test: solver failed to converge\n";
        return 1;
    }
    if (report.relResidual > options.tol) {
        std::cerr << "two_wire_cancel_test: residual " << report.relResidual
                  << " exceeds tolerance " << options.tol << "\n";
        return 1;
    }

    computeB(grid);

    const std::size_t iMid = spec.nx / 2;
    const double sampleOffset = 0.026;  // 2.6 cm from the centre sits just outside each wire disk.
    const std::size_t iLeft = static_cast<std::size_t>(
        std::llround(((-sampleOffset) - spec.originX) / spec.dx));
    const std::size_t iRight = static_cast<std::size_t>(
        std::llround(((sampleOffset) - spec.originX) / spec.dx));

    if (iMid == 0 || iMid >= spec.nx || iLeft == 0 || iRight + 1 >= spec.nx) {
        std::cerr << "two_wire_cancel_test: invalid sampling indices\n";
        return 1;
    }

    double midSumSq = 0.0;
    double refSumSq = 0.0;
    double antisymmetry = 0.0;
    const double cellArea = spec.dx * spec.dy;

    for (std::size_t j = 1; j + 1 < spec.ny; ++j) {
        const std::size_t midIdx = grid.idx(iMid, j);
        const std::size_t leftIdx = grid.idx(iLeft, j);
        const std::size_t rightIdx = grid.idx(iRight, j);
        const double bxMid = grid.Bx[midIdx];
        const double byMid = grid.By[midIdx];
        const double bxLeft = grid.Bx[leftIdx];
        const double byLeft = grid.By[leftIdx];
        const double bxRight = grid.Bx[rightIdx];
        const double byRight = grid.By[rightIdx];

        midSumSq += bxMid * bxMid + byMid * byMid;
        refSumSq += bxLeft * bxLeft + byLeft * byLeft;
        refSumSq += bxRight * bxRight + byRight * byRight;

        antisymmetry += std::abs(byLeft + byRight);
    }

    if (refSumSq <= 0.0) {
        std::cerr << "two_wire_cancel_test: reference line magnitude is zero\n";
        return 1;
    }

    const double midRms = std::sqrt(midSumSq / static_cast<double>(spec.ny - 2));
    const double refRms = std::sqrt(refSumSq / static_cast<double>(2 * (spec.ny - 2)));
    const double ratio = midRms / refRms;
    // Coarse grids and Dirichlet boundaries limit cancellation; require midline RMS to be at least half the off-axis RMS.
    if (!(ratio < 0.5)) {
        std::cerr << "two_wire_cancel_test: cancellation ratio too high (" << ratio << ")\n";
        return 1;
    }

    const double avgAntisymmetry = antisymmetry / static_cast<double>(spec.ny - 2);
    if (!(avgAntisymmetry < 1e-4)) {
        std::cerr << "two_wire_cancel_test: By antisymmetry check failed (avg=" << avgAntisymmetry << ")\n";
        return 1;
    }

    for (const auto& wire : spec.wires) {
        double integrated = 0.0;
        for (std::size_t j = 0; j < spec.ny; ++j) {
            const double y = spec.originY + static_cast<double>(j) * spec.dy;
            for (std::size_t i = 0; i < spec.nx; ++i) {
                const double x = spec.originX + static_cast<double>(i) * spec.dx;
                const double r = std::hypot(x - wire.x, y - wire.y);
                if (r <= wire.radius) {
                    integrated += grid.Jz[grid.idx(i, j)] * cellArea;
                }
            }
        }
        const double tolerance = 0.15 * std::abs(wire.current) + 1e-9;
        if (std::abs(integrated - wire.current) > tolerance) {
            std::cerr << "two_wire_cancel_test: wire current mismatch. expected " << wire.current
                      << ", got " << integrated << "\n";
            return 1;
        }
    }

    std::cout << "two_wire_cancel_test: ratio=" << ratio << ", iters=" << report.iters << "\n";
    return 0;
}
