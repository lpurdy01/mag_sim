// filename: analytic_wire_test.cpp
// part of 2D Electromagnetic Motor Simulator
// MIT License

#include "motorsim/grid.hpp"
#include "motorsim/io_csv.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <numeric>
#include <vector>

int main() {
    using namespace motorsim;

    const double L = 0.2;
    const std::size_t nx = 201;
    const std::size_t ny = 201;
    const double dx = L / static_cast<double>(nx - 1);
    const double dy = L / static_cast<double>(ny - 1);

    Grid2D grid(nx, ny, dx, dy);
    const double invMu0 = 1.0 / MU0;
    std::fill(grid.invMu.begin(), grid.invMu.end(), invMu0);

    const double x0 = -0.5 * L;
    const double y0 = -0.5 * L;
    const double xc = 0.0;
    const double yc = 0.0;

    const double current = 10.0;
    // Slightly enlarge the source disk (3 cells) to reduce discretization error
    // from approximating the singular wire current with a uniform patch.
    const double rc = 3.0 * dx;
    const double jzCore = current / (M_PI * rc * rc);

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const double x = x0 + static_cast<double>(i) * dx;
            const double y = y0 + static_cast<double>(j) * dy;
            const double r = std::hypot(x - xc, y - yc);
            if (r <= rc) {
                grid.Jz[grid.idx(i, j)] = jzCore;
            }
        }
    }

    SolveOptions options{};
    options.maxIters = 5000;
    options.tol = 1e-6;
    options.omega = 1.9;
    options.verbose = false;

    const SolveReport report = solveAz_GS_SOR(grid, options);
    if (!report.converged) {
        std::cerr << "Solver did not converge\n";
        return 1;
    }
    if (report.relResidual >= 1e-6) {
        std::cerr << "Residual too high: " << report.relResidual << "\n";
        return 1;
    }

    computeB(grid);

    const double rsample = 0.05;
    std::vector<double> sampleMagnitudes;
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            const double x = x0 + static_cast<double>(i) * dx;
            const double y = y0 + static_cast<double>(j) * dy;
            const double r = std::hypot(x - xc, y - yc);
            if (std::abs(r - rsample) <= 0.5 * dx) {
                const std::size_t idx = grid.idx(i, j);
                const double bx = grid.Bx[idx];
                const double by = grid.By[idx];
                sampleMagnitudes.push_back(std::hypot(bx, by));
            }
        }
    }

    if (sampleMagnitudes.size() < 10) {
        std::cerr << "Insufficient sampling points on validation ring\n";
        return 1;
    }

    const double avgB = std::accumulate(sampleMagnitudes.begin(), sampleMagnitudes.end(), 0.0) /
                        static_cast<double>(sampleMagnitudes.size());
    const double analyticB = MU0 * current / (2.0 * M_PI * rsample);
    const double relErr = std::abs(avgB - analyticB) / analyticB;

    if (relErr >= 0.25) {
        std::cerr << "Relative error too high: " << relErr << "\n";
        return 1;
    }

    const std::size_t midJ = ny / 2;
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> bmag;
    xs.reserve(nx);
    ys.reserve(nx);
    bmag.reserve(nx);

    for (std::size_t i = 0; i < nx; ++i) {
        const double x = x0 + static_cast<double>(i) * dx;
        xs.push_back(x);
        ys.push_back(0.0);
        const std::size_t idx = grid.idx(i, midJ);
        bmag.push_back(std::hypot(grid.Bx[idx], grid.By[idx]));
    }

    try {
        const std::filesystem::path outDir{"outputs"};
        std::filesystem::create_directories(outDir);
        const std::filesystem::path csvPath = outDir / "validation_wire_line.csv";
        write_csv_line_profile(csvPath.string(), xs, ys, bmag);
    } catch (const std::exception& ex) {
        std::cerr << "Warning: failed to write CSV profile: " << ex.what() << "\n";
    }

    std::cout << "Converged in " << report.iters << " iterations with relResidual="
              << report.relResidual << ", relErr=" << relErr << "\n";

    return 0;
}
