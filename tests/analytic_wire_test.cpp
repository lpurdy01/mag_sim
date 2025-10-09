// filename: analytic_wire_test.cpp
// part of 2D Electromagnetic Motor Simulator
// MIT License

#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"
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

    namespace fs = std::filesystem;
    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/analytic_wire_test.json").lexically_normal();
    ScenarioSpec spec = loadScenarioFromJson(scenarioPath.string());

    if (spec.wires.size() != 1) {
        std::cerr << "analytic_wire_test: scenario must contain exactly one wire\n";
        return 1;
    }

    Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    rasterizeScenarioToGrid(spec, grid);

    const double x0 = spec.originX;
    const double y0 = spec.originY;
    const double xc = spec.wires.front().x;
    const double yc = spec.wires.front().y;
    const double current = spec.wires.front().current;

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
    sampleMagnitudes.reserve(spec.nx);
    for (std::size_t j = 0; j < spec.ny; ++j) {
        for (std::size_t i = 0; i < spec.nx; ++i) {
            const double x = x0 + static_cast<double>(i) * spec.dx;
            const double y = y0 + static_cast<double>(j) * spec.dy;
            const double r = std::hypot(x - xc, y - yc);
            if (std::abs(r - rsample) <= 0.5 * spec.dx) {
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

    const std::size_t midJ = spec.ny / 2;
    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> bmag;
    xs.reserve(spec.nx);
    ys.reserve(spec.nx);
    bmag.reserve(spec.nx);

    for (std::size_t i = 0; i < spec.nx; ++i) {
        const double x = x0 + static_cast<double>(i) * spec.dx;
        xs.push_back(x);
        ys.push_back(0.0);
        const std::size_t idx = grid.idx(i, midJ);
        bmag.push_back(std::hypot(grid.Bx[idx], grid.By[idx]));
    }

    try {
        const std::filesystem::path outDir{"outputs"};
        std::filesystem::create_directories(outDir);

        const std::filesystem::path lineCsvPath = outDir / "validation_wire_line.csv";
        write_csv_line_profile(lineCsvPath.string(), xs, ys, bmag);

        std::vector<double> gridXs;
        std::vector<double> gridYs;
        std::vector<double> gridBx;
        std::vector<double> gridBy;
        std::vector<double> gridBmag;
        gridXs.reserve(spec.nx * spec.ny);
        gridYs.reserve(spec.nx * spec.ny);
        gridBx.reserve(spec.nx * spec.ny);
        gridBy.reserve(spec.nx * spec.ny);
        gridBmag.reserve(spec.nx * spec.ny);

        for (std::size_t j = 0; j < spec.ny; ++j) {
            for (std::size_t i = 0; i < spec.nx; ++i) {
                const double x = x0 + static_cast<double>(i) * spec.dx;
                const double y = y0 + static_cast<double>(j) * spec.dy;
                const std::size_t idx = grid.idx(i, j);
                const double bx = grid.Bx[idx];
                const double by = grid.By[idx];

                gridXs.push_back(x);
                gridYs.push_back(y);
                gridBx.push_back(bx);
                gridBy.push_back(by);
                gridBmag.push_back(std::hypot(bx, by));
            }
        }

        const std::filesystem::path fieldCsvPath = outDir / "validation_wire_field.csv";
        write_csv_field_map(fieldCsvPath.string(), gridXs, gridYs,
                            {{"Bx", &gridBx}, {"By", &gridBy}, {"Bmag", &gridBmag}});
    } catch (const std::exception& ex) {
        std::cerr << "Warning: failed to write CSV artifacts: " << ex.what() << "\n";
    }

    std::cout << "Converged in " << report.iters << " iterations with relResidual="
              << report.relResidual << ", relErr=" << relErr << "\n";

    return 0;
}
