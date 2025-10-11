// filename: analytic_wire_test.cpp
// part of 2D Electromagnetic Motor Simulator
// MIT License

#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"
#include "motorsim/io_csv.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <algorithm>
#include <array>
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

    double baselineResidual = 0.0;
    double baselineError = 0.0;
    bool baselineRecorded = false;

    std::vector<double> lineXs;
    std::vector<double> lineYs;
    std::vector<double> lineBmag;
    std::vector<double> fieldXs;
    std::vector<double> fieldYs;
    std::vector<double> fieldBx;
    std::vector<double> fieldBy;
    std::vector<double> fieldBmag;

    const double rsample = 0.05;
    const std::array<std::pair<SolverKind, const char*>, 2> solverCases = {
        std::make_pair(SolverKind::SOR, "SOR"), std::make_pair(SolverKind::CG, "CG")};

    for (const auto& [solverKind, label] : solverCases) {
        options.kind = solverKind;
        Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
        rasterizeScenarioToGrid(spec, grid);

        const SolveResult report = solveAz(grid, options);
        if (!report.converged) {
            std::cerr << label << " solver did not converge\n";
            return 1;
        }
        if (report.relResidual >= options.tol) {
            std::cerr << label << " residual too high: " << report.relResidual << "\n";
            return 1;
        }

        computeB(grid);

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
            std::cerr << label << ": insufficient sampling points on validation ring\n";
            return 1;
        }

        const double avgB = std::accumulate(sampleMagnitudes.begin(), sampleMagnitudes.end(), 0.0) /
                            static_cast<double>(sampleMagnitudes.size());
        const double analyticB = MU0 * current / (2.0 * M_PI * rsample);
        const double relErr = std::abs(avgB - analyticB) / analyticB;

        if (relErr >= 0.25) {
            std::cerr << label << ": relative error too high: " << relErr << "\n";
            return 1;
        }

        std::cout << label << " converged in " << report.iters
                  << " iterations with relResidual=" << report.relResidual
                  << ", relErr=" << relErr << "\n";

        if (!baselineRecorded) {
            baselineResidual = report.relResidual;
            baselineError = relErr;
            baselineRecorded = true;

            const std::size_t midJ = spec.ny / 2;
            lineXs.clear();
            lineYs.clear();
            lineBmag.clear();
            lineXs.reserve(spec.nx);
            lineYs.reserve(spec.nx);
            lineBmag.reserve(spec.nx);
            for (std::size_t i = 0; i < spec.nx; ++i) {
                const double x = x0 + static_cast<double>(i) * spec.dx;
                lineXs.push_back(x);
                lineYs.push_back(0.0);
                const std::size_t idx = grid.idx(i, midJ);
                lineBmag.push_back(std::hypot(grid.Bx[idx], grid.By[idx]));
            }

            fieldXs.clear();
            fieldYs.clear();
            fieldBx.clear();
            fieldBy.clear();
            fieldBmag.clear();
            fieldXs.reserve(spec.nx * spec.ny);
            fieldYs.reserve(spec.nx * spec.ny);
            fieldBx.reserve(spec.nx * spec.ny);
            fieldBy.reserve(spec.nx * spec.ny);
            fieldBmag.reserve(spec.nx * spec.ny);
            for (std::size_t j = 0; j < spec.ny; ++j) {
                for (std::size_t i = 0; i < spec.nx; ++i) {
                    const double x = x0 + static_cast<double>(i) * spec.dx;
                    const double y = y0 + static_cast<double>(j) * spec.dy;
                    const std::size_t idx = grid.idx(i, j);
                    const double bx = grid.Bx[idx];
                    const double by = grid.By[idx];
                    fieldXs.push_back(x);
                    fieldYs.push_back(y);
                    fieldBx.push_back(bx);
                    fieldBy.push_back(by);
                    fieldBmag.push_back(std::hypot(bx, by));
                }
            }
        } else {
            const double residualBound = baselineResidual * 1.1 + 1e-12;
            if (report.relResidual > residualBound) {
                std::cerr << label << ": residual " << report.relResidual
                          << " exceeds parity bound " << residualBound << "\n";
                return 1;
            }
            if (relErr > baselineError + 0.02) {
                std::cerr << label << ": analytic error " << relErr
                          << " exceeds allowed slack relative to baseline " << baselineError << "\n";
                return 1;
            }
        }
    }

    try {
        const std::filesystem::path outDir{"outputs"};
        std::filesystem::create_directories(outDir);

        const std::filesystem::path lineCsvPath = outDir / "validation_wire_line.csv";
        write_csv_line_profile(lineCsvPath.string(), lineXs, lineYs, lineBmag);

        const std::filesystem::path fieldCsvPath = outDir / "validation_wire_field.csv";
        write_csv_field_map(fieldCsvPath.string(), fieldXs, fieldYs,
                            {{"Bx", &fieldBx}, {"By", &fieldBy}, {"Bmag", &fieldBmag}});
    } catch (const std::exception& ex) {
        std::cerr << "Warning: failed to write CSV artifacts: " << ex.what() << "\n";
    }

    return 0;
}
