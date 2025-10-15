// filename: solver_benchmark.cpp
// part of 2D Electromagnetic Motor Simulator
// MIT License

#include "motorsim/grid.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct BenchmarkConfig {
    std::size_t nx{129};
    std::size_t ny{129};
    double lengthX{0.2};
    double lengthY{0.2};
    double current{10.0};
    double coreRadiusCells{3.0};
    std::size_t repeats{3};
    std::size_t maxIters{5000};
    double tol{1e-6};
    double omega{1.7};
    bool computeB{false};
    bool verbose{false};
    bool writeCsv{false};
    std::string csvPath{};
    std::string solver{"sor"};
};

void printUsage() {
    std::cout << "solver_benchmark options:\n"
              << "  --nx <int>             Number of cells in x-direction (default 129)\n"
              << "  --ny <int>             Number of cells in y-direction (default 129)\n"
              << "  --lx <float>           Domain length in x (meters, default 0.2)\n"
              << "  --ly <float>           Domain length in y (meters, default 0.2)\n"
              << "  --current <float>      Total current in wire core (A, default 10)\n"
              << "  --core-cells <float>   Wire core radius measured in cells (default 3)\n"
              << "  --max-iters <int>      Maximum solver iterations (default 5000)\n"
              << "  --tol <float>          Relative residual tolerance (default 1e-6)\n"
              << "  --omega <float>        SOR relaxation factor (default 1.7)\n"
              << "  --repeats <int>        Number of benchmark repeats (default 3)\n"
              << "  --compute-b            Include computeB() in timing\n"
              << "  --verbose              Print per-iteration residuals\n"
              << "  --solver {sor|cg}      Choose solver backend (default sor)\n"
              << "  --csv <path>           Append benchmark results to CSV file\n"
              << "  --help                 Show this message\n";
}

bool parseArgs(int argc, char** argv, BenchmarkConfig& cfg) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        try {
            if (arg == "--help") {
                printUsage();
                return false;
            } else if (arg == "--nx" && i + 1 < argc) {
                cfg.nx = static_cast<std::size_t>(std::stoul(argv[++i]));
            } else if (arg == "--ny" && i + 1 < argc) {
                cfg.ny = static_cast<std::size_t>(std::stoul(argv[++i]));
            } else if (arg == "--lx" && i + 1 < argc) {
                cfg.lengthX = std::stod(argv[++i]);
            } else if (arg == "--ly" && i + 1 < argc) {
                cfg.lengthY = std::stod(argv[++i]);
            } else if (arg == "--current" && i + 1 < argc) {
                cfg.current = std::stod(argv[++i]);
            } else if (arg == "--core-cells" && i + 1 < argc) {
                cfg.coreRadiusCells = std::stod(argv[++i]);
            } else if (arg == "--max-iters" && i + 1 < argc) {
                cfg.maxIters = static_cast<std::size_t>(std::stoul(argv[++i]));
            } else if (arg == "--tol" && i + 1 < argc) {
                cfg.tol = std::stod(argv[++i]);
            } else if (arg == "--omega" && i + 1 < argc) {
                cfg.omega = std::stod(argv[++i]);
            } else if (arg == "--repeats" && i + 1 < argc) {
                cfg.repeats = static_cast<std::size_t>(std::stoul(argv[++i]));
            } else if (arg == "--compute-b") {
                cfg.computeB = true;
            } else if (arg == "--verbose") {
                cfg.verbose = true;
            } else if (arg == "--solver" && i + 1 < argc) {
                cfg.solver = argv[++i];
            } else if (arg == "--csv" && i + 1 < argc) {
                cfg.writeCsv = true;
                cfg.csvPath = argv[++i];
            } else {
                std::cerr << "Unknown argument: " << arg << "\n";
                return false;
            }
        } catch (const std::exception& ex) {
            std::cerr << "Failed to parse argument " << arg << ": " << ex.what() << "\n";
            return false;
        }
    }
    return true;
}

void writeCsvResult(const BenchmarkConfig& cfg,
                    std::size_t cells,
                    const motorsim::SolveResult& report,
                    double avgMs,
                    double minMs,
                    double maxMs,
                    double updatesPerSecond) {
    namespace fs = std::filesystem;
    const fs::path csvPath{cfg.csvPath};
    const bool newFile = !fs::exists(csvPath);
    std::ofstream csv(csvPath, std::ios::app);
    if (!csv) {
        throw std::runtime_error("Failed to open CSV file: " + cfg.csvPath);
    }
    if (newFile) {
        csv << "nx,ny,length_x,length_y,cells,iters,tol,omega,avg_ms,min_ms,max_ms,updates_per_second\n";
    }
    csv << cfg.nx << ',' << cfg.ny << ',' << cfg.lengthX << ',' << cfg.lengthY << ',' << cells << ','
        << report.iters << ',' << cfg.tol << ',' << cfg.omega << ',' << avgMs << ',' << minMs << ','
        << maxMs << ',' << updatesPerSecond << '\n';
}

}  // namespace

int main(int argc, char** argv) {
    BenchmarkConfig cfg{};
    if (!parseArgs(argc, argv, cfg)) {
        return 1;
    }

    if (cfg.nx < 3 || cfg.ny < 3) {
        std::cerr << "Grid dimensions must be at least 3x3.\n";
        return 1;
    }

    const double dx = cfg.lengthX / static_cast<double>(cfg.nx - 1);
    const double dy = cfg.lengthY / static_cast<double>(cfg.ny - 1);
    motorsim::Grid2D grid(cfg.nx, cfg.ny, dx, dy);
    const double invMu0 = 1.0 / motorsim::MU0;
    std::fill(grid.invMu.begin(), grid.invMu.end(), invMu0);

    const double x0 = -0.5 * cfg.lengthX;
    const double y0 = -0.5 * cfg.lengthY;
    const double rc = cfg.coreRadiusCells * std::min(dx, dy);
    const double jzCore = cfg.current / (M_PI * rc * rc);

    for (std::size_t j = 0; j < cfg.ny; ++j) {
        for (std::size_t i = 0; i < cfg.nx; ++i) {
            const double x = x0 + static_cast<double>(i) * dx;
            const double y = y0 + static_cast<double>(j) * dy;
            const double r = std::hypot(x, y);
            if (r <= rc) {
                grid.Jz[grid.idx(i, j)] = jzCore;
            } else {
                grid.Jz[grid.idx(i, j)] = 0.0;
            }
        }
    }

    motorsim::SolveOptions options{};
    options.maxIters = cfg.maxIters;
    options.tol = cfg.tol;
    options.omega = cfg.omega;
    options.verbose = cfg.verbose;
    if (cfg.solver == "cg") {
        options.kind = motorsim::SolverKind::CG;
    } else {
        options.kind = motorsim::SolverKind::SOR;
    }

    std::vector<double> durationsMs;
    durationsMs.reserve(cfg.repeats);
    motorsim::SolveResult report{};

    for (std::size_t repeat = 0; repeat < cfg.repeats; ++repeat) {
        std::fill(grid.Az.begin(), grid.Az.end(), 0.0);

        const auto start = std::chrono::steady_clock::now();
        if (options.kind == motorsim::SolverKind::SOR) {
            report = motorsim::solveAz_GS_SOR(grid, options);
        } else {
            report = motorsim::solveAz_CG(grid, options);
        }
        if (cfg.computeB) {
            motorsim::computeB(grid);
        }
        const auto end = std::chrono::steady_clock::now();
        const double ms = std::chrono::duration<double, std::milli>(end - start).count();
        durationsMs.push_back(ms);
    }

    const std::size_t cells = cfg.nx * cfg.ny;
    const double avgMs = std::accumulate(durationsMs.begin(), durationsMs.end(), 0.0) /
                         static_cast<double>(durationsMs.size());
    const auto [minIt, maxIt] = std::minmax_element(durationsMs.begin(), durationsMs.end());
    const double minMs = *minIt;
    const double maxMs = *maxIt;

    const double cellIterations = static_cast<double>(cells) * static_cast<double>(report.iters);
    const double avgSeconds = avgMs / 1000.0;
    const double updatesPerSecond = cellIterations / avgSeconds;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Grid: " << cfg.nx << " x " << cfg.ny << " (" << cells << " cells)\n";
    std::cout << "Iterations: " << report.iters << " (tol=" << std::scientific << std::setprecision(2)
              << options.tol << std::fixed << std::setprecision(3)
              << ", omega=" << options.omega << ")\n";
    std::cout << "Average solve time: " << avgMs << " ms (min=" << minMs
              << " ms, max=" << maxMs << " ms)\n";
    std::cout << "Cell-iterations: " << cellIterations / 1.0e6 << "e6 total, "
              << (avgMs / cellIterations) * 1.0e6 << " ms per million cell-iterations\n";
    std::cout << "Throughput: " << updatesPerSecond / 1.0e6
              << "e6 cell-iterations/s\n";
    std::cout << "Converged: " << (report.converged ? "true" : "false")
              << ", relResidual=" << std::scientific << std::setprecision(2) << report.relResidual
              << std::fixed << std::setprecision(3) << '\n';

    if (cfg.writeCsv) {
        try {
            writeCsvResult(cfg, cells, report, avgMs, minMs, maxMs, updatesPerSecond);
        } catch (const std::exception& ex) {
            std::cerr << "Warning: " << ex.what() << "\n";
        }
    }

    return report.converged ? 0 : 2;
}

