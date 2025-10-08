#include "motorsim/ingest.hpp"
#include "motorsim/io_csv.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <optional>
#include <string>
#include <vector>

namespace {

void printUsage() {
    std::cout << "Usage: motor_sim [--scenario PATH] [--solve] [--write-midline]\n";
}

}  // namespace

int main(int argc, char** argv) {
    using namespace motorsim;

    std::optional<std::string> scenarioPath;
    bool runSolve = false;
    bool writeMidline = false;

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

    Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    try {
        rasterizeScenarioToGrid(spec, grid);
    } catch (const std::exception& ex) {
        std::cerr << "Rasterisation error: " << ex.what() << "\n";
        return 1;
    }

    if (!runSolve) {
        std::cout << "Scenario loaded. Use --solve to run the magnetostatic solver." << std::endl;
        return 0;
    }

    SolveOptions options{};
    options.maxIters = 5000;
    options.tol = 1e-6;
    options.omega = 1.7;

    const SolveReport report = solveAz_GS_SOR(grid, options);
    if (!report.converged) {
        std::cerr << "Solver did not converge. relResidual=" << report.relResidual << "\n";
        return 1;
    }

    computeB(grid);
    std::cout << "Solve complete. iters=" << report.iters << ", relResidual=" << report.relResidual << "\n";

    if (writeMidline) {
        const std::filesystem::path outDir{"outputs"};
        std::error_code ec;
        std::filesystem::create_directories(outDir, ec);

        const std::size_t iMid = spec.nx / 2;
        std::vector<double> xs;
        std::vector<double> ys;
        std::vector<double> bmag;
        xs.reserve(spec.ny);
        ys.reserve(spec.ny);
        bmag.reserve(spec.ny);
        for (std::size_t j = 0; j < spec.ny; ++j) {
            const double x = spec.originX + static_cast<double>(iMid) * spec.dx;
            const double y = spec.originY + static_cast<double>(j) * spec.dy;
            const std::size_t idx = grid.idx(iMid, j);
            xs.push_back(x);
            ys.push_back(y);
            bmag.push_back(std::hypot(grid.Bx[idx], grid.By[idx]));
        }

        const std::filesystem::path csvPath = outDir / "twowire_midline.csv";
        try {
            write_csv_line_profile(csvPath.string(), xs, ys, bmag);
            std::cout << "Wrote midline profile to " << csvPath << "\n";
        } catch (const std::exception& ex) {
            std::cerr << "Failed to write midline CSV: " << ex.what() << "\n";
        }
    }

    return 0;
}
