#include "motorsim/ingest.hpp"
#include "motorsim/io_csv.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <optional>
#include <sstream>
#include <system_error>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

void printUsage() {
    std::cout << "Usage: motor_sim [--scenario PATH] [--solve] [--write-midline]"
                 " [--list-outputs] [--outputs IDs] [--max-iters N] [--tol VALUE]"
                 " [--omega VALUE]\n";
}

std::vector<std::string> splitCommaSeparated(const std::string& text) {
    std::vector<std::string> items;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) {
            items.push_back(item);
        }
    }
    return items;
}

void ensureParentDirectory(const std::filesystem::path& path) {
    const std::filesystem::path parent = path.parent_path();
    if (!parent.empty()) {
        std::error_code ec;
        std::filesystem::create_directories(parent, ec);
    }
}

}  // namespace

int main(int argc, char** argv) {
    using namespace motorsim;

    std::optional<std::string> scenarioPath;
    bool runSolve = false;
    bool writeMidline = false;
    bool listOutputs = false;
    std::optional<std::string> outputsFilterArg;
    std::optional<int> maxItersOverride;
    std::optional<double> tolOverride;
    std::optional<double> omegaOverride;

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
        } else if (arg == "--list-outputs") {
            listOutputs = true;
        } else if (arg == "--outputs") {
            if (i + 1 >= argc) {
                std::cerr << "--outputs requires a comma-separated list or 'all'/'none'\n";
                printUsage();
                return 1;
            }
            outputsFilterArg = std::string(argv[++i]);
        } else if (arg == "--max-iters") {
            if (i + 1 >= argc) {
                std::cerr << "--max-iters requires an integer argument\n";
                printUsage();
                return 1;
            }
            int value = 0;
            try {
                value = std::stoi(argv[++i]);
            } catch (const std::exception&) {
                std::cerr << "--max-iters requires a valid integer argument\n";
                return 1;
            }
            if (value <= 0) {
                std::cerr << "--max-iters must be positive\n";
                return 1;
            }
            maxItersOverride = value;
        } else if (arg == "--tol") {
            if (i + 1 >= argc) {
                std::cerr << "--tol requires a positive floating-point argument\n";
                printUsage();
                return 1;
            }
            double value = 0.0;
            try {
                value = std::stod(argv[++i]);
            } catch (const std::exception&) {
                std::cerr << "--tol requires a valid floating-point argument\n";
                return 1;
            }
            if (!(value > 0.0)) {
                std::cerr << "--tol must be positive\n";
                return 1;
            }
            tolOverride = value;
        } else if (arg == "--omega") {
            if (i + 1 >= argc) {
                std::cerr << "--omega requires a floating-point argument\n";
                printUsage();
                return 1;
            }
            double value = 0.0;
            try {
                value = std::stod(argv[++i]);
            } catch (const std::exception&) {
                std::cerr << "--omega requires a valid floating-point argument\n";
                return 1;
            }
            if (!(value > 0.0)) {
                std::cerr << "--omega must be positive\n";
                return 1;
            }
            omegaOverride = value;
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

    std::unordered_set<std::string> availableOutputIds;
    availableOutputIds.reserve(spec.outputs.fieldMaps.size() + spec.outputs.lineProbes.size());
    for (const auto& request : spec.outputs.fieldMaps) {
        availableOutputIds.insert(request.id);
    }
    for (const auto& request : spec.outputs.lineProbes) {
        availableOutputIds.insert(request.id);
    }

    if (listOutputs) {
        std::cout << "Scenario outputs (" << availableOutputIds.size() << "):\n";
        for (const auto& request : spec.outputs.fieldMaps) {
            std::cout << "  - [field_map] id=" << request.id << ", quantity=" << request.quantity
                      << ", format=" << request.format << ", path=" << request.path << "\n";
        }
        for (const auto& request : spec.outputs.lineProbes) {
            std::cout << "  - [line_probe] id=" << request.id << ", axis=" << request.axis
                      << ", value=" << request.value << ", quantity=" << request.quantity
                      << ", format=" << request.format << ", path=" << request.path << "\n";
        }
    }

    std::unordered_set<std::string> requestedOutputIds;
    bool restrictOutputs = false;
    bool skipOutputs = false;
    if (outputsFilterArg) {
        if (outputsFilterArg->empty()) {
            std::cerr << "--outputs argument must not be empty\n";
            return 1;
        }
        if (*outputsFilterArg == "all") {
            restrictOutputs = false;
        } else if (*outputsFilterArg == "none") {
            skipOutputs = true;
        } else {
            if (availableOutputIds.empty()) {
                std::cerr << "Scenario defines no outputs but --outputs was specified\n";
                return 1;
            }
            const auto ids = splitCommaSeparated(*outputsFilterArg);
            if (ids.empty()) {
                std::cerr << "--outputs requires a comma-separated list of ids or 'all'/'none'\n";
                return 1;
            }
            restrictOutputs = true;
            for (const auto& id : ids) {
                if (availableOutputIds.find(id) == availableOutputIds.end()) {
                    std::cerr << "Requested output id not found in scenario: " << id << "\n";
                    return 1;
                }
                requestedOutputIds.insert(id);
            }
        }
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
    options.maxIters = maxItersOverride.value_or(20000);
    options.tol = tolOverride.value_or(1e-6);
    options.omega = omegaOverride.value_or(1.7);

    const SolveReport report = solveAz_GS_SOR(grid, options);
    if (!report.converged) {
        std::cerr << "Solver did not converge. relResidual=" << report.relResidual << "\n";
        return 1;
    }

    computeB(grid);
    std::cout << "Solve complete. iters=" << report.iters << ", relResidual=" << report.relResidual
              << "\n";

    const bool scenarioHasOutputs = !spec.outputs.fieldMaps.empty() || !spec.outputs.lineProbes.empty();
    const auto shouldEmit = [&](const std::string& id) {
        if (skipOutputs) {
            return false;
        }
        if (!restrictOutputs) {
            return true;
        }
        return requestedOutputIds.find(id) != requestedOutputIds.end();
    };

    bool fieldMapNeedsH = false;
    bool fieldMapNeedsEnergy = false;
    bool lineProbeNeedsH = false;
    bool lineProbeNeedsEnergy = false;
    if (scenarioHasOutputs && !skipOutputs) {
        for (const auto& request : spec.outputs.fieldMaps) {
            if (!shouldEmit(request.id)) {
                continue;
            }
            if (request.quantity == "H" || request.quantity == "BH") {
                fieldMapNeedsH = true;
            } else if (request.quantity == "energy_density") {
                fieldMapNeedsH = true;
                fieldMapNeedsEnergy = true;
            }
        }
        for (const auto& request : spec.outputs.lineProbes) {
            if (!shouldEmit(request.id)) {
                continue;
            }
            if (request.quantity == "Hx" || request.quantity == "Hy" || request.quantity == "Hmag") {
                lineProbeNeedsH = true;
            } else if (request.quantity == "energy_density") {
                lineProbeNeedsH = true;
                lineProbeNeedsEnergy = true;
            }
        }
    }

    const bool needHField = fieldMapNeedsH || fieldMapNeedsEnergy || lineProbeNeedsH ||
                            lineProbeNeedsEnergy;
    if (needHField) {
        computeH(grid);
    }

    bool outputError = false;
    if (scenarioHasOutputs && !skipOutputs) {
        std::vector<double> gridXs;
        std::vector<double> gridYs;
        std::vector<double> gridBx;
        std::vector<double> gridBy;
        std::vector<double> gridBmag;
        std::vector<double> gridHx;
        std::vector<double> gridHy;
        std::vector<double> gridHmag;
        std::vector<double> gridEnergyDensity;
        const bool fieldMapRequiresHData = fieldMapNeedsH || fieldMapNeedsEnergy;
        bool fieldPrepared = false;
        const auto prepareFieldVectors = [&]() {
            if (fieldPrepared) {
                return;
            }
            const std::size_t total = spec.nx * spec.ny;
            gridXs.resize(total);
            gridYs.resize(total);
            gridBx.resize(total);
            gridBy.resize(total);
            gridBmag.resize(total);
            if (fieldMapRequiresHData) {
                gridHx.resize(total);
                gridHy.resize(total);
                gridHmag.resize(total);
                if (fieldMapNeedsEnergy) {
                    gridEnergyDensity.resize(total);
                }
            }
            std::size_t idx = 0;
            for (std::size_t j = 0; j < spec.ny; ++j) {
                const double y = spec.originY + static_cast<double>(j) * spec.dy;
                for (std::size_t i = 0; i < spec.nx; ++i) {
                    const double x = spec.originX + static_cast<double>(i) * spec.dx;
                    const std::size_t gridIdx = grid.idx(i, j);
                    gridXs[idx] = x;
                    gridYs[idx] = y;
                    gridBx[idx] = grid.Bx[gridIdx];
                    gridBy[idx] = grid.By[gridIdx];
                    gridBmag[idx] = std::hypot(grid.Bx[gridIdx], grid.By[gridIdx]);
                    if (fieldMapRequiresHData) {
                        const double hx = grid.Hx[gridIdx];
                        const double hy = grid.Hy[gridIdx];
                        gridHx[idx] = hx;
                        gridHy[idx] = hy;
                        gridHmag[idx] = std::hypot(hx, hy);
                        if (fieldMapNeedsEnergy) {
                            gridEnergyDensity[idx] =
                                0.5 * (grid.Bx[gridIdx] * hx + grid.By[gridIdx] * hy);
                        }
                    }
                    ++idx;
                }
            }
            fieldPrepared = true;
        };

        for (const auto& request : spec.outputs.fieldMaps) {
            if (!shouldEmit(request.id)) {
                continue;
            }
            prepareFieldVectors();
            const std::filesystem::path outPath{request.path};
            ensureParentDirectory(outPath);
            try {
                std::vector<FieldColumnView> columns;
                if (request.quantity == "B") {
                    columns = {{"Bx", &gridBx}, {"By", &gridBy}, {"Bmag", &gridBmag}};
                } else if (request.quantity == "H") {
                    columns = {{"Hx", &gridHx}, {"Hy", &gridHy}, {"Hmag", &gridHmag}};
                } else if (request.quantity == "BH") {
                    columns = {{"Bx", &gridBx}, {"By", &gridBy}, {"Bmag", &gridBmag},
                               {"Hx", &gridHx}, {"Hy", &gridHy}, {"Hmag", &gridHmag}};
                } else {  // energy_density
                    columns = {{"Bx", &gridBx},        {"By", &gridBy},
                               {"Bmag", &gridBmag},   {"Hx", &gridHx},
                               {"Hy", &gridHy},       {"Hmag", &gridHmag},
                               {"EnergyDensity", &gridEnergyDensity}};
                }
                write_csv_field_map(outPath.string(), gridXs, gridYs, columns);
                std::cout << "Wrote field_map output '" << request.id << "' to " << outPath << "\n";
            } catch (const std::exception& ex) {
                std::cerr << "Failed to write field_map output '" << request.id << "': " << ex.what()
                          << "\n";
                outputError = true;
            }
        }

        for (const auto& request : spec.outputs.lineProbes) {
            if (!shouldEmit(request.id)) {
                continue;
            }

            std::vector<double> xs;
            std::vector<double> ys;
            std::vector<double> values;

            if (request.axis == "x") {
                const double coord = (request.value - spec.originX) / spec.dx;
                const auto idx = static_cast<long long>(std::llround(coord));
                if (idx < 0 || idx >= static_cast<long long>(spec.nx)) {
                    std::cerr << "line_probe '" << request.id << "' is outside the domain (x="
                              << request.value << ")\n";
                    outputError = true;
                    continue;
                }
                const double actual = spec.originX + static_cast<double>(idx) * spec.dx;
                if (std::abs(actual - request.value) > 0.5 * spec.dx + 1e-9) {
                    std::cerr << "line_probe '" << request.id
                              << "' does not align with grid column (requested x=" << request.value
                              << ", snapped to " << actual << ")\n";
                    outputError = true;
                    continue;
                }

                xs.resize(spec.ny, actual);
                ys.resize(spec.ny);
                values.resize(spec.ny);
                for (std::size_t j = 0; j < spec.ny; ++j) {
                    const double y = spec.originY + static_cast<double>(j) * spec.dy;
                    const std::size_t gridIdx = grid.idx(static_cast<std::size_t>(idx), j);
                    ys[j] = y;
                    if (request.quantity == "Bx") {
                        values[j] = grid.Bx[gridIdx];
                    } else if (request.quantity == "By") {
                        values[j] = grid.By[gridIdx];
                    } else if (request.quantity == "Bmag") {
                        values[j] = std::hypot(grid.Bx[gridIdx], grid.By[gridIdx]);
                    } else if (request.quantity == "Hx") {
                        values[j] = grid.Hx[gridIdx];
                    } else if (request.quantity == "Hy") {
                        values[j] = grid.Hy[gridIdx];
                    } else if (request.quantity == "Hmag") {
                        values[j] = std::hypot(grid.Hx[gridIdx], grid.Hy[gridIdx]);
                    } else {  // energy_density
                        values[j] =
                            0.5 * (grid.Bx[gridIdx] * grid.Hx[gridIdx] +
                                   grid.By[gridIdx] * grid.Hy[gridIdx]);
                    }
                }
            } else if (request.axis == "y") {
                const double coord = (request.value - spec.originY) / spec.dy;
                const auto idx = static_cast<long long>(std::llround(coord));
                if (idx < 0 || idx >= static_cast<long long>(spec.ny)) {
                    std::cerr << "line_probe '" << request.id << "' is outside the domain (y="
                              << request.value << ")\n";
                    outputError = true;
                    continue;
                }
                const double actual = spec.originY + static_cast<double>(idx) * spec.dy;
                if (std::abs(actual - request.value) > 0.5 * spec.dy + 1e-9) {
                    std::cerr << "line_probe '" << request.id
                              << "' does not align with grid row (requested y=" << request.value
                              << ", snapped to " << actual << ")\n";
                    outputError = true;
                    continue;
                }

                xs.resize(spec.nx);
                ys.resize(spec.nx, actual);
                values.resize(spec.nx);
                for (std::size_t i = 0; i < spec.nx; ++i) {
                    const double x = spec.originX + static_cast<double>(i) * spec.dx;
                    const std::size_t gridIdx = grid.idx(i, static_cast<std::size_t>(idx));
                    xs[i] = x;
                    if (request.quantity == "Bx") {
                        values[i] = grid.Bx[gridIdx];
                    } else if (request.quantity == "By") {
                        values[i] = grid.By[gridIdx];
                    } else if (request.quantity == "Bmag") {
                        values[i] = std::hypot(grid.Bx[gridIdx], grid.By[gridIdx]);
                    } else if (request.quantity == "Hx") {
                        values[i] = grid.Hx[gridIdx];
                    } else if (request.quantity == "Hy") {
                        values[i] = grid.Hy[gridIdx];
                    } else if (request.quantity == "Hmag") {
                        values[i] = std::hypot(grid.Hx[gridIdx], grid.Hy[gridIdx]);
                    } else {  // energy_density
                        values[i] =
                            0.5 * (grid.Bx[gridIdx] * grid.Hx[gridIdx] +
                                   grid.By[gridIdx] * grid.Hy[gridIdx]);
                    }
                }
            }

            const std::filesystem::path outPath{request.path};
            ensureParentDirectory(outPath);
            try {
                write_csv_line_profile(outPath.string(), xs, ys, values);
                std::cout << "Wrote line_probe output '" << request.id << "' to " << outPath << "\n";
            } catch (const std::exception& ex) {
                std::cerr << "Failed to write line_probe output '" << request.id << "': " << ex.what()
                          << "\n";
                outputError = true;
            }
        }
    }

    if (writeMidline) {
        const std::filesystem::path outDir{"outputs"};

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
        ensureParentDirectory(csvPath);
        try {
            write_csv_line_profile(csvPath.string(), xs, ys, bmag);
            std::cout << "Wrote midline profile to " << csvPath << "\n";
        } catch (const std::exception& ex) {
            std::cerr << "Failed to write midline CSV: " << ex.what() << "\n";
        }
    }

    return outputError ? 1 : 0;
}
