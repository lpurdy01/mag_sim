// filename: polygon_region_test.cpp
// part of 2D Electromagnetic Motor Simulator
// MIT License

#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"
#include "motorsim/types.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>

int main() {
    using namespace motorsim;

    namespace fs = std::filesystem;
    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/polygon_region_test.json").lexically_normal();

    ScenarioSpec spec;
    try {
        spec = loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to parse polygon scenario: " << ex.what() << "\n";
        return 1;
    }

    Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    try {
        rasterizeScenarioToGrid(spec, grid);
    } catch (const std::exception& ex) {
        std::cerr << "Rasterisation failed: " << ex.what() << "\n";
        return 1;
    }

    const double invMuAir = 1.0 / (MU0 * 1.0);
    const double invMuIron = 1.0 / (MU0 * 1500.0);

    auto cellIndex = [&](double x, double y) -> std::size_t {
        const double fi = std::round((x - spec.originX) / spec.dx);
        const double fj = std::round((y - spec.originY) / spec.dy);
        const long long i = static_cast<long long>(fi);
        const long long j = static_cast<long long>(fj);
        if (i < 0 || j < 0 ||
            i >= static_cast<long long>(spec.nx) ||
            j >= static_cast<long long>(spec.ny)) {
            throw std::runtime_error("Requested coordinate lies outside the grid");
        }
        return grid.idx(static_cast<std::size_t>(i), static_cast<std::size_t>(j));
    };

    try {
        const std::size_t boreIdx = cellIndex(0.0, 0.0);
        const std::size_t ringIdx = cellIndex(0.022, 0.0);
        const std::size_t outsideIdx = cellIndex(0.045, 0.0);

        const double boreInv = grid.invMu[boreIdx];
        const double ringInv = grid.invMu[ringIdx];
        const double outsideInv = grid.invMu[outsideIdx];

        if (std::abs(boreInv - invMuAir) > 1e-12) {
            std::cerr << "Inner bore did not retain air permeability" << "\n";
            return 1;
        }
        if (std::abs(ringInv - invMuIron) > 1e-12) {
            std::cerr << "Ring polygon did not assign iron permeability" << "\n";
            return 1;
        }
        if (std::abs(outsideInv - invMuAir) > 1e-12) {
            std::cerr << "Exterior reverted permeability unexpectedly" << "\n";
            return 1;
        }
    } catch (const std::exception& ex) {
        std::cerr << "Polygon raster validation failed: " << ex.what() << "\n";
        return 1;
    }

    std::cout << "Polygon region rasterisation validated successfully\n";
    return 0;
}
