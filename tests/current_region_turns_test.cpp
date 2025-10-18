// MIT License

#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>

int main() {
    using namespace motorsim;

    namespace fs = std::filesystem;
    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/current_region_turns_test.json").lexically_normal();

    ScenarioSpec spec;
    try {
        spec = loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to parse current-region scenario: " << ex.what() << "\n";
        return 1;
    }

    if (spec.currentRegions.empty()) {
        std::cerr << "Scenario did not define any current regions\n";
        return 1;
    }

    Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    try {
        rasterizeScenarioToGrid(spec, grid);
    } catch (const std::exception& ex) {
        std::cerr << "Rasterisation failed: " << ex.what() << "\n";
        return 1;
    }

    const auto& region = spec.currentRegions.front();
    if (region.copperArea <= 0.0) {
        std::cerr << "Copper area was not populated for the current region\n";
        return 1;
    }

    const double expectedCopperArea = region.area * region.fillFraction;
    if (std::abs(region.copperArea - expectedCopperArea) > 1e-9) {
        std::cerr << "Copper area mismatch: expected " << expectedCopperArea << " got " << region.copperArea << "\n";
        return 1;
    }

    double totalCurrent = 0.0;
    const double cellArea = spec.dx * spec.dy;
    for (std::size_t j = 0; j < spec.ny; ++j) {
        for (std::size_t i = 0; i < spec.nx; ++i) {
            totalCurrent += grid.Jz[grid.idx(i, j)] * cellArea;
        }
    }

    const double expectedCurrent = region.current * region.turns * region.fillFraction;
    const double tolerance = 0.05 * std::abs(expectedCurrent);
    if (std::abs(totalCurrent - expectedCurrent) > tolerance) {
        std::cerr << "Integrated current density " << totalCurrent << " differs from expected " << expectedCurrent
                  << " by more than 5%\n";
        return 1;
    }

    std::cout << "Current region turns scaling verified successfully\n";
    return 0;
}
