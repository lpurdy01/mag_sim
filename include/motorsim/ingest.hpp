#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "motorsim/grid.hpp"

namespace motorsim {

struct ScenarioSpec {
    struct Wire {
        double x{0.0};
        double y{0.0};
        double radius{0.0};
        double current{0.0};
    };

    std::string version;
    double Lx{0.0};
    double Ly{0.0};
    std::size_t nx{0};
    std::size_t ny{0};
    double originX{0.0};
    double originY{0.0};
    double dx{0.0};
    double dy{0.0};
    double mu_r_background{1.0};
    std::vector<Wire> wires;
};

ScenarioSpec loadScenarioFromJson(const std::string& path);

void rasterizeScenarioToGrid(const ScenarioSpec& spec, Grid2D& grid);

}  // namespace motorsim
