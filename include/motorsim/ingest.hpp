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

    struct Outputs {
        struct FieldMap {
            std::string id;
            std::string path;
            std::string quantity;  // e.g., "B"
            std::string format;    // e.g., "csv"
        };

        struct LineProbe {
            std::string id;
            std::string path;
            std::string axis;      // "x" or "y"
            double value{0.0};     // coordinate along the fixed axis
            std::string quantity;  // e.g., "Bmag", "Bx", "By"
            std::string format;    // e.g., "csv"
        };

        std::vector<FieldMap> fieldMaps;
        std::vector<LineProbe> lineProbes;
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
    Outputs outputs;
};

ScenarioSpec loadScenarioFromJson(const std::string& path);

void rasterizeScenarioToGrid(const ScenarioSpec& spec, Grid2D& grid);

}  // namespace motorsim
