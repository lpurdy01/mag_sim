#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "motorsim/grid.hpp"

namespace motorsim {

struct ScenarioSpec {
    enum class BoundaryType { Dirichlet, Neumann };

    struct Wire {
        double x{0.0};
        double y{0.0};
        double radius{0.0};
        double current{0.0};
    };

    struct Material {
        std::string name;
        double mu_r{1.0};
    };

    struct HalfspaceRegion {
        double normal_x{0.0};
        double normal_y{0.0};
        double offset{0.0};
        double mu_r{1.0};
        double inv_mu{0.0};
    };

    struct PolygonRegion {
        std::vector<double> xs;
        std::vector<double> ys;
        double mu_r{1.0};
        double inv_mu{0.0};
        double min_x{0.0};
        double max_x{0.0};
        double min_y{0.0};
        double max_y{0.0};
    };

    struct RegionMask {
        enum class Kind { Halfspace, Polygon };

        Kind kind{Kind::Halfspace};
        std::size_t index{0};
    };

    struct MagnetRegion {
        enum class Shape { Polygon, Rect };

        Shape shape{Shape::Polygon};
        std::vector<double> xs;
        std::vector<double> ys;
        double min_x{0.0};
        double max_x{0.0};
        double min_y{0.0};
        double max_y{0.0};
        double Mx{0.0};
        double My{0.0};
    };

    struct Outputs {
        struct FieldMap {
            std::string id;
            std::string path;
            std::string quantity;  // "B", "H", "BH", or "energy_density"
            std::string format;    // e.g., "csv"
        };

        struct LineProbe {
            std::string id;
            std::string path;
            std::string axis;      // "x" or "y"
            double value{0.0};     // coordinate along the fixed axis
            std::string quantity;  // "Bmag", "Bx", "By", "Hx", "Hy", "Hmag", or "energy_density"
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
    std::vector<Material> materials;
    std::vector<HalfspaceRegion> halfspaces;
    std::vector<PolygonRegion> polygons;
    std::vector<RegionMask> regionMasks;
    std::vector<Wire> wires;
    std::vector<MagnetRegion> magnetRegions;
    BoundaryType boundaryType{BoundaryType::Dirichlet};
    Outputs outputs;
};

ScenarioSpec loadScenarioFromJson(const std::string& path);

void rasterizeScenarioToGrid(const ScenarioSpec& spec, Grid2D& grid);

}  // namespace motorsim
