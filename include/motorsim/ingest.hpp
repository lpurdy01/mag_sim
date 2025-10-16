#pragma once

#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "motorsim/grid.hpp"
#include "motorsim/types.hpp"

namespace motorsim {

struct ScenarioSpec {
    enum class BoundaryType { Dirichlet, Neumann };

    struct Wire {
        double x{0.0};
        double y{0.0};
        double radius{0.0};
        double current{0.0};
    };

    struct CurrentRegion {
        std::string id;
        std::string phase;
        double orientation{1.0};
        std::vector<double> xs;
        std::vector<double> ys;
        double min_x{0.0};
        double max_x{0.0};
        double min_y{0.0};
        double max_y{0.0};
        double area{0.0};
        double current{0.0};
    };

    struct Rotor {
        std::string name;
        double pivotX{0.0};
        double pivotY{0.0};
        std::vector<std::size_t> polygonIndices;
        std::vector<std::size_t> magnetIndices;
        std::vector<std::size_t> wireIndices;
    };

    struct Material {
        std::string name;
        double mu_r{1.0};
        double sigma{0.0};
    };

    struct HalfspaceRegion {
        double normal_x{0.0};
        double normal_y{0.0};
        double offset{0.0};
        double mu_r{1.0};
        double inv_mu{0.0};
        double sigma{0.0};
    };

    struct PolygonRegion {
        std::vector<double> xs;
        std::vector<double> ys;
        double mu_r{1.0};
        double inv_mu{0.0};
        double sigma{0.0};
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

        struct VtkSeries {
            std::string id;
            std::string directory;
            std::string basename;
            bool includeB{true};
            bool includeH{true};
            bool includeEnergy{false};
        };

        struct LineProbe {
            std::string id;
            std::string path;
            std::string axis;      // "x" or "y"
            double value{0.0};     // coordinate along the fixed axis
            std::string quantity;  // "Bmag", "Bx", "By", "Hx", "Hy", "Hmag", or "energy_density"
            std::string format;    // e.g., "csv"
        };

        struct Probe {
            enum class Quantity { Force, Torque, ForceAndTorque };
            enum class Method { StressTensor };

            std::string id;
            std::string path;
            Quantity quantity{Quantity::Torque};
            Method method{Method::StressTensor};
            std::vector<double> loopXs;
            std::vector<double> loopYs;
        };

        struct BackEmfProbe {
            enum class RegionShape { Polygon, Rect };

            std::string id;
            std::string path;
            FluxComponent component{FluxComponent::Bmag};
            RegionShape shape{RegionShape::Polygon};
            std::vector<double> xs;
            std::vector<double> ys;
            double minX{0.0};
            double maxX{0.0};
            double minY{0.0};
            double maxY{0.0};
            std::vector<std::size_t> frameIndices;
        };

        struct PolylineOutlines {
            std::string id;
            std::string path;
        };

        struct BoreAverageProbe {
            std::string id;
            std::vector<double> xs;
            std::vector<double> ys;
            std::string path;
        };

        struct MechanicalTrace {
            std::string id;
            std::string path;
            std::vector<std::string> rotors;
        };

        std::vector<FieldMap> fieldMaps;
        std::vector<VtkSeries> vtkSeries;
        std::vector<LineProbe> lineProbes;
        std::vector<Probe> probes;
        std::vector<BackEmfProbe> backEmfProbes;
        std::vector<PolylineOutlines> polylineOutlines;
        std::vector<BoreAverageProbe> boreProbes;
        std::vector<MechanicalTrace> mechanicalTraces;
    };

    struct MechanicalSystem {
        struct Rotor {
            std::string name;
            std::size_t rotorIndex{std::numeric_limits<std::size_t>::max()};
            double inertia{0.0};
            double damping{0.0};
            double loadTorque{0.0};
            double initialAngleDeg{0.0};
            bool hasInitialAngle{false};
            double initialSpeedRadPerSec{0.0};
            bool hasInitialSpeed{false};
            std::string torqueProbeId;
            bool hasTorqueProbe{false};
        };

        std::vector<Rotor> rotors;
    };

    struct Circuit {
        struct Resistor {
            std::string id;
            std::size_t nodePos{0};
            std::size_t nodeNeg{0};
            double resistance{0.0};
        };

        struct Inductor {
            std::string id;
            std::size_t nodePos{0};
            std::size_t nodeNeg{0};
            double inductance{0.0};
            double initialCurrent{0.0};
        };

        struct VoltageSource {
            std::string id;
            std::size_t nodePos{0};
            std::size_t nodeNeg{0};
            double value{0.0};
        };

        struct CoilLink {
            std::string id;
            std::size_t inductorIndex{0};
            std::size_t regionIndex{0};
            double turns{0.0};
        };

        std::string id;
        std::vector<std::string> nodes;
        std::vector<Resistor> resistors;
        std::vector<Inductor> inductors;
        std::vector<VoltageSource> voltageSources;
        std::vector<CoilLink> coilLinks;
    };

    struct TimelineFrame {
        struct WireOverride {
            std::size_t index{0};
            double current{0.0};
        };

        struct CurrentRegionOverride {
            std::size_t index{0};
            double current{0.0};
        };

        struct MagnetOverride {
            std::size_t index{0};
            bool hasAngle{false};
            double angleDegrees{0.0};
            bool hasVector{false};
            double magnetizationX{0.0};
            double magnetizationY{0.0};
        };

        struct RotorAngleOverride {
            std::size_t index{0};
            double angleDegrees{0.0};
        };

        struct VoltageSourceOverride {
            std::size_t circuitIndex{0};
            std::size_t sourceIndex{0};
            double value{0.0};
        };

        double time{0.0};
        bool hasRotorAngle{false};
        double rotorAngleDeg{0.0};
        std::vector<RotorAngleOverride> rotorAngles;
        std::vector<WireOverride> wireOverrides;
        std::vector<CurrentRegionOverride> currentRegionOverrides;
        std::vector<MagnetOverride> magnetOverrides;
        std::vector<VoltageSourceOverride> voltageSourceOverrides;
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
    double sigma_background{0.0};
    std::vector<Material> materials;
    std::vector<HalfspaceRegion> halfspaces;
    std::vector<PolygonRegion> polygons;
    std::vector<RegionMask> regionMasks;
    std::vector<Wire> wires;
    std::vector<CurrentRegion> currentRegions;
    std::vector<Rotor> rotors;
    std::vector<MagnetRegion> magnetRegions;
    std::vector<Circuit> circuits;
    std::optional<MechanicalSystem> mechanical;
    BoundaryType boundaryType{BoundaryType::Dirichlet};
    Outputs outputs;
    std::vector<TimelineFrame> timeline;

    struct SolverSettings {
        bool solverSpecified{false};
        std::string solverId{"sor"};
        bool warmStartSpecified{false};
        bool warmStart{false};
        bool prolongationSpecified{false};
        bool prolongationEnabled{false};
        std::optional<std::size_t> prolongationNx;
        std::optional<std::size_t> prolongationNy;
        bool progressEverySecSpecified{false};
        double progressEverySec{0.0};
        bool snapshotEveryItersSpecified{false};
        std::size_t snapshotEveryIters{0};
        bool quietSpecified{false};
        bool quiet{false};
        bool harmonicFrequencySpecified{false};
        double harmonicFrequencyHz{0.0};
        bool harmonicOmegaSpecified{false};
        double harmonicOmega{0.0};
    } solverSettings;
};

ScenarioSpec loadScenarioFromJson(const std::string& path);

void rasterizeScenarioToGrid(const ScenarioSpec& spec, Grid2D& grid);

struct ScenarioFrame {
    std::size_t index{0};
    double time{0.0};
    ScenarioSpec spec;
};

std::vector<ScenarioFrame> expandScenarioTimeline(const ScenarioSpec& spec);

}  // namespace motorsim
