#include "motorsim/ingest.hpp"

#include "motorsim/types.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <nlohmann/json.hpp>

namespace motorsim {
namespace {
constexpr double kTiny = 1e-12;
constexpr double kPi = 3.14159265358979323846;

bool pointInPolygonCoords(double x, double y, const std::vector<double>& xs,
                          const std::vector<double>& ys) {
    const std::size_t count = xs.size();
    if (count == 0U) {
        return false;
    }

    bool inside = false;
    for (std::size_t i = 0, j = count - 1; i < count; j = i++) {
        const double xi = xs[i];
        const double yi = ys[i];
        const double xj = xs[j];
        const double yj = ys[j];

        const double denom = yj - yi;
        if (std::abs(denom) < kTiny) {
            continue;
        }
        const bool intersects = ((yi > y) != (yj > y)) &&
                                (x < (xj - xi) * (y - yi) / denom + xi);
        if (intersects) {
            inside = !inside;
        }
    }

    return inside;
}

bool pointInPolygon(double x, double y, const ScenarioSpec::PolygonRegion& poly) {
    return pointInPolygonCoords(x, y, poly.xs, poly.ys);
}

bool pointInRect(double x, double y, double minX, double maxX, double minY, double maxY) {
    return x >= minX && x <= maxX && y >= minY && y <= maxY;
}

ScenarioSpec::BoundaryType parseBoundaryType(const std::string& value) {
    if (value == "dirichlet") {
        return ScenarioSpec::BoundaryType::Dirichlet;
    }
    if (value == "neumann") {
        return ScenarioSpec::BoundaryType::Neumann;
    }
    throw std::runtime_error("Unsupported boundary type: " + value);
}

double requirePositive(const std::string& field, double value) {
    if (!(value > 0.0)) {
        throw std::runtime_error(field + " must be positive");
    }
    return value;
}

std::size_t requireGridSize(const std::string& field, std::size_t value) {
    if (value < 2) {
        throw std::runtime_error(field + " must be at least 2");
    }
    return value;
}

std::string requireNonEmpty(const std::string& field, const std::string& value) {
    if (value.empty()) {
        throw std::runtime_error(field + " must be a non-empty string");
    }
    return value;
}

double polygonArea(const std::vector<double>& xs, const std::vector<double>& ys) {
    if (xs.size() != ys.size() || xs.size() < 3) {
        return 0.0;
    }
    double area = 0.0;
    const std::size_t count = xs.size();
    for (std::size_t i = 0, j = count - 1; i < count; j = i++) {
        area += (xs[j] + xs[i]) * (ys[j] - ys[i]);
    }
    return 0.5 * std::abs(area);
}
}  // namespace

ScenarioSpec loadScenarioFromJson(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Failed to open scenario JSON: " + path);
    }

    nlohmann::json json;
    input >> json;

    ScenarioSpec spec{};
    spec.version = json.value("version", std::string{});
    if (spec.version.empty()) {
        throw std::runtime_error("Scenario JSON missing required field: version");
    }
    if (spec.version != "0.1" && spec.version != "0.2") {
        throw std::runtime_error("Unsupported scenario version: " + spec.version);
    }
    const bool allowAdvancedRegions = (spec.version != "0.1");

    const std::string units = json.value("units", std::string{"SI"});
    if (units != "SI") {
        throw std::runtime_error("Unsupported units: " + units + ". Only SI is supported.");
    }

    const auto& domain = json.at("domain");
    spec.Lx = requirePositive("domain.Lx", domain.at("Lx").get<double>());
    spec.Ly = requirePositive("domain.Ly", domain.at("Ly").get<double>());
    spec.nx = requireGridSize("domain.nx", domain.at("nx").get<std::size_t>());
    spec.ny = requireGridSize("domain.ny", domain.at("ny").get<std::size_t>());

    spec.dx = spec.Lx / static_cast<double>(spec.nx - 1);
    spec.dy = spec.Ly / static_cast<double>(spec.ny - 1);
    spec.originX = -0.5 * spec.Lx;
    spec.originY = -0.5 * spec.Ly;

    spec.boundaryType = ScenarioSpec::BoundaryType::Dirichlet;
    if (json.contains("boundary")) {
        const auto& boundary = json.at("boundary");
        const std::string type = boundary.value("type", std::string{"dirichlet"});
        spec.boundaryType = parseBoundaryType(type);
    }

    const auto& materials = json.at("materials");
    if (!materials.is_array() || materials.empty()) {
        throw std::runtime_error("Scenario must define at least one material");
    }

    spec.materials.clear();
    spec.halfspaces.clear();
    spec.polygons.clear();
    spec.regionMasks.clear();

    std::unordered_map<std::string, double> muMap;
    for (const auto& material : materials) {
        const std::string name = material.at("name").get<std::string>();
        if (muMap.count(name) != 0U) {
            throw std::runtime_error("Duplicate material name: " + name);
        }
        const double mu_r = requirePositive("material.mu_r", material.at("mu_r").get<double>());
        muMap.emplace(name, mu_r);
        ScenarioSpec::Material entry{};
        entry.name = name;
        entry.mu_r = mu_r;
        spec.materials.push_back(entry);
    }

    if (!spec.materials.empty()) {
        spec.mu_r_background = spec.materials.front().mu_r;
    } else {
        spec.mu_r_background = 1.0;
    }

    const auto& regions = json.at("regions");
    if (!regions.is_array() || regions.empty()) {
        throw std::runtime_error("Scenario must define at least one region");
    }

    bool foundUniform = false;
    for (const auto& region : regions) {
        const std::string type = region.at("type").get<std::string>();
        if (type == "uniform") {
            const std::string materialName = region.at("material").get<std::string>();
            const auto it = muMap.find(materialName);
            if (it == muMap.end()) {
                throw std::runtime_error("Region references unknown material: " + materialName);
            }
            if (spec.version == "0.1" && foundUniform) {
                throw std::runtime_error(
                    "Multiple uniform regions defined; only one is supported in v0.1");
            }
            spec.mu_r_background = it->second;
            foundUniform = true;
        } else if (type == "halfspace") {
            if (!allowAdvancedRegions) {
                throw std::runtime_error("Region type 'halfspace' requires scenario version 0.2");
            }
            const std::string materialName = region.at("material").get<std::string>();
            const auto it = muMap.find(materialName);
            if (it == muMap.end()) {
                throw std::runtime_error("Region references unknown material: " + materialName);
            }
            const auto& normal = region.at("normal");
            if (!normal.is_array() || normal.size() != 2) {
                throw std::runtime_error("halfspace region requires a 2-element 'normal' array");
            }
            double nx = normal.at(0).get<double>();
            double ny = normal.at(1).get<double>();
            const double length = std::hypot(nx, ny);
            if (!(length > 0.0)) {
                throw std::runtime_error("halfspace region normal vector must be non-zero");
            }
            const double invLength = 1.0 / length;
            ScenarioSpec::HalfspaceRegion hs{};
            hs.normal_x = nx * invLength;
            hs.normal_y = ny * invLength;
            hs.offset = region.at("offset").get<double>() * invLength;
            hs.mu_r = it->second;
            hs.inv_mu = 1.0 / (MU0 * hs.mu_r);
            spec.halfspaces.push_back(hs);
            ScenarioSpec::RegionMask mask{};
            mask.kind = ScenarioSpec::RegionMask::Kind::Halfspace;
            mask.index = spec.halfspaces.size() - 1;
            spec.regionMasks.push_back(mask);
        } else if (type == "polygon") {
            if (!allowAdvancedRegions) {
                throw std::runtime_error("Region type 'polygon' requires scenario version 0.2");
            }
            const std::string materialName = region.at("material").get<std::string>();
            const auto it = muMap.find(materialName);
            if (it == muMap.end()) {
                throw std::runtime_error("Region references unknown material: " + materialName);
            }
            const auto& vertices = region.at("vertices");
            if (!vertices.is_array() || vertices.size() < 3) {
                throw std::runtime_error("polygon region requires an array of at least three vertices");
            }
            ScenarioSpec::PolygonRegion poly{};
            poly.mu_r = it->second;
            poly.inv_mu = 1.0 / (MU0 * poly.mu_r);
            poly.xs.reserve(vertices.size());
            poly.ys.reserve(vertices.size());
            double minX = std::numeric_limits<double>::infinity();
            double maxX = -std::numeric_limits<double>::infinity();
            double minY = std::numeric_limits<double>::infinity();
            double maxY = -std::numeric_limits<double>::infinity();
            for (const auto& vertex : vertices) {
                if (!vertex.is_array() || vertex.size() != 2) {
                    throw std::runtime_error("polygon vertices must be [x, y] arrays");
                }
                const double vx = vertex.at(0).get<double>();
                const double vy = vertex.at(1).get<double>();
                poly.xs.push_back(vx);
                poly.ys.push_back(vy);
                minX = std::min(minX, vx);
                maxX = std::max(maxX, vx);
                minY = std::min(minY, vy);
                maxY = std::max(maxY, vy);
            }
            poly.min_x = minX;
            poly.max_x = maxX;
            poly.min_y = minY;
            poly.max_y = maxY;
            spec.polygons.push_back(std::move(poly));
            ScenarioSpec::RegionMask mask{};
            mask.kind = ScenarioSpec::RegionMask::Kind::Polygon;
            mask.index = spec.polygons.size() - 1;
            spec.regionMasks.push_back(mask);
        } else {
            throw std::runtime_error("Unsupported region type: " + type);
        }
    }

    if (spec.version == "0.1") {
        if (!foundUniform) {
            throw std::runtime_error("Scenario must define a uniform region in v0.1");
        }
    }

    const auto& sources = json.at("sources");
    if (!sources.is_array()) {
        throw std::runtime_error("Scenario sources must be an array");
    }

    spec.wires.clear();
    spec.currentRegions.clear();
    spec.wires.reserve(sources.size());
    spec.currentRegions.reserve(sources.size());
    std::unordered_map<std::string, std::size_t> currentRegionIdToIndex;
    std::unordered_map<std::string, std::vector<std::size_t>> phaseToCurrentRegions;
    for (const auto& source : sources) {
        const std::string type = source.at("type").get<std::string>();
        if (type == "wire") {
            ScenarioSpec::Wire wire{};
            wire.x = source.at("x").get<double>();
            wire.y = source.at("y").get<double>();
            wire.radius = requirePositive("wire.radius", source.at("radius").get<double>());
            wire.current = source.value("I", source.value("current", 0.0));
            spec.wires.push_back(wire);
        } else if (type == "current_region" || type == "current-region") {
            ScenarioSpec::CurrentRegion region{};
            if (source.contains("id")) {
                region.id = requireNonEmpty("current_region.id", source.at("id").get<std::string>());
                if (currentRegionIdToIndex.find(region.id) != currentRegionIdToIndex.end()) {
                    throw std::runtime_error("Duplicate current_region id: " + region.id);
                }
            }
            region.phase = source.value("phase", std::string{});
            region.orientation = source.value("orientation", source.value("polarity", 1.0));
            if (!(std::abs(region.orientation) > 0.0)) {
                throw std::runtime_error("current_region orientation must be non-zero");
            }
            const auto& vertices = source.at("vertices");
            if (!vertices.is_array() || vertices.size() < 3) {
                throw std::runtime_error("current_region requires an array of at least three vertices");
            }
            region.xs.reserve(vertices.size());
            region.ys.reserve(vertices.size());
            double minX = std::numeric_limits<double>::infinity();
            double maxX = -std::numeric_limits<double>::infinity();
            double minY = std::numeric_limits<double>::infinity();
            double maxY = -std::numeric_limits<double>::infinity();
            for (const auto& vertex : vertices) {
                if (!vertex.is_array() || vertex.size() != 2) {
                    throw std::runtime_error("current_region vertices must be [x, y] arrays");
                }
                const double vx = vertex.at(0).get<double>();
                const double vy = vertex.at(1).get<double>();
                region.xs.push_back(vx);
                region.ys.push_back(vy);
                minX = std::min(minX, vx);
                maxX = std::max(maxX, vx);
                minY = std::min(minY, vy);
                maxY = std::max(maxY, vy);
            }
            region.min_x = minX;
            region.max_x = maxX;
            region.min_y = minY;
            region.max_y = maxY;
            region.area = polygonArea(region.xs, region.ys);
            if (!(region.area > kTiny)) {
                throw std::runtime_error("current_region polygon area must be positive");
            }
            region.current = source.value("I", source.value("current", 0.0));
            const std::size_t index = spec.currentRegions.size();
            spec.currentRegions.push_back(region);
            if (!region.id.empty()) {
                currentRegionIdToIndex.emplace(region.id, index);
            }
            if (!region.phase.empty()) {
                phaseToCurrentRegions[region.phase].push_back(index);
            }
        } else {
            throw std::runtime_error("Unsupported source type: " + type);
        }
    }

    spec.magnetRegions.clear();
    if (json.contains("magnet_regions")) {
        const auto& magnets = json.at("magnet_regions");
        if (!magnets.is_array()) {
            throw std::runtime_error("magnet_regions must be an array");
        }
        for (const auto& magnet : magnets) {
            ScenarioSpec::MagnetRegion region{};
            const std::string shape = magnet.at("type").get<std::string>();

            const auto& magnetisation = magnet.at("magnetization");
            if (magnetisation.is_array()) {
                if (magnetisation.size() != 2) {
                    throw std::runtime_error("magnetization array must contain exactly two elements");
                }
                region.Mx = magnetisation.at(0).get<double>();
                region.My = magnetisation.at(1).get<double>();
            } else if (magnetisation.is_object()) {
                region.Mx = magnetisation.value("Mx", 0.0);
                region.My = magnetisation.value("My", 0.0);
            } else {
                throw std::runtime_error("magnetization must be an array [Mx, My] or object with Mx/My");
            }

            if (shape == "polygon") {
                region.shape = ScenarioSpec::MagnetRegion::Shape::Polygon;
                const auto& vertices = magnet.at("vertices");
                if (!vertices.is_array() || vertices.size() < 3) {
                    throw std::runtime_error(
                        "magnet polygon region requires an array of at least three vertices");
                }
                double minX = std::numeric_limits<double>::infinity();
                double maxX = -std::numeric_limits<double>::infinity();
                double minY = std::numeric_limits<double>::infinity();
                double maxY = -std::numeric_limits<double>::infinity();
                region.xs.reserve(vertices.size());
                region.ys.reserve(vertices.size());
                for (const auto& vertex : vertices) {
                    if (!vertex.is_array() || vertex.size() != 2) {
                        throw std::runtime_error("magnet polygon vertices must be [x, y] arrays");
                    }
                    const double vx = vertex.at(0).get<double>();
                    const double vy = vertex.at(1).get<double>();
                    region.xs.push_back(vx);
                    region.ys.push_back(vy);
                    minX = std::min(minX, vx);
                    maxX = std::max(maxX, vx);
                    minY = std::min(minY, vy);
                    maxY = std::max(maxY, vy);
                }
                region.min_x = minX;
                region.max_x = maxX;
                region.min_y = minY;
                region.max_y = maxY;
            } else if (shape == "rect") {
                region.shape = ScenarioSpec::MagnetRegion::Shape::Rect;
                const auto& xrange = magnet.at("x_range");
                const auto& yrange = magnet.at("y_range");
                if (!xrange.is_array() || xrange.size() != 2) {
                    throw std::runtime_error("rect magnet region requires x_range [xmin, xmax]");
                }
                if (!yrange.is_array() || yrange.size() != 2) {
                    throw std::runtime_error("rect magnet region requires y_range [ymin, ymax]");
                }
                region.min_x = xrange.at(0).get<double>();
                region.max_x = xrange.at(1).get<double>();
                region.min_y = yrange.at(0).get<double>();
                region.max_y = yrange.at(1).get<double>();
                if (!(region.min_x < region.max_x) || !(region.min_y < region.max_y)) {
                    throw std::runtime_error("rect magnet region requires min < max in both ranges");
                }
                region.xs = {region.min_x, region.max_x, region.max_x, region.min_x};
                region.ys = {region.min_y, region.min_y, region.max_y, region.max_y};
            } else {
                throw std::runtime_error("Unsupported magnet region type: " + shape);
            }

            spec.magnetRegions.push_back(std::move(region));
        }
    }

    spec.rotors.clear();
    std::unordered_map<std::string, std::size_t> rotorNameToIndex;
    if (json.contains("rotors")) {
        const auto& rotorsNode = json.at("rotors");
        if (!rotorsNode.is_array()) {
            throw std::runtime_error("rotors must be an array when provided");
        }
        spec.rotors.reserve(rotorsNode.size());

        const auto parsePivot = [](const nlohmann::json& node, double& px, double& py) {
            if (node.is_array()) {
                if (node.size() != 2) {
                    throw std::runtime_error("rotor pivot array must contain exactly two values");
                }
                px = node.at(0).get<double>();
                py = node.at(1).get<double>();
            } else if (node.is_object()) {
                px = node.value("x", node.value("px", 0.0));
                py = node.value("y", node.value("py", 0.0));
            } else {
                throw std::runtime_error("rotor pivot must be an array [x, y] or object with x/y fields");
            }
        };

        for (const auto& rotorNode : rotorsNode) {
            if (!rotorNode.is_object()) {
                throw std::runtime_error("rotor entries must be JSON objects");
            }
            ScenarioSpec::Rotor rotor{};
            rotor.name = rotorNode.value("name", std::string{});
            if (rotor.name.empty()) {
                rotor.name = "rotor_" + std::to_string(spec.rotors.size());
            }
            if (!rotorNameToIndex.emplace(rotor.name, spec.rotors.size()).second) {
                throw std::runtime_error("Duplicate rotor name: " + rotor.name);
            }

            if (rotorNode.contains("pivot")) {
                parsePivot(rotorNode.at("pivot"), rotor.pivotX, rotor.pivotY);
            } else if (rotorNode.contains("center")) {
                parsePivot(rotorNode.at("center"), rotor.pivotX, rotor.pivotY);
            } else {
                rotor.pivotX = rotorNode.value("pivot_x", rotorNode.value("center_x", 0.0));
                rotor.pivotY = rotorNode.value("pivot_y", rotorNode.value("center_y", 0.0));
            }

            const auto parseIndexList = [&](const nlohmann::json& arr, std::size_t count,
                                            const std::string& field, std::vector<std::size_t>& dest) {
                if (count == 0) {
                    throw std::runtime_error("rotor '" + rotor.name + "' references " + field +
                                             " but scenario defines none");
                }
                if (!arr.is_array()) {
                    throw std::runtime_error("rotor field '" + field + "' must be an array");
                }
                for (const auto& value : arr) {
                    std::size_t index = 0;
                    if (value.is_number_integer()) {
                        index = value.get<std::size_t>();
                    } else if (value.is_number_unsigned()) {
                        index = value.get<std::size_t>();
                    } else if (value.is_object()) {
                        if (value.contains("index")) {
                            index = value.at("index").get<std::size_t>();
                        } else if (value.contains("idx")) {
                            index = value.at("idx").get<std::size_t>();
                        } else {
                            throw std::runtime_error(
                                "rotor field '" + field + "' object entries must include 'index' or 'idx'");
                        }
                    } else {
                        throw std::runtime_error("rotor field '" + field + "' entries must be indices or objects");
                    }
                    if (index >= count) {
                        throw std::runtime_error("rotor field '" + field + "' index out of range");
                    }
                    dest.push_back(index);
                }
            };

            if (rotorNode.contains("polygon_indices")) {
                parseIndexList(rotorNode.at("polygon_indices"), spec.polygons.size(), "polygon_indices",
                               rotor.polygonIndices);
            }
            if (rotorNode.contains("polygons")) {
                parseIndexList(rotorNode.at("polygons"), spec.polygons.size(), "polygons", rotor.polygonIndices);
            }

            if (rotorNode.contains("magnet_indices")) {
                parseIndexList(rotorNode.at("magnet_indices"), spec.magnetRegions.size(), "magnet_indices",
                               rotor.magnetIndices);
            }
            if (rotorNode.contains("magnets")) {
                parseIndexList(rotorNode.at("magnets"), spec.magnetRegions.size(), "magnets",
                               rotor.magnetIndices);
            }

            if (rotorNode.contains("wire_indices")) {
                parseIndexList(rotorNode.at("wire_indices"), spec.wires.size(), "wire_indices",
                               rotor.wireIndices);
            }
            if (rotorNode.contains("wires")) {
                parseIndexList(rotorNode.at("wires"), spec.wires.size(), "wires", rotor.wireIndices);
            }

            spec.rotors.push_back(std::move(rotor));
        }
    }

    spec.outputs = ScenarioSpec::Outputs{};
    if (json.contains("outputs")) {
        const auto& outputs = json.at("outputs");
        if (!outputs.is_array()) {
            throw std::runtime_error("Scenario outputs must be an array");
        }

        std::unordered_set<std::string> ids;
        ids.reserve(outputs.size());
        for (const auto& output : outputs) {
            const std::string type = output.at("type").get<std::string>();
            const std::string id = requireNonEmpty("outputs.id", output.at("id").get<std::string>());
            if (!ids.insert(id).second) {
                throw std::runtime_error("Duplicate output id: " + id);
            }

            if (type == "field_map") {
                ScenarioSpec::Outputs::FieldMap request{};
                request.id = id;
                request.quantity = output.value("quantity", std::string{"B"});
                if (request.quantity != "B" && request.quantity != "H" &&
                    request.quantity != "BH" && request.quantity != "energy_density") {
                    throw std::runtime_error(
                        "Unsupported field_map quantity for output '" + id + "': " + request.quantity);
                }
                request.format = output.value("format", std::string{"csv"});
                if (request.format != "csv" && request.format != "vti") {
                    throw std::runtime_error(
                        "Unsupported field_map format for output '" + id + "': " + request.format);
                }
                if (output.contains("path")) {
                    request.path =
                        requireNonEmpty("outputs.path", output.at("path").get<std::string>());
                } else {
                    request.path = "outputs/" + id + (request.format == "csv" ? ".csv" : ".vti");
                }
                spec.outputs.fieldMaps.push_back(std::move(request));
            } else if (type == "vtk_field_series") {
                ScenarioSpec::Outputs::VtkSeries request{};
                request.id = id;
                request.directory = output.value("dir", std::string{"outputs"});
                request.basename = output.value("basename", id);
                if (request.basename.empty()) {
                    throw std::runtime_error("vtk_field_series basename must be non-empty for output '" + id + "'");
                }
                bool includeB = false;
                bool includeH = false;
                bool includeEnergy = false;
                if (output.contains("quantities")) {
                    const auto& quantities = output.at("quantities");
                    if (!quantities.is_array() || quantities.empty()) {
                        throw std::runtime_error(
                            "vtk_field_series quantities must be a non-empty array for output '" + id + "'");
                    }
                    for (const auto& quantity : quantities) {
                        const std::string value = quantity.get<std::string>();
                        if (value == "B") {
                            includeB = true;
                        } else if (value == "H") {
                            includeH = true;
                        } else if (value == "BH") {
                            includeB = true;
                            includeH = true;
                        } else if (value == "energy" || value == "energy_density" || value == "BH_energy") {
                            includeEnergy = true;
                            includeB = true;
                            includeH = true;
                        } else {
                            throw std::runtime_error(
                                "Unsupported vtk_field_series quantity for output '" + id + "': " + value);
                        }
                    }
                } else {
                    includeB = true;
                    includeH = true;
                }
                request.includeB = includeB;
                request.includeH = includeH;
                request.includeEnergy = includeEnergy;
                spec.outputs.vtkSeries.push_back(std::move(request));
            } else if (type == "line_probe") {
                ScenarioSpec::Outputs::LineProbe request{};
                request.id = id;
                request.quantity = output.value("quantity", std::string{"Bmag"});
                if (request.quantity != "Bmag" && request.quantity != "Bx" &&
                    request.quantity != "By" && request.quantity != "Hx" &&
                    request.quantity != "Hy" && request.quantity != "Hmag" &&
                    request.quantity != "energy_density") {
                    throw std::runtime_error(
                        "Unsupported line_probe quantity for output '" + id + "': " + request.quantity);
                }
                request.axis = requireNonEmpty("outputs.axis", output.value("axis", std::string{}));
                if (request.axis != "x" && request.axis != "y") {
                    throw std::runtime_error(
                        "line_probe axis must be 'x' or 'y' for output '" + id + "'");
                }
                request.value = output.at("value").get<double>();
                request.format = output.value("format", std::string{"csv"});
                if (request.format != "csv") {
                    throw std::runtime_error(
                        "Unsupported line_probe format for output '" + id + "': " + request.format);
                }
                if (output.contains("path")) {
                    request.path =
                        requireNonEmpty("outputs.path", output.at("path").get<std::string>());
                } else {
                    request.path = "outputs/" + id + ".csv";
                }
                spec.outputs.lineProbes.push_back(std::move(request));
            } else if (type == "probe") {
                ScenarioSpec::Outputs::Probe request{};
                request.id = id;
                const std::string quantity = output.value("probe_type", std::string{"torque"});
                if (quantity == "torque") {
                    request.quantity = ScenarioSpec::Outputs::Probe::Quantity::Torque;
                } else if (quantity == "force") {
                    request.quantity = ScenarioSpec::Outputs::Probe::Quantity::Force;
                } else if (quantity == "force_and_torque") {
                    request.quantity = ScenarioSpec::Outputs::Probe::Quantity::ForceAndTorque;
                } else {
                    throw std::runtime_error("Unsupported probe_type for output '" + id + "': " + quantity);
                }

                const std::string method = output.value("method", std::string{"stress_tensor"});
                if (method == "stress_tensor") {
                    request.method = ScenarioSpec::Outputs::Probe::Method::StressTensor;
                } else {
                    throw std::runtime_error("Unsupported probe method for output '" + id + "': " + method);
                }

                const auto parseVertices = [&](const nlohmann::json& vertices) {
                    if (!vertices.is_array() || vertices.size() < 3) {
                        throw std::runtime_error("Probe loop for output '" + id + "' requires at least three vertices");
                    }
                    request.loopXs.clear();
                    request.loopYs.clear();
                    request.loopXs.reserve(vertices.size());
                    request.loopYs.reserve(vertices.size());
                    for (const auto& vertex : vertices) {
                        if (!vertex.is_array() || vertex.size() != 2) {
                            throw std::runtime_error(
                                "Probe loop vertices must be [x, y] pairs for output '" + id + "'");
                        }
                        request.loopXs.push_back(vertex.at(0).get<double>());
                        request.loopYs.push_back(vertex.at(1).get<double>());
                    }
                };

                if (!output.contains("loop")) {
                    throw std::runtime_error("Probe output '" + id + "' is missing required loop definition");
                }
                const auto& loop = output.at("loop");
                if (loop.is_array()) {
                    parseVertices(loop);
                } else if (loop.is_object()) {
                    const std::string loopType = loop.value("type", std::string{"polygon"});
                    if (loopType != "polygon") {
                        throw std::runtime_error(
                            "Unsupported probe loop type for output '" + id + "': " + loopType);
                    }
                    if (!loop.contains("vertices")) {
                        throw std::runtime_error(
                            "Probe loop object for output '" + id + "' must contain a 'vertices' array");
                    }
                    parseVertices(loop.at("vertices"));
                } else {
                    throw std::runtime_error("Probe loop for output '" + id + "' must be an array or object");
                }

                if (output.contains("path")) {
                    request.path =
                        requireNonEmpty("outputs.path", output.at("path").get<std::string>());
                } else {
                    request.path = "outputs/" + id + ".csv";
                }

                spec.outputs.probes.push_back(std::move(request));
            } else if (type == "back_emf_probe") {
                ScenarioSpec::Outputs::BackEmfProbe request{};
                request.id = id;

                std::string component = output.value("component", std::string{"Bmag"});
                std::transform(component.begin(), component.end(), component.begin(),
                               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
                if (component == "bx") {
                    request.component = FluxComponent::Bx;
                } else if (component == "by") {
                    request.component = FluxComponent::By;
                } else if (component == "b" || component == "bmag" || component == "mag") {
                    request.component = FluxComponent::Bmag;
                } else {
                    throw std::runtime_error("Unsupported back_emf_probe component for output '" + id +
                                             "': " + component);
                }

                const auto parseFrames = [&](const nlohmann::json& frameNode) {
                    if (!frameNode.is_array()) {
                        throw std::runtime_error(
                            "back_emf_probe frames must be an array of frame indices for output '" + id + "'");
                    }
                    request.frameIndices.clear();
                    request.frameIndices.reserve(frameNode.size());
                    for (const auto& entry : frameNode) {
                        if (!entry.is_number_integer()) {
                            throw std::runtime_error(
                                "back_emf_probe frames entries must be integers for output '" + id + "'");
                        }
                        const int value = entry.get<int>();
                        if (value < 0) {
                            throw std::runtime_error(
                                "back_emf_probe frames entries must be non-negative for output '" + id + "'");
                        }
                        request.frameIndices.push_back(static_cast<std::size_t>(value));
                    }
                    if (request.frameIndices.size() < 2) {
                        throw std::runtime_error(
                            "back_emf_probe output '" + id + "' requires at least two frame indices");
                    }
                };

                if (output.contains("frames")) {
                    parseFrames(output.at("frames"));
                } else if (output.contains("frame_indices")) {
                    parseFrames(output.at("frame_indices"));
                }

                const auto parsePolygonVertices = [&](const nlohmann::json& vertices) {
                    if (!vertices.is_array() || vertices.size() < 3) {
                        throw std::runtime_error("back_emf_probe region for output '" + id +
                                                 "' requires at least three vertices");
                    }
                    request.shape = ScenarioSpec::Outputs::BackEmfProbe::RegionShape::Polygon;
                    request.xs.clear();
                    request.ys.clear();
                    request.xs.reserve(vertices.size());
                    request.ys.reserve(vertices.size());
                    double minX = std::numeric_limits<double>::infinity();
                    double maxX = -std::numeric_limits<double>::infinity();
                    double minY = std::numeric_limits<double>::infinity();
                    double maxY = -std::numeric_limits<double>::infinity();
                    for (const auto& vertex : vertices) {
                        if (!vertex.is_array() || vertex.size() != 2) {
                            throw std::runtime_error(
                                "back_emf_probe vertices must be [x, y] pairs for output '" + id + "'");
                        }
                        const double vx = vertex.at(0).get<double>();
                        const double vy = vertex.at(1).get<double>();
                        request.xs.push_back(vx);
                        request.ys.push_back(vy);
                        minX = std::min(minX, vx);
                        maxX = std::max(maxX, vx);
                        minY = std::min(minY, vy);
                        maxY = std::max(maxY, vy);
                    }
                    request.minX = minX;
                    request.maxX = maxX;
                    request.minY = minY;
                    request.maxY = maxY;
                };

                const auto parseRect = [&](const nlohmann::json& region) {
                    if (!region.contains("x_range") || !region.contains("y_range")) {
                        throw std::runtime_error(
                            "back_emf_probe rect region requires x_range and y_range for output '" + id + "'");
                    }
                    const auto& xr = region.at("x_range");
                    const auto& yr = region.at("y_range");
                    if (!xr.is_array() || xr.size() != 2 || !yr.is_array() || yr.size() != 2) {
                        throw std::runtime_error(
                            "back_emf_probe rect ranges must contain two entries for output '" + id + "'");
                    }
                    const double x0 = xr.at(0).get<double>();
                    const double x1 = xr.at(1).get<double>();
                    const double y0 = yr.at(0).get<double>();
                    const double y1 = yr.at(1).get<double>();
                    const double minX = std::min(x0, x1);
                    const double maxX = std::max(x0, x1);
                    const double minY = std::min(y0, y1);
                    const double maxY = std::max(y0, y1);
                    if (!(maxX > minX) || !(maxY > minY)) {
                        throw std::runtime_error("back_emf_probe rect ranges must span a positive area for output '" + id +
                                                 "'");
                    }
                    request.shape = ScenarioSpec::Outputs::BackEmfProbe::RegionShape::Rect;
                    request.minX = minX;
                    request.maxX = maxX;
                    request.minY = minY;
                    request.maxY = maxY;
                    request.xs.clear();
                    request.ys.clear();
                };

                bool regionParsed = false;
                if (output.contains("region")) {
                    const auto& region = output.at("region");
                    if (!region.is_object()) {
                        throw std::runtime_error("back_emf_probe region must be an object for output '" + id + "'");
                    }
                    const std::string regionType = region.value("type", std::string{"polygon"});
                    if (regionType == "polygon") {
                        if (!region.contains("vertices")) {
                            throw std::runtime_error(
                                "back_emf_probe polygon region requires vertices for output '" + id + "'");
                        }
                        parsePolygonVertices(region.at("vertices"));
                        regionParsed = true;
                    } else if (regionType == "rect" || regionType == "rectangle") {
                        parseRect(region);
                        regionParsed = true;
                    } else {
                        throw std::runtime_error("Unsupported back_emf_probe region type for output '" + id +
                                                 "': " + regionType);
                    }
                }

                if (!regionParsed && output.contains("loop")) {
                    const auto& loop = output.at("loop");
                    if (loop.is_array()) {
                        parsePolygonVertices(loop);
                        regionParsed = true;
                    } else if (loop.is_object()) {
                        const std::string loopType = loop.value("type", std::string{"polygon"});
                        if (loopType != "polygon") {
                            throw std::runtime_error(
                                "Unsupported back_emf_probe loop type for output '" + id + "': " + loopType);
                        }
                        if (!loop.contains("vertices")) {
                            throw std::runtime_error(
                                "back_emf_probe loop object for output '" + id + "' requires vertices");
                        }
                        parsePolygonVertices(loop.at("vertices"));
                        regionParsed = true;
                    } else {
                        throw std::runtime_error(
                            "back_emf_probe loop must be an array or object for output '" + id + "'");
                    }
                }

                if (!regionParsed) {
                    throw std::runtime_error(
                        "back_emf_probe output '" + id + "' requires a region or loop definition");
                }

                if (output.contains("path")) {
                    request.path = requireNonEmpty("outputs.path", output.at("path").get<std::string>());
                } else {
                    request.path = "outputs/" + id + "_emf.csv";
                }

                spec.outputs.backEmfProbes.push_back(std::move(request));
            } else if (type == "polyline_outlines") {
                ScenarioSpec::Outputs::PolylineOutlines request{};
                request.id = id;
                if (output.contains("path")) {
                    request.path = requireNonEmpty("outputs.path", output.at("path").get<std::string>());
                } else {
                    std::string dir = output.value("dir", std::string{"outputs"});
                    if (!dir.empty() && dir.back() != '/' && dir.back() != '\\') {
                        dir.push_back('/');
                    }
                    const std::string filename = output.value("filename", id + "_outlines.vtp");
                    request.path = dir + filename;
                }
                spec.outputs.polylineOutlines.push_back(std::move(request));
            } else if (type == "bore_avg_B" || type == "bore_avg_b") {
                ScenarioSpec::Outputs::BoreAverageProbe request{};
                request.id = id;
                const auto parseVertices = [&](const nlohmann::json& vertices) {
                    if (!vertices.is_array() || vertices.size() < 3) {
                        throw std::runtime_error(
                            "bore_avg_B vertices must contain at least three entries for output '" + id + "'");
                    }
                    request.xs.clear();
                    request.ys.clear();
                    request.xs.reserve(vertices.size());
                    request.ys.reserve(vertices.size());
                    for (const auto& vertex : vertices) {
                        if (!vertex.is_array() || vertex.size() != 2) {
                            throw std::runtime_error(
                                "bore_avg_B vertices must be [x, y] pairs for output '" + id + "'");
                        }
                        request.xs.push_back(vertex.at(0).get<double>());
                        request.ys.push_back(vertex.at(1).get<double>());
                    }
                };
                if (output.contains("vertices")) {
                    parseVertices(output.at("vertices"));
                } else if (output.contains("polygon")) {
                    const auto& node = output.at("polygon");
                    if (node.is_array()) {
                        parseVertices(node);
                    } else if (node.is_object() && node.contains("vertices")) {
                        parseVertices(node.at("vertices"));
                    } else {
                        throw std::runtime_error(
                            "bore_avg_B polygon definition must be an array or object with vertices for output '" + id + "'");
                    }
                } else {
                    throw std::runtime_error(
                        "bore_avg_B output '" + id + "' requires 'vertices' or 'polygon' definition");
                }
                if (request.xs.size() != request.ys.size() || request.xs.size() < 3) {
                    throw std::runtime_error("bore_avg_B output '" + id + "' requires a valid polygon region");
                }
                if (output.contains("path")) {
                    request.path = requireNonEmpty("outputs.path", output.at("path").get<std::string>());
                } else {
                    request.path = "outputs/" + id + "_bore.csv";
                }
                spec.outputs.boreProbes.push_back(std::move(request));
            } else {
                throw std::runtime_error("Unsupported output type: " + type);
            }
        }
    }

    spec.timeline.clear();
    if (json.contains("timeline")) {
        const auto& timeline = json.at("timeline");
        if (!timeline.is_array() || timeline.empty()) {
            throw std::runtime_error("timeline must be a non-empty array when provided");
        }
        spec.timeline.reserve(timeline.size());
        std::size_t frameIndex = 0;
        for (const auto& frame : timeline) {
            if (!frame.is_object()) {
                throw std::runtime_error("timeline entries must be JSON objects");
            }
            ScenarioSpec::TimelineFrame entry{};
            if (frame.contains("t")) {
                entry.time = frame.at("t").get<double>();
            } else if (frame.contains("time")) {
                entry.time = frame.at("time").get<double>();
            } else {
                entry.time = static_cast<double>(frameIndex);
            }

            if (frame.contains("rotor_angle")) {
                entry.hasRotorAngle = true;
                entry.rotorAngleDeg = frame.at("rotor_angle").get<double>();
            } else if (frame.contains("rotor_angle_deg")) {
                entry.hasRotorAngle = true;
                entry.rotorAngleDeg = frame.at("rotor_angle_deg").get<double>();
            }

            const auto addRotorAngle = [&](std::size_t rotorIndex, double angleDeg) {
                if (rotorIndex >= spec.rotors.size()) {
                    throw std::runtime_error("timeline rotor angle index out of range");
                }
                ScenarioSpec::TimelineFrame::RotorAngleOverride override{};
                override.index = rotorIndex;
                override.angleDegrees = angleDeg;
                entry.rotorAngles.push_back(override);
            };

            if (frame.contains("rotor_angles")) {
                if (spec.rotors.empty()) {
                    throw std::runtime_error(
                        "timeline rotor_angles provided but scenario defines no rotors");
                }
                const auto& rotorAnglesNode = frame.at("rotor_angles");
                if (rotorAnglesNode.is_array()) {
                    if (!rotorAnglesNode.empty() && rotorAnglesNode.front().is_number()) {
                        const std::size_t count =
                            std::min<std::size_t>(rotorAnglesNode.size(), spec.rotors.size());
                        for (std::size_t i = 0; i < count; ++i) {
                            addRotorAngle(i, rotorAnglesNode.at(i).get<double>());
                        }
                    } else {
                        for (const auto& item : rotorAnglesNode) {
                            if (!item.is_object()) {
                                throw std::runtime_error(
                                    "timeline rotor_angles array entries must be objects when not numeric");
                            }
                            std::size_t rotorIndex = 0;
                            if (item.contains("index")) {
                                rotorIndex = item.at("index").get<std::size_t>();
                            } else if (item.contains("idx")) {
                                rotorIndex = item.at("idx").get<std::size_t>();
                            } else if (item.contains("name")) {
                                const std::string name = item.at("name").get<std::string>();
                                const auto it = rotorNameToIndex.find(name);
                                if (it == rotorNameToIndex.end()) {
                                    throw std::runtime_error(
                                        "timeline rotor_angles references unknown rotor name: " + name);
                                }
                                rotorIndex = it->second;
                            } else {
                                throw std::runtime_error(
                                    "timeline rotor_angles entries must include 'index', 'idx', or 'name'");
                            }

                            double angleDeg = 0.0;
                            if (item.contains("angle_deg")) {
                                angleDeg = item.at("angle_deg").get<double>();
                            } else if (item.contains("angle")) {
                                angleDeg = item.at("angle").get<double>();
                            } else if (item.contains("degrees")) {
                                angleDeg = item.at("degrees").get<double>();
                            } else if (item.contains("value")) {
                                angleDeg = item.at("value").get<double>();
                            } else {
                                throw std::runtime_error(
                                    "timeline rotor_angles entries must include an angle field");
                            }
                            addRotorAngle(rotorIndex, angleDeg);
                        }
                    }
                } else if (rotorAnglesNode.is_object()) {
                    for (const auto& kv : rotorAnglesNode.items()) {
                        const auto it = rotorNameToIndex.find(kv.key());
                        if (it == rotorNameToIndex.end()) {
                            throw std::runtime_error("timeline rotor_angles references unknown rotor name: " +
                                                     kv.key());
                        }
                        if (!kv.value().is_number()) {
                            throw std::runtime_error(
                                "timeline rotor_angles object values must be numeric angles in degrees");
                        }
                        addRotorAngle(it->second, kv.value().get<double>());
                    }
                } else {
                    throw std::runtime_error("timeline rotor_angles must be an array or object");
                }
            }

            const auto parseWireOverrides = [&](const nlohmann::json& wiresNode) {
                if (!wiresNode.is_array()) {
                    throw std::runtime_error("timeline wires entry must be an array");
                }
                if (spec.wires.empty()) {
                    throw std::runtime_error(
                        "timeline wires overrides require the scenario to define wires");
                }
                if (wiresNode.empty()) {
                    return;
                }

                if (wiresNode.front().is_number()) {
                    if (wiresNode.size() != spec.wires.size()) {
                        throw std::runtime_error(
                            "timeline wire_currents array must match the number of wires in the scenario");
                    }
                    for (std::size_t i = 0; i < wiresNode.size(); ++i) {
                        ScenarioSpec::TimelineFrame::WireOverride override{};
                        override.index = i;
                        override.current = wiresNode.at(i).get<double>();
                        entry.wireOverrides.push_back(override);
                    }
                    return;
                }

                for (const auto& wireEntry : wiresNode) {
                    if (!wireEntry.is_object()) {
                        throw std::runtime_error(
                            "timeline wires overrides must be numeric array or objects with index/current");
                    }
                    ScenarioSpec::TimelineFrame::WireOverride override{};
                    override.index = wireEntry.at("index").get<std::size_t>();
                    if (override.index >= spec.wires.size()) {
                        throw std::runtime_error("timeline wire override index out of range");
                    }
                    if (wireEntry.contains("current")) {
                        override.current = wireEntry.at("current").get<double>();
                    } else if (wireEntry.contains("I")) {
                        override.current = wireEntry.at("I").get<double>();
                    } else if (wireEntry.contains("scale")) {
                        const double scale = wireEntry.at("scale").get<double>();
                        override.current = spec.wires[override.index].current * scale;
                    } else {
                        throw std::runtime_error(
                            "timeline wire override must provide 'current', 'I', or 'scale'");
                    }
                    entry.wireOverrides.push_back(override);
                }
            };

            if (frame.contains("wire_currents")) {
                parseWireOverrides(frame.at("wire_currents"));
            }
            if (frame.contains("wires")) {
                parseWireOverrides(frame.at("wires"));
            }

            const auto parseCurrentRegionOverrides = [&](const nlohmann::json& regionsNode) {
                if (!regionsNode.is_array()) {
                    throw std::runtime_error("timeline current_regions entry must be an array");
                }
                if (spec.currentRegions.empty()) {
                    throw std::runtime_error(
                        "timeline current region overrides require the scenario to define current_regions");
                }
                if (regionsNode.empty()) {
                    return;
                }
                if (regionsNode.front().is_number()) {
                    if (regionsNode.size() != spec.currentRegions.size()) {
                        throw std::runtime_error(
                            "timeline current region array must match number of current_regions in scenario");
                    }
                    for (std::size_t i = 0; i < regionsNode.size(); ++i) {
                        ScenarioSpec::TimelineFrame::CurrentRegionOverride override{};
                        override.index = i;
                        override.current = regionsNode.at(i).get<double>();
                        entry.currentRegionOverrides.push_back(override);
                    }
                    return;
                }

                for (const auto& regionEntry : regionsNode) {
                    if (!regionEntry.is_object()) {
                        throw std::runtime_error(
                            "timeline current region overrides must be numeric array or objects with index/id");
                    }
                    ScenarioSpec::TimelineFrame::CurrentRegionOverride override{};
                    if (regionEntry.contains("index")) {
                        override.index = regionEntry.at("index").get<std::size_t>();
                    } else if (regionEntry.contains("idx")) {
                        override.index = regionEntry.at("idx").get<std::size_t>();
                    } else if (regionEntry.contains("id")) {
                        const std::string regionId = regionEntry.at("id").get<std::string>();
                        const auto it = currentRegionIdToIndex.find(regionId);
                        if (it == currentRegionIdToIndex.end()) {
                            throw std::runtime_error(
                                "timeline current region override references unknown id: " + regionId);
                        }
                        override.index = it->second;
                    } else {
                        throw std::runtime_error(
                            "timeline current region overrides must specify 'index', 'idx', or 'id'");
                    }
                    if (override.index >= spec.currentRegions.size()) {
                        throw std::runtime_error("timeline current region override index out of range");
                    }
                    if (regionEntry.contains("current")) {
                        override.current = regionEntry.at("current").get<double>();
                    } else if (regionEntry.contains("I")) {
                        override.current = regionEntry.at("I").get<double>();
                    } else if (regionEntry.contains("scale")) {
                        const double scale = regionEntry.at("scale").get<double>();
                        override.current = spec.currentRegions[override.index].current * scale;
                    } else {
                        throw std::runtime_error(
                            "timeline current region override must provide 'current', 'I', or 'scale'");
                    }
                    entry.currentRegionOverrides.push_back(override);
                }
            };

            if (frame.contains("current_regions")) {
                parseCurrentRegionOverrides(frame.at("current_regions"));
            }
            if (frame.contains("slot_currents")) {
                parseCurrentRegionOverrides(frame.at("slot_currents"));
            }

            if (frame.contains("phase_currents")) {
                if (phaseToCurrentRegions.empty()) {
                    throw std::runtime_error(
                        "timeline phase_currents provided but no current regions declare a phase label");
                }
                const auto& phasesNode = frame.at("phase_currents");
                if (!phasesNode.is_object()) {
                    throw std::runtime_error("timeline phase_currents must be an object mapping phase name to current");
                }
                for (const auto& item : phasesNode.items()) {
                    const auto mapIt = phaseToCurrentRegions.find(item.key());
                    if (mapIt == phaseToCurrentRegions.end()) {
                        throw std::runtime_error("timeline phase_currents references unknown phase: " + item.key());
                    }
                    const double value = item.value().get<double>();
                    for (std::size_t index : mapIt->second) {
                        ScenarioSpec::TimelineFrame::CurrentRegionOverride override{};
                        override.index = index;
                        override.current = spec.currentRegions[index].orientation * value;
                        entry.currentRegionOverrides.push_back(override);
                    }
                }
            }

            const auto parseMagnetOverrides = [&](const nlohmann::json& magnetNode) {
                if (spec.magnetRegions.empty()) {
                    throw std::runtime_error(
                        "timeline magnet overrides require the scenario to define magnet_regions");
                }
                if (!magnetNode.is_array()) {
                    throw std::runtime_error("timeline magnets entry must be an array");
                }
                for (const auto& magnetEntry : magnetNode) {
                    if (!magnetEntry.is_object()) {
                        throw std::runtime_error(
                            "timeline magnets overrides must be objects with index/angle/magnetization");
                    }
                    ScenarioSpec::TimelineFrame::MagnetOverride override{};
                    override.index = magnetEntry.at("index").get<std::size_t>();
                    if (override.index >= spec.magnetRegions.size()) {
                        throw std::runtime_error("timeline magnet override index out of range");
                    }
                    if (magnetEntry.contains("angle_deg")) {
                        override.hasAngle = true;
                        override.angleDegrees = magnetEntry.at("angle_deg").get<double>();
                    } else if (magnetEntry.contains("angle")) {
                        override.hasAngle = true;
                        override.angleDegrees = magnetEntry.at("angle").get<double>();
                    }
                    if (magnetEntry.contains("magnetization")) {
                        const auto& mag = magnetEntry.at("magnetization");
                        if (mag.is_array()) {
                            if (mag.size() != 2) {
                                throw std::runtime_error(
                                    "timeline magnetization array must contain exactly two entries");
                            }
                            override.hasVector = true;
                            override.magnetizationX = mag.at(0).get<double>();
                            override.magnetizationY = mag.at(1).get<double>();
                        } else if (mag.is_object()) {
                            override.hasVector = true;
                            override.magnetizationX = mag.value("Mx", 0.0);
                            override.magnetizationY = mag.value("My", 0.0);
                        } else {
                            throw std::runtime_error(
                                "timeline magnetization must be an array or object with Mx/My");
                        }
                    }
                    entry.magnetOverrides.push_back(override);
                }
            };

            if (frame.contains("magnets")) {
                parseMagnetOverrides(frame.at("magnets"));
            }
            if (frame.contains("magnet_overrides")) {
                parseMagnetOverrides(frame.at("magnet_overrides"));
            }

            spec.timeline.push_back(std::move(entry));
            ++frameIndex;
        }
    }

    return spec;
}

void rasterizeScenarioToGrid(const ScenarioSpec& spec, Grid2D& grid) {
    if (grid.nx != spec.nx || grid.ny != spec.ny) {
        throw std::runtime_error("Grid dimensions do not match scenario domain");
    }

    grid.dx = spec.dx;
    grid.dy = spec.dy;
    grid.boundaryCondition = (spec.boundaryType == ScenarioSpec::BoundaryType::Dirichlet)
                                 ? Grid2D::BoundaryKind::Dirichlet
                                 : Grid2D::BoundaryKind::Neumann;

    const double invMuBackground = 1.0 / (MU0 * spec.mu_r_background);
    std::fill(grid.Jz.begin(), grid.Jz.end(), 0.0);
    const std::size_t cellCount = grid.nx * grid.ny;
    if (grid.Mx.size() != cellCount) {
        grid.Mx.assign(cellCount, 0.0);
    } else {
        std::fill(grid.Mx.begin(), grid.Mx.end(), 0.0);
    }
    if (grid.My.size() != cellCount) {
        grid.My.assign(cellCount, 0.0);
    } else {
        std::fill(grid.My.begin(), grid.My.end(), 0.0);
    }
    grid.Hx.clear();
    grid.Hy.clear();

    const double x0 = spec.originX;
    const double y0 = spec.originY;

    if (spec.regionMasks.empty()) {
        std::fill(grid.invMu.begin(), grid.invMu.end(), invMuBackground);
    } else {
        for (std::size_t j = 0; j < grid.ny; ++j) {
            const double y = y0 + static_cast<double>(j) * spec.dy;
            for (std::size_t i = 0; i < grid.nx; ++i) {
                const double x = x0 + static_cast<double>(i) * spec.dx;
                double invMu = invMuBackground;
                for (const auto& mask : spec.regionMasks) {
                    if (mask.kind == ScenarioSpec::RegionMask::Kind::Halfspace) {
                        const auto& hs = spec.halfspaces[mask.index];
                        if (hs.normal_x * x + hs.normal_y * y + hs.offset < 0.0) {
                            invMu = hs.inv_mu;
                        }
                    } else {
                        const auto& poly = spec.polygons[mask.index];
                        if (x < poly.min_x || x > poly.max_x || y < poly.min_y || y > poly.max_y) {
                            continue;
                        }
                        if (pointInPolygon(x, y, poly)) {
                            invMu = poly.inv_mu;
                        }
                    }
                }
                grid.invMu[grid.idx(i, j)] = invMu;
            }
        }
    }

    if (!spec.magnetRegions.empty()) {
        for (std::size_t j = 0; j < grid.ny; ++j) {
            const double y = y0 + static_cast<double>(j) * spec.dy;
            for (std::size_t i = 0; i < grid.nx; ++i) {
                const double x = x0 + static_cast<double>(i) * spec.dx;
                const std::size_t idx = grid.idx(i, j);
                for (const auto& region : spec.magnetRegions) {
                    const bool hasPolygon = region.xs.size() == region.ys.size() && region.xs.size() >= 3;
                    if (hasPolygon) {
                        if (x < region.min_x || x > region.max_x || y < region.min_y || y > region.max_y) {
                            continue;
                        }
                        if (pointInPolygonCoords(x, y, region.xs, region.ys)) {
                            grid.Mx[idx] += region.Mx;
                            grid.My[idx] += region.My;
                        }
                    } else {
                        if (pointInRect(x, y, region.min_x, region.max_x, region.min_y, region.max_y)) {
                            grid.Mx[idx] += region.Mx;
                            grid.My[idx] += region.My;
                        }
                    }
                }
            }
        }
    }

    if (!spec.currentRegions.empty()) {
        for (const auto& region : spec.currentRegions) {
            if (!(region.area > kTiny)) {
                continue;
            }
            const double density = region.current / region.area;
            if (std::abs(density) < kTiny) {
                continue;
            }
            const double minX = region.min_x - 1e-12;
            const double maxX = region.max_x + 1e-12;
            const double minY = region.min_y - 1e-12;
            const double maxY = region.max_y + 1e-12;
            for (std::size_t j = 0; j < grid.ny; ++j) {
                const double y = y0 + static_cast<double>(j) * spec.dy;
                if (y < minY || y > maxY) {
                    continue;
                }
                for (std::size_t i = 0; i < grid.nx; ++i) {
                    const double x = x0 + static_cast<double>(i) * spec.dx;
                    if (x < minX || x > maxX) {
                        continue;
                    }
                    if (pointInPolygonCoords(x, y, region.xs, region.ys)) {
                        grid.Jz[grid.idx(i, j)] += density;
                    }
                }
            }
        }
    }

    for (const auto& wire : spec.wires) {
        const double area = kPi * wire.radius * wire.radius;
        if (area <= kTiny) {
            continue;
        }
        const double jzDeposit = wire.current / area;
        for (std::size_t j = 0; j < grid.ny; ++j) {
            const double y = y0 + static_cast<double>(j) * spec.dy;
            for (std::size_t i = 0; i < grid.nx; ++i) {
                const double x = x0 + static_cast<double>(i) * spec.dx;
                const double r = std::hypot(x - wire.x, y - wire.y);
                if (r <= wire.radius) {
                    grid.Jz[grid.idx(i, j)] += jzDeposit;
                }
            }
        }
    }

    if (!spec.magnetRegions.empty()) {
        const double invDx = (grid.nx > 1) ? 1.0 / grid.dx : 0.0;
        const double invDy = (grid.ny > 1) ? 1.0 / grid.dy : 0.0;
        const double inv2Dx = (grid.nx > 2) ? 0.5 * invDx : invDx;
        const double inv2Dy = (grid.ny > 2) ? 0.5 * invDy : invDy;
        for (std::size_t j = 0; j < grid.ny; ++j) {
            for (std::size_t i = 0; i < grid.nx; ++i) {
                const std::size_t idx = grid.idx(i, j);

                double dMy_dx = 0.0;
                if (grid.nx > 1) {
                    if (i == 0) {
                        dMy_dx = (grid.My[grid.idx(i + 1, j)] - grid.My[idx]) * invDx;
                    } else if (i + 1 == grid.nx) {
                        dMy_dx = (grid.My[idx] - grid.My[grid.idx(i - 1, j)]) * invDx;
                    } else {
                        dMy_dx = (grid.My[grid.idx(i + 1, j)] - grid.My[grid.idx(i - 1, j)]) * inv2Dx;
                    }
                }

                double dMx_dy = 0.0;
                if (grid.ny > 1) {
                    if (j == 0) {
                        dMx_dy = (grid.Mx[grid.idx(i, j + 1)] - grid.Mx[idx]) * invDy;
                    } else if (j + 1 == grid.ny) {
                        dMx_dy = (grid.Mx[idx] - grid.Mx[grid.idx(i, j - 1)]) * invDy;
                    } else {
                        dMx_dy = (grid.Mx[grid.idx(i, j + 1)] - grid.Mx[grid.idx(i, j - 1)]) * inv2Dy;
                    }
                }

                grid.Jz[idx] += dMy_dx - dMx_dy;
            }
        }
    }
}

namespace {

void applyTimelineOverrides(const ScenarioSpec& baseSpec, const ScenarioSpec::TimelineFrame& frame,
                            ScenarioSpec& workingSpec) {
    if (!baseSpec.rotors.empty()) {
        std::vector<double> rotorAngles(baseSpec.rotors.size(), 0.0);
        std::vector<bool> rotorActive(baseSpec.rotors.size(), false);

        if (frame.hasRotorAngle && !baseSpec.rotors.empty()) {
            rotorAngles[0] = frame.rotorAngleDeg;
            rotorActive[0] = true;
        }

        for (const auto& override : frame.rotorAngles) {
            if (override.index >= baseSpec.rotors.size()) {
                throw std::runtime_error("timeline rotor angle index out of range");
            }
            rotorAngles[override.index] = override.angleDegrees;
            rotorActive[override.index] = true;
        }

        for (std::size_t r = 0; r < baseSpec.rotors.size(); ++r) {
            if (!rotorActive[r]) {
                continue;
            }
            const auto& rotor = baseSpec.rotors[r];
            const double angleRad = rotorAngles[r] * kPi / 180.0;
            const double cosA = std::cos(angleRad);
            const double sinA = std::sin(angleRad);

            for (std::size_t idx : rotor.polygonIndices) {
                if (idx >= baseSpec.polygons.size() || idx >= workingSpec.polygons.size()) {
                    throw std::runtime_error("rotor polygon index out of range");
                }
                const auto& basePoly = baseSpec.polygons[idx];
                auto& poly = workingSpec.polygons[idx];
                if (basePoly.xs.size() != basePoly.ys.size() || basePoly.xs.empty()) {
                    continue;
                }
                poly.xs.resize(basePoly.xs.size());
                poly.ys.resize(basePoly.ys.size());
                double minX = std::numeric_limits<double>::infinity();
                double maxX = -std::numeric_limits<double>::infinity();
                double minY = std::numeric_limits<double>::infinity();
                double maxY = -std::numeric_limits<double>::infinity();
                for (std::size_t v = 0; v < basePoly.xs.size(); ++v) {
                    const double baseX = basePoly.xs[v];
                    const double baseY = basePoly.ys[v];
                    const double dx = baseX - rotor.pivotX;
                    const double dy = baseY - rotor.pivotY;
                    const double rx = rotor.pivotX + dx * cosA - dy * sinA;
                    const double ry = rotor.pivotY + dx * sinA + dy * cosA;
                    poly.xs[v] = rx;
                    poly.ys[v] = ry;
                    minX = std::min(minX, rx);
                    maxX = std::max(maxX, rx);
                    minY = std::min(minY, ry);
                    maxY = std::max(maxY, ry);
                }
                poly.min_x = minX;
                poly.max_x = maxX;
                poly.min_y = minY;
                poly.max_y = maxY;
            }

            for (std::size_t idx : rotor.magnetIndices) {
                if (idx >= baseSpec.magnetRegions.size() || idx >= workingSpec.magnetRegions.size()) {
                    throw std::runtime_error("rotor magnet index out of range");
                }
                const auto& baseMagnet = baseSpec.magnetRegions[idx];
                auto& magnet = workingSpec.magnetRegions[idx];

                std::vector<double> baseXs = baseMagnet.xs;
                std::vector<double> baseYs = baseMagnet.ys;
                if (baseXs.size() != baseYs.size() || baseXs.size() < 3) {
                    baseXs = {baseMagnet.min_x, baseMagnet.max_x, baseMagnet.max_x, baseMagnet.min_x};
                    baseYs = {baseMagnet.min_y, baseMagnet.min_y, baseMagnet.max_y, baseMagnet.max_y};
                }
                magnet.xs.resize(baseXs.size());
                magnet.ys.resize(baseYs.size());
                double minX = std::numeric_limits<double>::infinity();
                double maxX = -std::numeric_limits<double>::infinity();
                double minY = std::numeric_limits<double>::infinity();
                double maxY = -std::numeric_limits<double>::infinity();
                for (std::size_t v = 0; v < baseXs.size(); ++v) {
                    const double dx = baseXs[v] - rotor.pivotX;
                    const double dy = baseYs[v] - rotor.pivotY;
                    const double rx = rotor.pivotX + dx * cosA - dy * sinA;
                    const double ry = rotor.pivotY + dx * sinA + dy * cosA;
                    magnet.xs[v] = rx;
                    magnet.ys[v] = ry;
                    minX = std::min(minX, rx);
                    maxX = std::max(maxX, rx);
                    minY = std::min(minY, ry);
                    maxY = std::max(maxY, ry);
                }
                magnet.min_x = minX;
                magnet.max_x = maxX;
                magnet.min_y = minY;
                magnet.max_y = maxY;

                double mx = baseMagnet.Mx;
                double my = baseMagnet.My;
                magnet.Mx = mx * cosA - my * sinA;
                magnet.My = mx * sinA + my * cosA;
            }

            for (std::size_t idx : rotor.wireIndices) {
                if (idx >= baseSpec.wires.size() || idx >= workingSpec.wires.size()) {
                    throw std::runtime_error("rotor wire index out of range");
                }
                const auto& baseWire = baseSpec.wires[idx];
                auto& wire = workingSpec.wires[idx];
                const double dx = baseWire.x - rotor.pivotX;
                const double dy = baseWire.y - rotor.pivotY;
                wire.x = rotor.pivotX + dx * cosA - dy * sinA;
                wire.y = rotor.pivotY + dx * sinA + dy * cosA;
            }
        }
    } else if (frame.hasRotorAngle) {
        if (baseSpec.magnetRegions.empty()) {
            throw std::runtime_error(
                "timeline rotor_angle specified but scenario defines no magnet_regions");
        }
        const double angleRad = frame.rotorAngleDeg * kPi / 180.0;
        const double cosA = std::cos(angleRad);
        const double sinA = std::sin(angleRad);
        for (std::size_t i = 0; i < workingSpec.magnetRegions.size(); ++i) {
            const auto& baseRegion = baseSpec.magnetRegions[i];
            double mx = baseRegion.Mx;
            double my = baseRegion.My;
            workingSpec.magnetRegions[i].Mx = mx * cosA - my * sinA;
            workingSpec.magnetRegions[i].My = mx * sinA + my * cosA;
        }
    }

    for (const auto& override : frame.magnetOverrides) {
        if (override.index >= workingSpec.magnetRegions.size()) {
            throw std::runtime_error("timeline magnet override index out of range");
        }
        auto& region = workingSpec.magnetRegions[override.index];
        const auto& baseRegion = baseSpec.magnetRegions[override.index];
        if (override.hasVector) {
            region.Mx = override.magnetizationX;
            region.My = override.magnetizationY;
        }
        if (override.hasAngle) {
            const double angleRad = override.angleDegrees * kPi / 180.0;
            const double cosA = std::cos(angleRad);
            const double sinA = std::sin(angleRad);
            double mx = baseRegion.Mx;
            double my = baseRegion.My;
            region.Mx = mx * cosA - my * sinA;
            region.My = mx * sinA + my * cosA;
        }
    }

    for (const auto& override : frame.wireOverrides) {
        if (override.index >= workingSpec.wires.size()) {
            throw std::runtime_error("timeline wire override index out of range");
        }
        workingSpec.wires[override.index].current = override.current;
    }

    for (const auto& override : frame.currentRegionOverrides) {
        if (override.index >= workingSpec.currentRegions.size()) {
            throw std::runtime_error("timeline current region override index out of range");
        }
        workingSpec.currentRegions[override.index].current = override.current;
    }
}

}  // namespace

std::vector<ScenarioFrame> expandScenarioTimeline(const ScenarioSpec& spec) {
    std::vector<ScenarioFrame> frames;
    if (spec.timeline.empty()) {
        ScenarioFrame frame{};
        frame.index = 0;
        frame.time = 0.0;
        frame.spec = spec;
        frame.spec.timeline.clear();
        frames.push_back(std::move(frame));
        return frames;
    }

    frames.reserve(spec.timeline.size());
    for (std::size_t idx = 0; idx < spec.timeline.size(); ++idx) {
        ScenarioFrame frame{};
        frame.index = idx;
        frame.time = spec.timeline[idx].time;
        frame.spec = spec;
        frame.spec.timeline.clear();
        applyTimelineOverrides(spec, spec.timeline[idx], frame.spec);
        frames.push_back(std::move(frame));
    }

    return frames;
}

}  // namespace motorsim
