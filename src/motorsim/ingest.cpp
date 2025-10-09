#include "motorsim/ingest.hpp"

#include "motorsim/types.hpp"

#include <algorithm>
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
    spec.wires.reserve(sources.size());
    for (const auto& source : sources) {
        const std::string type = source.at("type").get<std::string>();
        if (type != "wire") {
            throw std::runtime_error("Unsupported source type: " + type);
        }
        ScenarioSpec::Wire wire{};
        wire.x = source.at("x").get<double>();
        wire.y = source.at("y").get<double>();
        wire.radius = requirePositive("wire.radius", source.at("radius").get<double>());
        wire.current = source.at("I").get<double>();
        spec.wires.push_back(wire);
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
            } else {
                throw std::runtime_error("Unsupported magnet region type: " + shape);
            }

            spec.magnetRegions.push_back(std::move(region));
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
                if (request.format != "csv") {
                    throw std::runtime_error(
                        "Unsupported field_map format for output '" + id + "': " + request.format);
                }
                if (output.contains("path")) {
                    request.path =
                        requireNonEmpty("outputs.path", output.at("path").get<std::string>());
                } else {
                    request.path = "outputs/" + id + ".csv";
                }
                spec.outputs.fieldMaps.push_back(std::move(request));
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
            } else {
                throw std::runtime_error("Unsupported output type: " + type);
            }
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
                    if (region.shape == ScenarioSpec::MagnetRegion::Shape::Polygon) {
                        if (x < region.min_x || x > region.max_x || y < region.min_y || y > region.max_y) {
                            continue;
                        }
                        if (pointInRect(x, y, region.min_x, region.max_x, region.min_y, region.max_y) &&
                            pointInPolygonCoords(x, y, region.xs, region.ys)) {
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

}  // namespace motorsim
