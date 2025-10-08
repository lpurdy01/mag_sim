#include "motorsim/ingest.hpp"

#include "motorsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <nlohmann/json.hpp>

namespace motorsim {
namespace {
constexpr double kTiny = 1e-12;
constexpr double kPi = 3.14159265358979323846;

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
    if (spec.version != "0.1") {
        throw std::runtime_error("Unsupported scenario version: " + spec.version);
    }

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

    const auto& materials = json.at("materials");
    if (!materials.is_array() || materials.empty()) {
        throw std::runtime_error("Scenario must define at least one material");
    }

    std::unordered_map<std::string, double> muMap;
    for (const auto& material : materials) {
        const std::string name = material.at("name").get<std::string>();
        if (muMap.count(name) != 0U) {
            throw std::runtime_error("Duplicate material name: " + name);
        }
        const double mu_r = requirePositive("material.mu_r", material.at("mu_r").get<double>());
        muMap.emplace(name, mu_r);
    }

    const auto& regions = json.at("regions");
    if (!regions.is_array() || regions.empty()) {
        throw std::runtime_error("Scenario must define at least one region");
    }

    bool foundUniform = false;
    for (const auto& region : regions) {
        const std::string type = region.at("type").get<std::string>();
        if (type != "uniform") {
            throw std::runtime_error("Unsupported region type: " + type);
        }
        const std::string materialName = region.at("material").get<std::string>();
        const auto it = muMap.find(materialName);
        if (it == muMap.end()) {
            throw std::runtime_error("Region references unknown material: " + materialName);
        }
        if (foundUniform) {
            throw std::runtime_error("Multiple uniform regions defined; only one is supported in v0.1");
        }
        spec.mu_r_background = it->second;
        foundUniform = true;
    }

    if (!foundUniform) {
        throw std::runtime_error("Scenario must define a uniform region in v0.1");
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

    return spec;
}

void rasterizeScenarioToGrid(const ScenarioSpec& spec, Grid2D& grid) {
    if (grid.nx != spec.nx || grid.ny != spec.ny) {
        throw std::runtime_error("Grid dimensions do not match scenario domain");
    }

    grid.dx = spec.dx;
    grid.dy = spec.dy;

    const double invMuBackground = 1.0 / (MU0 * spec.mu_r_background);
    std::fill(grid.invMu.begin(), grid.invMu.end(), invMuBackground);
    std::fill(grid.Jz.begin(), grid.Jz.end(), 0.0);

    const double x0 = spec.originX;
    const double y0 = spec.originY;

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
}

}  // namespace motorsim
