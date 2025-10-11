#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"
#include "motorsim/probes.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;

struct SolveOutcome {
    motorsim::Grid2D grid;
    bool converged{false};
    double relResidual{0.0};
};

motorsim::ScenarioSpec rotate_magnet(const motorsim::ScenarioSpec& base, double angleDeg) {
    motorsim::ScenarioSpec rotated = base;
    const double angleRad = angleDeg * kPi / 180.0;
    const double cosA = std::cos(angleRad);
    const double sinA = std::sin(angleRad);
    for (std::size_t i = 0; i < rotated.magnetRegions.size(); ++i) {
        const double mx = base.magnetRegions[i].Mx;
        const double my = base.magnetRegions[i].My;
        rotated.magnetRegions[i].Mx = mx * cosA - my * sinA;
        rotated.magnetRegions[i].My = mx * sinA + my * cosA;
    }
    return rotated;
}

SolveOutcome solve_scenario(const motorsim::ScenarioSpec& spec) {
    motorsim::Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    rasterizeScenarioToGrid(spec, grid);

    motorsim::SolveOptions options{};
    options.maxIters = 40000;
    options.tol = 2.5e-6;
    options.omega = 1.7;
    options.verbose = false;

    const motorsim::SolveReport report = solveAz_GS_SOR(grid, options);
    SolveOutcome outcome{std::move(grid), report.converged, report.relResidual};
    if (!report.converged || !(report.relResidual < options.tol)) {
        return outcome;
    }
    computeB(outcome.grid);
    computeH(outcome.grid);
    return outcome;
}

double compute_interaction_energy(const motorsim::Grid2D& grid, const motorsim::ScenarioSpec& spec) {
    const std::size_t count = grid.nx * grid.ny;
    const double cellArea = spec.dx * spec.dy;
    double energy = 0.0;
    for (std::size_t idx = 0; idx < count; ++idx) {
        const double mx = (idx < grid.Mx.size()) ? grid.Mx[idx] : 0.0;
        const double my = (idx < grid.My.size()) ? grid.My[idx] : 0.0;
        if (mx == 0.0 && my == 0.0) {
            continue;
        }
        energy -= (mx * grid.Bx[idx] + my * grid.By[idx]) * cellArea;
    }
    return energy;
}

}  // namespace

int main() {
    using namespace motorsim;
    namespace fs = std::filesystem;

    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/torque_validation_test.json")
            .lexically_normal();

    ScenarioSpec baseSpec;
    try {
        baseSpec = loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load torque validation scenario: " << ex.what() << '\n';
        return 1;
    }

    if (baseSpec.magnetRegions.size() != 1) {
        std::cerr << "Expected a single magnet region" << '\n';
        return 1;
    }
    const ScenarioSpec::MagnetRegion magnet = baseSpec.magnetRegions.front();
    const double magnetArea = (magnet.max_x - magnet.min_x) * (magnet.max_y - magnet.min_y);
    if (!(magnetArea > 0.0)) {
        std::cerr << "Magnet area is non-positive" << '\n';
        return 1;
    }

    SolveOutcome baseOutcome = solve_scenario(baseSpec);
    if (!baseOutcome.converged || !(baseOutcome.relResidual < 2.5e-6)) {
        std::cerr << "Base frame failed to converge (relResidual=" << baseOutcome.relResidual << ")" << '\n';
        return 1;
    }

    const double margin = 0.004;
    const double loopMinX = magnet.min_x - margin;
    const double loopMaxX = magnet.max_x + margin;
    const double loopMinY = magnet.min_y - margin;
    const double loopMaxY = magnet.max_y + margin;

    std::vector<double> loopXs{loopMinX, loopMaxX, loopMaxX, loopMinX};
    std::vector<double> loopYs{loopMinY, loopMinY, loopMaxY, loopMaxY};

    StressTensorResult stress;
    try {
        stress = evaluate_maxwell_stress_probe(baseOutcome.grid, baseSpec.originX, baseSpec.originY,
                                               baseSpec.dx, baseSpec.dy, loopXs, loopYs);
    } catch (const std::exception& ex) {
        std::cerr << "Failed to evaluate Maxwell stress probe: " << ex.what() << '\n';
        return 1;
    }

    double sumBy = 0.0;
    std::size_t magnetSamples = 0;
    for (std::size_t j = 0; j < baseSpec.ny; ++j) {
        const double y = baseSpec.originY + static_cast<double>(j) * baseSpec.dy;
        if (y < magnet.min_y || y > magnet.max_y) {
            continue;
        }
        for (std::size_t i = 0; i < baseSpec.nx; ++i) {
            const double x = baseSpec.originX + static_cast<double>(i) * baseSpec.dx;
            if (x < magnet.min_x || x > magnet.max_x) {
                continue;
            }
            const std::size_t idx = baseOutcome.grid.idx(i, j);
            sumBy += baseOutcome.grid.By[idx];
            ++magnetSamples;
        }
    }
    if (magnetSamples == 0) {
        std::cerr << "No grid samples found inside magnet" << '\n';
        return 1;
    }

    const double avgBy = sumBy / static_cast<double>(magnetSamples);
    const double dipoleTorque = magnetArea * magnet.Mx * avgBy;
    if (!(dipoleTorque > 0.0)) {
        std::cerr << "Expected positive dipole torque estimate, got " << dipoleTorque << '\n';
        return 1;
    }

    const double mstTorque = stress.torqueZ;
    if (!(mstTorque > 0.0)) {
        std::cerr << "Expected positive Maxwell stress torque, got " << mstTorque << '\n';
        return 1;
    }

    const double relDiffDipole = std::abs(mstTorque - dipoleTorque) / dipoleTorque;
    if (!(relDiffDipole < 0.2)) {
        std::cerr << "MST torque " << mstTorque << " differs from dipole estimate " << dipoleTorque
                  << " by relDiff=" << relDiffDipole << '\n';
        return 1;
    }

    const double deltaAngleDeg = 5.0;
    const ScenarioSpec specMinus = rotate_magnet(baseSpec, -deltaAngleDeg);
    const ScenarioSpec specPlus = rotate_magnet(baseSpec, deltaAngleDeg);

    SolveOutcome minusOutcome = solve_scenario(specMinus);
    SolveOutcome plusOutcome = solve_scenario(specPlus);
    if (!minusOutcome.converged || !(minusOutcome.relResidual < 2.5e-6) || !plusOutcome.converged ||
        !(plusOutcome.relResidual < 2.5e-6)) {
        std::cerr << "Virtual-work frames failed to converge" << '\n';
        return 1;
    }

    const double energyMinus = compute_interaction_energy(minusOutcome.grid, specMinus);
    const double energyPlus = compute_interaction_energy(plusOutcome.grid, specPlus);
    const double deltaAngleRad = deltaAngleDeg * kPi / 180.0;
    const double virtualTorque = -(energyPlus - energyMinus) / (2.0 * deltaAngleRad);

    const double relDiffVirtual = std::abs(mstTorque - virtualTorque) / std::max(1e-9, std::abs(virtualTorque));
    if (!(relDiffVirtual < 0.25)) {
        std::cerr << "Virtual-work torque " << virtualTorque << " differs from MST torque " << mstTorque
                  << " by relDiff=" << relDiffVirtual << '\n';
        return 1;
    }

    std::cout << "TorqueValidation: mst=" << mstTorque << " dipole=" << dipoleTorque
              << " relDiffDipole=" << relDiffDipole << " virtual=" << virtualTorque
              << " relDiffVirtual=" << relDiffVirtual << '\n';

    return 0;
}
