#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"
#include "motorsim/probes.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846;

struct SolveOutcome {
    motorsim::Grid2D grid;
    motorsim::SolveResult result{};
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

SolveOutcome solve_scenario(const motorsim::ScenarioSpec& spec, motorsim::SolverKind solverKind) {
    motorsim::Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    rasterizeScenarioToGrid(spec, grid);

    motorsim::SolveOptions options{};
    options.maxIters = 40000;
    options.tol = 2.5e-6;
    options.omega = 1.7;
    options.verbose = false;

    options.kind = solverKind;

    SolveOutcome outcome{std::move(grid), {}};
    outcome.result = solveAz(outcome.grid, options);
    if (outcome.result.converged && outcome.result.relResidual < options.tol) {
        computeB(outcome.grid);
        computeH(outcome.grid);
    }
    return outcome;
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

    const double margin = 0.004;
    const double loopMinX = magnet.min_x - margin;
    const double loopMaxX = magnet.max_x + margin;
    const double loopMinY = magnet.min_y - margin;
    const double loopMaxY = magnet.max_y + margin;

    std::vector<double> loopXs{loopMinX, loopMaxX, loopMaxX, loopMinX};
    std::vector<double> loopYs{loopMinY, loopMinY, loopMaxY, loopMaxY};

    struct SolverMetrics {
        SolveOutcome base;
        SolveOutcome minus;
        SolveOutcome plus;
        double mstTorque{0.0};
        double dipoleTorque{0.0};
        double relDiffDipole{0.0};
        double virtualTorque{0.0};
        double relDiffVirtual{0.0};
        double coenergyMinus{0.0};
        double coenergyPlus{0.0};
    };

    const std::array<std::pair<SolverKind, const char*>, 2> solverCases = {
        std::make_pair(SolverKind::SOR, "SOR"),
        std::make_pair(SolverKind::CG, "CG")};

    std::array<SolverMetrics, 2> metrics{};

    const double deltaAngleDeg = 2.0;
    const ScenarioSpec specMinus = rotate_magnet(baseSpec, -deltaAngleDeg);
    const ScenarioSpec specPlus = rotate_magnet(baseSpec, deltaAngleDeg);
    const double deltaAngleRad = deltaAngleDeg * kPi / 180.0;

    for (std::size_t idx = 0; idx < solverCases.size(); ++idx) {
        const auto [solverKind, label] = solverCases[idx];
        SolverMetrics& m = metrics[idx];

        m.base = solve_scenario(baseSpec, solverKind);
        if (!m.base.result.converged || !(m.base.result.relResidual < 2.5e-6)) {
            std::cerr << label << " base frame failed to converge (relResidual=" << m.base.result.relResidual
                      << ")" << '\n';
            return 1;
        }

        StressTensorResult stress;
        try {
            stress = evaluate_maxwell_stress_probe(m.base.grid, baseSpec.originX, baseSpec.originY,
                                                   baseSpec.dx, baseSpec.dy, loopXs, loopYs);
        } catch (const std::exception& ex) {
            std::cerr << label << ": failed to evaluate Maxwell stress probe: " << ex.what() << '\n';
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
                const std::size_t gidx = m.base.grid.idx(i, j);
                sumBy += m.base.grid.By[gidx];
                ++magnetSamples;
            }
        }
        if (magnetSamples == 0) {
            std::cerr << label << ": no grid samples found inside magnet" << '\n';
            return 1;
        }

        const double avgBy = sumBy / static_cast<double>(magnetSamples);
        m.dipoleTorque = magnetArea * magnet.Mx * avgBy;
        if (!(m.dipoleTorque > 0.0)) {
            std::cerr << label << ": expected positive dipole torque estimate, got " << m.dipoleTorque << '\n';
            return 1;
        }

        m.mstTorque = stress.torqueZ;
        if (!(m.mstTorque > 0.0)) {
            std::cerr << label << ": expected positive Maxwell stress torque, got " << m.mstTorque << '\n';
            return 1;
        }

        m.relDiffDipole = std::abs(m.mstTorque - m.dipoleTorque) / m.dipoleTorque;
        if (!(m.relDiffDipole < 0.2)) {
            std::cerr << label << ": MST torque " << m.mstTorque << " differs from dipole estimate "
                      << m.dipoleTorque << " by relDiff=" << m.relDiffDipole << '\n';
            return 1;
        }

        m.minus = solve_scenario(specMinus, solverKind);
        m.plus = solve_scenario(specPlus, solverKind);
        if (!m.minus.result.converged || !(m.minus.result.relResidual < 2.5e-6) || !m.plus.result.converged ||
            !(m.plus.result.relResidual < 2.5e-6)) {
            std::cerr << label << ": virtual-work frames failed to converge" << '\n';
            return 1;
        }

        try {
            m.coenergyMinus = motorsim::compute_magnetic_coenergy(m.minus.grid, specMinus.dx, specMinus.dy);
            m.coenergyPlus = motorsim::compute_magnetic_coenergy(m.plus.grid, specPlus.dx, specPlus.dy);
        } catch (const std::exception& ex) {
            std::cerr << label << ": failed to compute magnetic co-energy: " << ex.what() << '\n';
            return 1;
        }

        m.virtualTorque = (m.coenergyPlus - m.coenergyMinus) / (2.0 * deltaAngleRad);

        m.relDiffVirtual = std::abs(m.mstTorque - m.virtualTorque) /
                           std::max(1e-9, std::abs(m.virtualTorque));
        if (!(m.relDiffVirtual < 0.10)) {
            std::cerr << label << ": virtual-work torque " << m.virtualTorque
                      << " differs from MST torque " << m.mstTorque
                      << " by relDiff=" << m.relDiffVirtual << '\n';
            return 1;
        }
    }

    const SolverMetrics& sorMetrics = metrics[0];
    const SolverMetrics& cgMetrics = metrics[1];

    const double torqueDenom = std::max({std::abs(sorMetrics.mstTorque), 1.0, std::abs(cgMetrics.mstTorque)});
    const double torqueRelDiff = std::abs(sorMetrics.mstTorque - cgMetrics.mstTorque) / torqueDenom;
    if (torqueRelDiff > 0.02) {
        std::cerr << "CG MST torque deviates from SOR (relDiff=" << torqueRelDiff << ")" << '\n';
        return 1;
    }

    const double dipoleDenom = std::max({std::abs(sorMetrics.dipoleTorque), 1.0, std::abs(cgMetrics.dipoleTorque)});
    if (std::abs(sorMetrics.dipoleTorque - cgMetrics.dipoleTorque) / dipoleDenom > 0.02) {
        std::cerr << "CG dipole torque deviates from SOR" << '\n';
        return 1;
    }

    if (std::abs(sorMetrics.virtualTorque - cgMetrics.virtualTorque) /
            std::max({std::abs(sorMetrics.virtualTorque), 1.0, std::abs(cgMetrics.virtualTorque)}) > 0.05) {
        std::cerr << "CG virtual torque deviates from SOR" << '\n';
        return 1;
    }

    const double sorResidualBound = sorMetrics.base.result.relResidual * 1.1 + 1e-12;
    if (cgMetrics.base.result.relResidual > sorResidualBound) {
        std::cerr << "CG residual exceeds SOR baseline" << '\n';
        return 1;
    }

    std::cout << "TorqueValidation: mst=" << sorMetrics.mstTorque << " dipole=" << sorMetrics.dipoleTorque
              << " relDiffDipole=" << sorMetrics.relDiffDipole << " virtual=" << sorMetrics.virtualTorque
              << " relDiffVirtual=" << sorMetrics.relDiffVirtual << '\n';

    return 0;
}
