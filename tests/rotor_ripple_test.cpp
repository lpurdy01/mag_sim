#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"
#include "motorsim/probes.hpp"
#include "motorsim/solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

namespace {

struct SolveOutcome {
    motorsim::Grid2D grid;
    motorsim::SolveResult result{};
};

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
        (fs::path(__FILE__).parent_path() / "../inputs/tests/rotor_ripple_test.json").lexically_normal();

    ScenarioSpec spec;
    try {
        spec = loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load rotor ripple scenario: " << ex.what() << '\n';
        return 1;
    }

    const ScenarioSpec::Outputs::Probe* torqueProbe = nullptr;
    for (const auto& probe : spec.outputs.probes) {
        if (probe.method == ScenarioSpec::Outputs::Probe::Method::StressTensor &&
            probe.quantity != ScenarioSpec::Outputs::Probe::Quantity::Force) {
            torqueProbe = &probe;
            if (probe.id == "torque_loop") {
                break;
            }
        }
    }
    if (torqueProbe == nullptr) {
        std::cerr << "Scenario must define a stress-tensor torque probe" << '\n';
        return 1;
    }
    if (torqueProbe->loopXs.size() != torqueProbe->loopYs.size() || torqueProbe->loopXs.size() < 3) {
        std::cerr << "Torque probe loop is malformed" << '\n';
        return 1;
    }

    std::vector<ScenarioFrame> frames;
    try {
        frames = expandScenarioTimeline(spec);
    } catch (const std::exception& ex) {
        std::cerr << "Timeline expansion failed: " << ex.what() << '\n';
        return 1;
    }

    if (frames.size() != 4) {
        std::cerr << "Expected four frames, got " << frames.size() << '\n';
        return 1;
    }

    const std::vector<double> loopXs = torqueProbe->loopXs;
    const std::vector<double> loopYs = torqueProbe->loopYs;

    const std::array<std::pair<motorsim::SolverKind, const char*>, 2> solverCases = {
        std::make_pair(motorsim::SolverKind::SOR, "SOR"),
        std::make_pair(motorsim::SolverKind::CG, "CG")};

    std::array<std::vector<double>, 2> solverTorques{};
    std::array<std::vector<motorsim::SolveResult>, 2> solverReports{};

    for (std::size_t solverIdx = 0; solverIdx < solverCases.size(); ++solverIdx) {
        const auto [solverKind, label] = solverCases[solverIdx];
        std::vector<double>& torques = solverTorques[solverIdx];
        std::vector<motorsim::SolveResult>& reports = solverReports[solverIdx];
        torques.clear();
        reports.clear();
        torques.reserve(frames.size());
        reports.reserve(frames.size());

        for (const ScenarioFrame& frame : frames) {
            SolveOutcome outcome = solve_scenario(frame.spec, solverKind);
            reports.push_back(outcome.result);
            if (!outcome.result.converged || !(outcome.result.relResidual < 2.5e-6)) {
                std::cerr << label << " frame " << frame.index
                          << " failed to converge (relResidual=" << outcome.result.relResidual << ")" << '\n';
                return 1;
            }

            StressTensorResult stress = evaluate_maxwell_stress_probe(outcome.grid, frame.spec.originX,
                                                                      frame.spec.originY, frame.spec.dx,
                                                                      frame.spec.dy, loopXs, loopYs);
            torques.push_back(stress.torqueZ);
        }
    }

    const std::vector<double>& sorTorques = solverTorques[0];
    const std::vector<double>& cgTorques = solverTorques[1];

    const auto magnitude = [](double value) { return std::abs(value); };

    if (!(magnitude(sorTorques[0]) < 5.0)) {
        std::cerr << "Torque at 0 degrees should be near zero: " << sorTorques[0] << '\n';
        return 1;
    }
    if (!(magnitude(sorTorques[3]) < 5.0)) {
        std::cerr << "Torque at 180 degrees should be near zero: " << sorTorques[3] << '\n';
        return 1;
    }
    if (!(sorTorques[1] * sorTorques[2] < 0.0)) {
        std::cerr << "Frames 60 and 120 degrees should have opposing torque signs: " << sorTorques[1] << ", "
                  << sorTorques[2] << '\n';
        return 1;
    }
    if (!(magnitude(sorTorques[1]) > 50.0 && magnitude(sorTorques[2]) > 50.0)) {
        std::cerr << "Peak torque magnitudes too small: " << sorTorques[1] << ", " << sorTorques[2] << '\n';
        return 1;
    }
    const double dominant = std::max(magnitude(sorTorques[1]), magnitude(sorTorques[2]));
    if (!(dominant > magnitude(sorTorques[0]) * 2.0 && dominant > magnitude(sorTorques[3]) * 2.0)) {
        std::cerr << "Peak torque does not dominate baseline frames" << '\n';
        return 1;
    }

    const auto minmax = std::minmax_element(sorTorques.begin(), sorTorques.end());
    const double peakToPeak = *minmax.second - *minmax.first;
    if (!(peakToPeak > 120.0)) {
        std::cerr << "Torque ripple too small: peak-to-peak=" << peakToPeak << '\n';
        return 1;
    }

    const auto oppositionError = [](double a, double b) {
        const double denom = std::max(std::max(std::abs(a), std::abs(b)), 1e-9);
        return std::abs(a + b) / denom;
    };
    const auto repeatError = [](double a, double b) {
        const double denom = std::max(std::max(std::abs(a), std::abs(b)), 1e-9);
        return std::abs(a - b) / denom;
    };

    const double poleOppositionError = oppositionError(sorTorques[1], sorTorques[2]);
    const double repeatError0 = repeatError(sorTorques[0], sorTorques[3]);

    std::cout << "RotorRipple: torque_deg0=" << sorTorques[0] << " torque_deg60=" << sorTorques[1]
              << " torque_deg120=" << sorTorques[2] << " torque_deg180=" << sorTorques[3]
              << " peak_to_peak=" << peakToPeak << " opposition60_120=" << poleOppositionError
              << " repeat0_180=" << repeatError0 << '\n';

    for (std::size_t frameIdx = 0; frameIdx < frames.size(); ++frameIdx) {
        const double sorTorque = sorTorques[frameIdx];
        const double cgTorque = cgTorques[frameIdx];
        const double denom = std::max({std::abs(sorTorque), 1.0, std::abs(cgTorque)});
        const double relDiff = std::abs(sorTorque - cgTorque) / denom;
        if (relDiff > 0.02) {
            std::cerr << "CG torque deviates from SOR on frame " << frames[frameIdx].index
                      << " (relDiff=" << relDiff << ")" << '\n';
            return 1;
        }

        const motorsim::SolveResult& sorReport = solverReports[0][frameIdx];
        const motorsim::SolveResult& cgReport = solverReports[1][frameIdx];
        const double sorBound = sorReport.relResidual * 1.1 + 1e-12;
        if (cgReport.relResidual > sorBound) {
            std::cerr << "CG residual exceeds SOR baseline on frame " << frames[frameIdx].index
                      << " (cg=" << cgReport.relResidual << ", bound=" << sorBound << ")" << '\n';
            return 1;
        }
    }

    return 0;
}
