#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"
#include "motorsim/probes.hpp"
#include "motorsim/solver.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

namespace {

struct SolveOutcome {
    motorsim::Grid2D grid;
    bool converged{false};
    double relResidual{0.0};
};

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

    std::vector<double> torques;
    torques.reserve(frames.size());

    for (const ScenarioFrame& frame : frames) {
        SolveOutcome outcome = solve_scenario(frame.spec);
        if (!outcome.converged || !(outcome.relResidual < 2.5e-6)) {
            std::cerr << "Frame " << frame.index << " failed to converge (relResidual=" << outcome.relResidual
                      << ")" << '\n';
            return 1;
        }

        StressTensorResult stress = evaluate_maxwell_stress_probe(outcome.grid, frame.spec.originX,
                                                                  frame.spec.originY, frame.spec.dx, frame.spec.dy,
                                                                  loopXs, loopYs);
        torques.push_back(stress.torqueZ);
    }

    const auto magnitude = [](double value) { return std::abs(value); };

    if (!(magnitude(torques[0]) > 5.0)) {
        std::cerr << "Torque at 0 degrees too small: " << torques[0] << '\n';
        return 1;
    }
    if (!(magnitude(torques[3]) > 5.0)) {
        std::cerr << "Torque at 180 degrees too small: " << torques[3] << '\n';
        return 1;
    }
    if (!(torques[1] * torques[2] < 0.0)) {
        std::cerr << "Frames 60 and 120 degrees should have opposing torque signs: " << torques[1] << ", "
                  << torques[2] << '\n';
        return 1;
    }
    if (!(magnitude(torques[1]) > 40.0 && magnitude(torques[2]) > 40.0)) {
        std::cerr << "Peak torque magnitudes too small: " << torques[1] << ", " << torques[2] << '\n';
        return 1;
    }
    const double dominant = std::max(magnitude(torques[1]), magnitude(torques[2]));
    if (!(dominant > magnitude(torques[0]) * 2.0 && dominant > magnitude(torques[3]) * 2.0)) {
        std::cerr << "Peak torque does not dominate baseline frames" << '\n';
        return 1;
    }

    const auto minmax = std::minmax_element(torques.begin(), torques.end());
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

    const double poleOppositionError = oppositionError(torques[1], torques[2]);
    const double repeatError0 = repeatError(torques[0], torques[3]);

    std::cout << "RotorRipple: torque_deg0=" << torques[0] << " torque_deg60=" << torques[1]
              << " torque_deg120=" << torques[2] << " torque_deg180=" << torques[3]
              << " peak_to_peak=" << peakToPeak << " opposition60_120=" << poleOppositionError
              << " repeat0_180=" << repeatError0 << '\n';

    return 0;
}
