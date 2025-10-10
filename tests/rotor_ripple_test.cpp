#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"
#include "motorsim/probes.hpp"
#include "motorsim/solver.hpp"

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

    if (spec.magnetRegions.size() != 1) {
        std::cerr << "Scenario requires one magnet region" << '\n';
        return 1;
    }
    const ScenarioSpec::MagnetRegion magnet = spec.magnetRegions.front();

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

    const double margin = 0.004;
    std::vector<double> loopXs{magnet.min_x - margin, magnet.max_x + margin, magnet.max_x + margin,
                               magnet.min_x - margin};
    std::vector<double> loopYs{magnet.min_y - margin, magnet.min_y - margin, magnet.max_y + margin,
                               magnet.max_y + margin};

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

    if (!(torques[0] > 0.05)) {
        std::cerr << "Expected positive torque at 0 degrees, got " << torques[0] << '\n';
        return 1;
    }
    if (!(torques[1] > 0.0 && torques[1] < torques[0])) {
        std::cerr << "Expected reduced positive torque at 60 degrees, got " << torques[1] << '\n';
        return 1;
    }
    if (!(torques[2] < -0.02)) {
        std::cerr << "Expected negative torque at 120 degrees, got " << torques[2] << '\n';
        return 1;
    }
    if (!(torques[3] < torques[2])) {
        std::cerr << "Expected strongest negative torque at 180 degrees" << '\n';
        return 1;
    }

    const double peakToPeak = torques[0] - torques[3];
    if (!(peakToPeak > 0.3)) {
        std::cerr << "Torque ripple too small: peak-to-peak=" << peakToPeak << '\n';
        return 1;
    }

    std::cout << "RotorRipple: torque_deg0=" << torques[0] << " torque_deg60=" << torques[1]
              << " torque_deg120=" << torques[2] << " torque_deg180=" << torques[3]
              << " peak_to_peak=" << peakToPeak << '\n';

    return 0;
}
