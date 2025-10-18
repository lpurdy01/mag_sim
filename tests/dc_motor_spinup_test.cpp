#include "motorsim/circuit.hpp"
#include "motorsim/grid.hpp"
#include "motorsim/ingest.hpp"
#include "motorsim/mechanical.hpp"
#include "motorsim/probes.hpp"
#include "motorsim/solver.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <unordered_map>
#include <vector>

namespace {
constexpr double kPi = 3.14159265358979323846;
}

int main() {
    using namespace motorsim;
    namespace fs = std::filesystem;

    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/dc_motor_spinup_test.json").lexically_normal();

    ScenarioSpec spec;
    try {
        spec = loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load DC motor spin-up scenario: " << ex.what() << '\n';
        return 1;
    }

    std::vector<ScenarioFrame> frames = expandScenarioTimeline(spec);
    if (frames.empty()) {
        std::cerr << "Spin-up scenario produced no timeline frames\n";
        return 1;
    }

    CircuitSimulator circuit(frames);
    MechanicalSimulator mechanical;
    mechanical.initialize(spec, frames);

    if (!mechanical.is_active()) {
        std::cerr << "Mechanical simulator should be active for DC spin-up test\n";
        return 1;
    }

    const bool circuitsActive = circuit.is_active();

    SolveOptions options{};
    options.kind = SolverKind::CG;
    options.maxIters = 40000;
    options.tol = 6.5e-6;

    std::vector<double> warmAz;
    bool haveWarmStart = false;

    for (std::size_t idx = 0; idx < frames.size(); ++idx) {
        mechanical.apply_state(frames[idx]);
        circuit.update_for_frame(frames[idx], &mechanical);
        if (circuitsActive) {
            circuit.apply_currents(frames[idx]);
        }

        Grid2D grid(frames[idx].spec.nx, frames[idx].spec.ny, frames[idx].spec.dx, frames[idx].spec.dy);
        try {
            rasterizeScenarioToGrid(frames[idx].spec, grid);
        } catch (const std::exception& ex) {
            std::cerr << "Frame " << idx << ": rasterisation error: " << ex.what() << '\n';
            return 1;
        }

        InitialGuess guess{};
        const InitialGuess* guessPtr = nullptr;
        if (haveWarmStart && warmAz.size() == grid.Az.size()) {
            guess.Az0 = &warmAz;
            guessPtr = &guess;
        }

        const SolveResult result = solveAz(grid, options, guessPtr, nullptr);
        if (!result.converged) {
            std::cerr << "Frame " << idx << ": solver did not converge (relResidual=" << result.relResidual << ")\n";
            return 1;
        }

        computeB(grid);

        std::unordered_map<std::string, StressTensorResult> torqueSamples;
        for (const auto& probe : frames[idx].spec.outputs.probes) {
            try {
                const StressTensorResult sample = evaluate_maxwell_stress_probe(
                    grid, frames[idx].spec.originX, frames[idx].spec.originY, frames[idx].spec.dx, frames[idx].spec.dy,
                    probe.loopXs, probe.loopYs);
                torqueSamples.emplace(probe.id, sample);
            } catch (const std::exception& ex) {
                std::cerr << "Frame " << idx << ": torque probe evaluation failed: " << ex.what() << '\n';
                return 1;
            }
        }

        ScenarioFrame* nextFrame = (idx + 1 < frames.size()) ? &frames[idx + 1] : nullptr;
        mechanical.handle_solved_frame(frames[idx], torqueSamples, nextFrame);

        if (circuitsActive) {
            circuit.record_solved_frame(frames[idx], grid);
            if (nextFrame) {
                circuit.prepare_next_frame(idx, frames[idx], nextFrame);
            }
        }

        warmAz = grid.Az;
        haveWarmStart = true;
    }

    const auto& histories = mechanical.history();
    const auto it = histories.find("dc_rotor");
    if (it == histories.end() || it->second.size() < 2) {
        std::cerr << "Rotor history missing or too short for DC spin-up test\n";
        return 1;
    }

    const auto& samples = it->second;
    const double initialAngle = samples.front().angleRad;
    const double finalAngle = samples.back().angleRad;
    const double initialOmega = samples.front().omega;
    const double finalOmega = samples.back().omega;

    const double angleRiseDeg = (finalAngle - initialAngle) * 180.0 / kPi;
    const double speedRise = finalOmega - initialOmega;

    if (angleRiseDeg < 8.0) {
        std::cerr << "Rotor angle rise too small: " << angleRiseDeg << " deg\n";
        return 1;
    }
    if (speedRise < 5.0) {
        std::cerr << "Rotor speed rise too small: " << speedRise << " rad/s\n";
        return 1;
    }
    if (std::abs(finalOmega) <= 0.0) {
        std::cerr << "Final rotor speed magnitude too small: " << finalOmega << " rad/s\n";
        return 1;
    }

    const double finalRpm = std::abs(finalOmega) * 60.0 / (2.0 * kPi);
    std::cout << "DCMotorSpinupTest: angleRiseDeg=" << angleRiseDeg << " speedRise=" << speedRise
              << " final|omega|=" << std::abs(finalOmega) << " rad/s (" << finalRpm << " rpm)\n";

    return 0;
}
