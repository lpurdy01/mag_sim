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
        (fs::path(__FILE__).parent_path() / "../inputs/tests/pm_motor_spinup_test.json").lexically_normal();

    ScenarioSpec spec;
    try {
        spec = loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load PM motor spin-up scenario: " << ex.what() << '\n';
        return 1;
    }

    std::vector<ScenarioFrame> frames = expandScenarioTimeline(spec);
    if (frames.empty()) {
        std::cerr << "Spin-up scenario produced no timeline frames\n";
        return 1;
    }

    CircuitSimulator circuitSim(frames);
    MechanicalSimulator mechanicalSim;
    mechanicalSim.initialize(spec, frames);

    if (!mechanicalSim.is_active()) {
        std::cerr << "Mechanical simulator should be active for spin-up test\n";
        return 1;
    }

    const bool circuitsActive = circuitSim.is_active();

    SolveOptions options{};
    options.kind = SolverKind::CG;
    options.maxIters = 20000;
    options.tol = 2e-6;

    std::vector<double> warmAz;
    bool haveWarmStart = false;

    for (std::size_t idx = 0; idx < frames.size(); ++idx) {
        if (mechanicalSim.is_active()) {
            mechanicalSim.apply_state(frames[idx]);
        }
        if (circuitsActive) {
            circuitSim.update_for_frame(frames[idx], mechanicalSim.is_active() ? &mechanicalSim : nullptr);
        }
        if (circuitsActive) {
            circuitSim.apply_currents(frames[idx]);
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

        const SolveResult report = solveAz(grid, options, guessPtr, nullptr);
        if (!report.converged) {
            std::cerr << "Frame " << idx << ": solver did not converge (relResidual=" << report.relResidual << ")\n";
            return 1;
        }

        motorsim::computeB(grid);

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

        if (mechanicalSim.is_active()) {
            ScenarioFrame* nextFrame = (idx + 1 < frames.size()) ? &frames[idx + 1] : nullptr;
            mechanicalSim.handle_solved_frame(frames[idx], torqueSamples, nextFrame);
        }

        if (circuitsActive) {
            circuitSim.record_solved_frame(frames[idx], grid);
            if (idx + 1 < frames.size()) {
                circuitSim.prepare_next_frame(idx, frames[idx], &frames[idx + 1]);
            }
        }

        warmAz = grid.Az;
        haveWarmStart = true;
    }

    const auto& histories = mechanicalSim.history();
    const auto histIt = histories.find("pm_rotor");
    if (histIt == histories.end() || histIt->second.size() < 2) {
        std::cerr << "Rotor history missing or too short for spin-up test\n";
        return 1;
    }

    const auto& samples = histIt->second;
    const double initialAngleRad = samples.front().angleRad;
    const double finalAngleRad = samples.back().angleRad;
    const double initialOmega = samples.front().omega;
    const double finalOmega = samples.back().omega;

    const double angleRiseDeg = (finalAngleRad - initialAngleRad) * 180.0 / kPi;
    const double speedRise = finalOmega - initialOmega;

    if (angleRiseDeg < 15.0) {
        std::cerr << "Rotor angle rise too small: " << angleRiseDeg << " deg\n";
        return 1;
    }
    if (speedRise < 10.0) {
        std::cerr << "Rotor speed rise too small: " << speedRise << " rad/s\n";
        return 1;
    }
    if (finalOmega <= 0.0) {
        std::cerr << "Final rotor speed is non-positive: " << finalOmega << " rad/s\n";
        return 1;
    }

    const double finalRpm = finalOmega * 60.0 / (2.0 * kPi);
    std::cout << "PMMotorSpinupTest: angleRiseDeg=" << angleRiseDeg << " speedRise=" << speedRise
              << " finalOmega=" << finalOmega << " rad/s (" << finalRpm << " rpm)\n";

    return 0;
}
