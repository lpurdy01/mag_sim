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
constexpr double kElectricalHz = 60.0;
}

int main() {
    using namespace motorsim;
    namespace fs = std::filesystem;

    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/induction_spinup_test.json").lexically_normal();

    ScenarioSpec spec;
    try {
        spec = loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load induction motor scenario: " << ex.what() << '\n';
        return 1;
    }

    std::vector<ScenarioFrame> frames = expandScenarioTimeline(spec);
    if (frames.empty()) {
        std::cerr << "Induction scenario produced no timeline frames\n";
        return 1;
    }

    if (!spec.transient.enabled) {
        std::cerr << "Induction scenario must enable transient stepping\n";
        return 1;
    }

    const double dt = spec.transient.dt > 0.0 ? spec.transient.dt : (frames.size() > 1 ? frames[1].time - frames[0].time : 0.0);
    if (!(dt > 0.0)) {
        std::cerr << "Transient timestep must be positive for induction spin-up test\n";
        return 1;
    }

    MechanicalSimulator mechanicalSim;
    mechanicalSim.initialize(spec, frames);
    if (!mechanicalSim.is_active()) {
        std::cerr << "Mechanical simulator should be active for induction spin-up test\n";
        return 1;
    }

    SolveOptions options{};
    options.kind = SolverKind::CG;
    options.maxIters = 35000;
    options.tol = 5e-6;
    options.preconditioner = PreconditionerKind::None;

    std::vector<double> transientState;
    bool haveTransientState = false;

    Grid2D grid(frames.front().spec.nx, frames.front().spec.ny, frames.front().spec.dx, frames.front().spec.dy);

    double torqueSamples = 0.0;
    double positiveTorqueSamples = 0.0;
    double negativeTorqueSamples = 0.0;

    for (std::size_t idx = 0; idx < frames.size(); ++idx) {
        ScenarioFrame& frame = frames[idx];
        mechanicalSim.apply_state(frame);

        if (grid.nx != frame.spec.nx || grid.ny != frame.spec.ny) {
            grid.resize(frame.spec.nx, frame.spec.ny, frame.spec.dx, frame.spec.dy);
        }

        try {
            rasterizeScenarioToGrid(frame.spec, grid);
        } catch (const std::exception& ex) {
            std::cerr << "Frame " << idx << ": rasterisation error: " << ex.what() << '\n';
            return 1;
        }

        if (!haveTransientState || transientState.size() != grid.Az.size()) {
            transientState.assign(grid.Az.size(), 0.0);
            haveTransientState = true;
        }

        const SolveResult report = solveTransientStep(grid, options, dt, transientState, nullptr);
        if (!report.converged) {
            std::cerr << "Frame " << idx << ": transient solve did not converge (relResidual=" << report.relResidual
                      << ")\n";
            return 1;
        }

        computeB(grid);

        std::unordered_map<std::string, StressTensorResult> torqueResults;
        for (const auto& probe : frame.spec.outputs.probes) {
            try {
                const StressTensorResult sample = evaluate_maxwell_stress_probe(
                    grid, frame.spec.originX, frame.spec.originY, frame.spec.dx, frame.spec.dy, probe.loopXs, probe.loopYs);
                torqueResults.emplace(probe.id, sample);
                ++torqueSamples;
                if (sample.torqueZ > 0.0) {
                    positiveTorqueSamples += 1.0;
                } else if (sample.torqueZ < 0.0) {
                    negativeTorqueSamples += 1.0;
                }
            } catch (const std::exception& ex) {
                std::cerr << "Frame " << idx << ": torque probe evaluation failed: " << ex.what() << '\n';
                return 1;
            }
        }

        ScenarioFrame* nextFrame = (idx + 1 < frames.size()) ? &frames[idx + 1] : nullptr;
        mechanicalSim.handle_solved_frame(frame, torqueResults, nextFrame);

        transientState = grid.Az;
    }

    const auto& histories = mechanicalSim.history();
    const auto histIt = histories.find("induction_rotor");
    if (histIt == histories.end() || histIt->second.size() < 2) {
        std::cerr << "Rotor history missing or too short for induction spin-up test\n";
        return 1;
    }

    const auto& samples = histIt->second;
    const double initialAngleRad = samples.front().angleRad;
    const double finalAngleRad = samples.back().angleRad;
    const double initialOmega = samples.front().omega;
    const double finalOmega = samples.back().omega;

    const double angleDeltaDeg = (finalAngleRad - initialAngleRad) * 180.0 / kPi;
    const double angleRiseDeg = std::abs(angleDeltaDeg);
    const double initialOmegaMag = std::abs(initialOmega);
    const double finalOmegaMag = std::abs(finalOmega);
    const double speedDelta = finalOmega - initialOmega;
    const double speedRise = finalOmegaMag - initialOmegaMag;

    if (angleRiseDeg < 5.0) {
        std::cerr << "Rotor angle change too small: |Δθ|=" << angleRiseDeg << " deg\n";
        return 1;
    }
    if (speedRise < 5.0) {
        std::cerr << "Rotor speed change too small: |Δω|=" << speedRise << " rad/s\n";
        return 1;
    }
    if (!(finalOmegaMag > 0.0)) {
        std::cerr << "Final rotor speed magnitude must be positive\n";
        return 1;
    }

    const double synchronousOmega = 2.0 * kPi * kElectricalHz;
    if (!(finalOmegaMag < synchronousOmega * 0.8)) {
        std::cerr << "Final rotor speed magnitude too close to synchronous speed: " << finalOmegaMag << " rad/s\n";
        return 1;
    }

    const double slip = 1.0 - (finalOmegaMag / synchronousOmega);
    if (!(slip > 0.05 && slip < 0.95)) {
        std::cerr << "Slip outside expected bounds: " << slip << '\n';
        return 1;
    }

    if (torqueSamples > 0.0) {
        const double dominant = std::max(positiveTorqueSamples, negativeTorqueSamples);
        const double dominantFraction = dominant / torqueSamples;
        if (dominantFraction < 0.5) {
            std::cerr << "Torque samples do not exhibit a dominant sign: fraction=" << dominantFraction << '\n';
            return 1;
        }
    } else {
        std::cerr << "Torque probe did not yield any samples\n";
        return 1;
    }

    const double finalRpm = finalOmegaMag * 60.0 / (2.0 * kPi);
    const char* directionLabel = (finalOmega >= 0.0) ? "positive" : "negative";
    std::cout << "InductionSpinupTest: angleDeltaDeg=" << angleDeltaDeg << " |Δθ|=" << angleRiseDeg
              << " speedDelta=" << speedDelta << " |Δω|=" << speedRise << " finalOmegaMag=" << finalOmegaMag
              << " rad/s (" << finalRpm << " rpm) direction=" << directionLabel << " slip=" << slip << '\n';

    return 0;
}
