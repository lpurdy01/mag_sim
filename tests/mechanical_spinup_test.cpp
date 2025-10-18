#include "motorsim/ingest.hpp"
#include "motorsim/mechanical.hpp"
#include "motorsim/probes.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <unordered_map>

namespace {
constexpr double kPi = 3.14159265358979323846;
}

int main() {
    using namespace motorsim;
    namespace fs = std::filesystem;

    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/mechanical_spinup_test.json").lexically_normal();

    ScenarioSpec spec;
    try {
        spec = loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load mechanical spin-up scenario: " << ex.what() << '\n';
        return 1;
    }

    std::vector<ScenarioFrame> frames = expandScenarioTimeline(spec);
    if (frames.size() != 3) {
        std::cerr << "Expected three timeline frames\n";
        return 1;
    }

    MechanicalSimulator sim;
    sim.initialize(spec, frames);
    if (!sim.is_active()) {
        std::cerr << "Mechanical simulator should be active for test scenario\n";
        return 1;
    }

    sim.apply_state(frames[0]);

    StressTensorResult torqueSample{};
    torqueSample.torqueZ = 1.0;  // N·m constant torque
    std::unordered_map<std::string, StressTensorResult> torqueMap;
    torqueMap.emplace("test_torque", torqueSample);

    sim.handle_solved_frame(frames[0], torqueMap, &frames[1]);

    auto angleDeg = sim.rotor_angle_deg("test_rotor");
    auto omegaRad = sim.rotor_speed_rad_s("test_rotor");
    if (!angleDeg || !omegaRad) {
        std::cerr << "Rotor diagnostics unavailable after first step\n";
        return 1;
    }

    const double expectedOmega1 = 0.1;  // rad/s
    const double expectedAngle1 = 0.5 * (1.0 / 0.01) * 1e-6;  // 0.5 * a * dt^2
    if (std::abs(*omegaRad - expectedOmega1) > 1e-6) {
        std::cerr << "Unexpected omega after first step: " << *omegaRad << " expected " << expectedOmega1 << '\n';
        return 1;
    }
    const double angleRad1 = *angleDeg * kPi / 180.0;
    if (std::abs(angleRad1 - expectedAngle1) > 1e-7) {
        std::cerr << "Unexpected angle after first step: " << angleRad1 << " expected " << expectedAngle1 << '\n';
        return 1;
    }

    // Advance one more frame with the same torque
    sim.apply_state(frames[1]);
    sim.handle_solved_frame(frames[1], torqueMap, &frames[2]);

    angleDeg = sim.rotor_angle_deg("test_rotor");
    omegaRad = sim.rotor_speed_rad_s("test_rotor");
    if (!angleDeg || !omegaRad) {
        std::cerr << "Rotor diagnostics unavailable after second step\n";
        return 1;
    }

    const double expectedOmega2 = 0.2;  // rad/s
    if (std::abs(*omegaRad - expectedOmega2) > 1e-6) {
        std::cerr << "Unexpected omega after second step: " << *omegaRad << " expected " << expectedOmega2 << '\n';
        return 1;
    }
    const double angleRad2 = *angleDeg * kPi / 180.0;
    const double expectedAngle2 = expectedAngle1 + (expectedOmega1 * 0.001 + 0.5 * (1.0 / 0.01) * 1e-6);
    if (std::abs(angleRad2 - expectedAngle2) > 1e-7) {
        std::cerr << "Unexpected angle after second step: " << angleRad2 << " expected " << expectedAngle2 << '\n';
        return 1;
    }

    std::cout << "MechanicalSpinupTest: omega1=" << expectedOmega1 << " omega2=" << expectedOmega2
              << " angle1_deg=" << (*angleDeg) << '\n';

    return 0;
}
