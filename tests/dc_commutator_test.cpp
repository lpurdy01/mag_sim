#include "motorsim/circuit.hpp"
#include "motorsim/ingest.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {
constexpr double kPi = 3.14159265358979323846;

double expected_factor(double angle_deg) {
    const double cos_val = std::cos(angle_deg * kPi / 180.0);
    return (cos_val >= 0.0) ? 1.0 : -1.0;
}

}  // namespace

int main() {
    using namespace motorsim;
    namespace fs = std::filesystem;

    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/dc_commutator_test.json").lexically_normal();

    ScenarioSpec spec;
    try {
        spec = loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load DC commutator scenario: " << ex.what() << '\n';
        return 1;
    }

    std::vector<ScenarioFrame> frames = expandScenarioTimeline(spec);
    if (frames.empty()) {
        std::cerr << "Commutator scenario produced no frames\n";
        return 1;
    }

    CircuitSimulator circuit(frames);
    if (!circuit.is_active()) {
        std::cerr << "Circuit simulator should be active for commutator test\n";
        return 1;
    }

    std::unordered_map<std::string, std::size_t> regionIndex;
    for (std::size_t idx = 0; idx < spec.currentRegions.size(); ++idx) {
        regionIndex.emplace(spec.currentRegions[idx].id, idx);
    }
    auto findIndex = [&](const std::string& id) -> std::size_t {
        auto it = regionIndex.find(id);
        if (it == regionIndex.end()) {
            std::cerr << "Missing current region id '" << id << "' in scenario" << '\n';
            std::exit(1);
        }
        return it->second;
    };

    const std::size_t regionA = findIndex("armature_a");
    const std::size_t regionB = findIndex("armature_b");

    const double baseOrientationA = spec.currentRegions[regionA].orientation;
    const double baseOrientationB = spec.currentRegions[regionB].orientation;

    bool sawPositive = false;
    bool sawNegative = false;
    int signTransitions = 0;
    double previousFactor = 0.0;

    for (auto& frame : frames) {
        circuit.update_for_frame(frame, nullptr);
        circuit.apply_currents(frame);

        if (frame.spec.rotors.empty()) {
            std::cerr << "Frame missing rotor definition for commutator test\n";
            return 1;
        }

        const double angleDeg = frame.spec.rotors.front().hasCurrentAngle
                                     ? frame.spec.rotors.front().currentAngleDeg
                                     : frame.spec.rotors.front().initialAngleDeg;
        const double factor = expected_factor(angleDeg);
        if (factor > 0.0) {
            sawPositive = true;
        } else {
            sawNegative = true;
        }
        if (previousFactor != 0.0 && factor * previousFactor < 0.0) {
            ++signTransitions;
        }
        previousFactor = factor;

        const double orientationA = frame.spec.currentRegions[regionA].orientation;
        const double orientationB = frame.spec.currentRegions[regionB].orientation;

        if (std::abs(orientationA - baseOrientationA * factor) > 1e-9) {
            std::cerr << "Orientation mismatch for armature_a at angle " << angleDeg << " deg\n";
            return 1;
        }
        if (std::abs(orientationB - baseOrientationB * factor) > 1e-9) {
            std::cerr << "Orientation mismatch for armature_b at angle " << angleDeg << " deg\n";
            return 1;
        }
    }

    if (!sawPositive || !sawNegative) {
        std::cerr << "Commutator test did not cover both polarities" << '\n';
        return 1;
    }
    if (signTransitions == 0) {
        std::cerr << "Commutator test did not cross a commutation boundary" << '\n';
        return 1;
    }

    std::cout << "DCCommutatorTest: verified orientation across " << frames.size() << " frames\n";
    return 0;
}
