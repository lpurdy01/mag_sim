#include "motorsim/ingest.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>

namespace {

bool approxEqual(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) <= tol;
}

}

int main() {
    using namespace motorsim;
    namespace fs = std::filesystem;

    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/timeline_frames_test.json")
            .lexically_normal();

    ScenarioSpec spec;
    try {
        spec = loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load scenario: " << ex.what() << '\n';
        return 1;
    }

    if (spec.wires.size() != 2) {
        std::cerr << "Expected exactly two wires in base scenario\n";
        return 1;
    }
    if (spec.magnetRegions.size() != 1) {
        std::cerr << "Expected a single magnet region in base scenario\n";
        return 1;
    }
    if (!approxEqual(spec.wires[0].current, 10.0) || !approxEqual(spec.wires[1].current, -10.0)) {
        std::cerr << "Base wire currents were not parsed correctly\n";
        return 1;
    }
    if (!approxEqual(spec.magnetRegions[0].Mx, 1.0) ||
        !approxEqual(spec.magnetRegions[0].My, 0.0)) {
        std::cerr << "Base magnetisation vector unexpected\n";
        return 1;
    }

    std::vector<ScenarioFrame> frames;
    try {
        frames = expandScenarioTimeline(spec);
    } catch (const std::exception& ex) {
        std::cerr << "Timeline expansion failed: " << ex.what() << '\n';
        return 1;
    }

    if (frames.size() != 3) {
        std::cerr << "Expected three frames after expansion, got " << frames.size() << '\n';
        return 1;
    }

    for (const auto& frame : frames) {
        if (!frame.spec.timeline.empty()) {
            std::cerr << "Frame spec should not retain nested timeline data\n";
            return 1;
        }
        if (frame.spec.wires.size() != spec.wires.size()) {
            std::cerr << "Frame wire count mismatch\n";
            return 1;
        }
        if (frame.spec.magnetRegions.size() != spec.magnetRegions.size()) {
            std::cerr << "Frame magnet region count mismatch\n";
            return 1;
        }
    }

    const ScenarioFrame& frame0 = frames[0];
    if (!approxEqual(frame0.time, 0.0) || !approxEqual(frame0.spec.wires[0].current, 10.0) ||
        !approxEqual(frame0.spec.wires[1].current, -10.0)) {
        std::cerr << "Frame 0 currents/time incorrect\n";
        return 1;
    }
    if (!approxEqual(frame0.spec.magnetRegions[0].Mx, 1.0) ||
        !approxEqual(frame0.spec.magnetRegions[0].My, 0.0)) {
        std::cerr << "Frame 0 magnetisation mismatch\n";
        return 1;
    }

    const ScenarioFrame& frame1 = frames[1];
    if (!approxEqual(frame1.time, 1.0e-3) || !approxEqual(frame1.spec.wires[0].current, 5.0) ||
        !approxEqual(frame1.spec.wires[1].current, -5.0)) {
        std::cerr << "Frame 1 overrides not applied\n";
        return 1;
    }
    if (!approxEqual(frame1.spec.magnetRegions[0].Mx, 0.0, 1e-9) ||
        !approxEqual(frame1.spec.magnetRegions[0].My, 1.0, 1e-9)) {
        std::cerr << "Frame 1 rotor rotation incorrect\n";
        return 1;
    }

    const ScenarioFrame& frame2 = frames[2];
    if (!approxEqual(frame2.time, 2.0e-3) || !approxEqual(frame2.spec.wires[0].current, 0.0) ||
        !approxEqual(frame2.spec.wires[1].current, 2.5)) {
        std::cerr << "Frame 2 wire overrides incorrect\n";
        return 1;
    }

    const double expectedMx = std::cos(45.0 * 3.14159265358979323846 / 180.0);
    const double expectedMy = std::sin(45.0 * 3.14159265358979323846 / 180.0);
    if (!approxEqual(frame2.spec.magnetRegions[0].Mx, expectedMx, 1e-9) ||
        !approxEqual(frame2.spec.magnetRegions[0].My, expectedMy, 1e-9)) {
        std::cerr << "Frame 2 magnet angle override incorrect\n";
        return 1;
    }

    if (!approxEqual(spec.wires[0].current, 10.0) || !approxEqual(spec.wires[1].current, -10.0) ||
        !approxEqual(spec.magnetRegions[0].Mx, 1.0) || !approxEqual(spec.magnetRegions[0].My, 0.0)) {
        std::cerr << "Base scenario mutated by timeline expansion\n";
        return 1;
    }

    return 0;
}
