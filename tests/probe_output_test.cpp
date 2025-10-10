#include "motorsim/ingest.hpp"
#include "motorsim/probes.hpp"
#include "motorsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

int main() {
    using namespace motorsim;
    namespace fs = std::filesystem;

    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/probe_output_test.json").lexically_normal();

    ScenarioSpec spec;
    try {
        spec = loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load probe scenario: " << ex.what() << '\n';
        return 1;
    }

    if (spec.outputs.probes.size() != 1) {
        std::cerr << "Expected a single probe output in scenario" << '\n';
        return 1;
    }

    const auto& probe = spec.outputs.probes.front();
    if (probe.quantity != ScenarioSpec::Outputs::Probe::Quantity::ForceAndTorque) {
        std::cerr << "Probe quantity was not parsed as force_and_torque" << '\n';
        return 1;
    }
    if (probe.method != ScenarioSpec::Outputs::Probe::Method::StressTensor) {
        std::cerr << "Probe method was not parsed as stress_tensor" << '\n';
        return 1;
    }
    if (probe.loopXs.size() != probe.loopYs.size() || probe.loopXs.size() < 3) {
        std::cerr << "Probe loop coordinates not populated correctly" << '\n';
        return 1;
    }

    Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    const std::size_t count = grid.nx * grid.ny;
    grid.Bx.assign(count, 0.0);
    grid.By.assign(count, 0.0);

    const double B0 = 0.2;
    for (std::size_t j = 0; j < grid.ny; ++j) {
        const double y = spec.originY + static_cast<double>(j) * spec.dy;
        const double value = (y >= 0.0) ? B0 : 0.0;
        for (std::size_t i = 0; i < grid.nx; ++i) {
            const std::size_t idx = grid.idx(i, j);
            grid.By[idx] = value;
        }
    }

    StressTensorResult result;
    try {
        result = evaluate_maxwell_stress_probe(grid, spec.originX, spec.originY, spec.dx, spec.dy,
                                               probe.loopXs, probe.loopYs);
    } catch (const std::exception& ex) {
        std::cerr << "Probe evaluation failed: " << ex.what() << '\n';
        return 1;
    }

    const double bottomLength = std::hypot(probe.loopXs[1] - probe.loopXs[0], probe.loopYs[1] - probe.loopYs[0]);
    const double expectedFy = 0.5 * (B0 * B0 / MU0) * bottomLength;

    const auto approxEqual = [](double a, double b, double tol) {
        return std::abs(a - b) <= tol;
    };

    const auto relError = [](double value, double reference) {
        const double denom = std::max(std::abs(reference), 1e-12);
        return std::abs(value - reference) / denom;
    };

    if (!approxEqual(result.forceX, 0.0, 1e-6)) {
        std::cerr << "Expected Fx ~ 0, got " << result.forceX << '\n';
        return 1;
    }
    if (!approxEqual(result.torqueZ, 0.0, 1e-6)) {
        std::cerr << "Expected Tz ~ 0, got " << result.torqueZ << '\n';
        return 1;
    }
    if (!approxEqual(result.forceY, expectedFy, 5.0)) {
        std::cerr << "Unexpected Fy value: got " << result.forceY << ", expected " << expectedFy << '\n';
        return 1;
    }

    std::vector<double> reversedXs(probe.loopXs.rbegin(), probe.loopXs.rend());
    std::vector<double> reversedYs(probe.loopYs.rbegin(), probe.loopYs.rend());
    StressTensorResult reversed = evaluate_maxwell_stress_probe(grid, spec.originX, spec.originY, spec.dx, spec.dy,
                                                                reversedXs, reversedYs);
    if (!approxEqual(reversed.forceX, result.forceX, 1e-6) ||
        !approxEqual(reversed.forceY, result.forceY, 5.0) ||
        !approxEqual(reversed.torqueZ, result.torqueZ, 1e-6)) {
        std::cerr << "Reversed loop produced inconsistent results" << '\n';
        return 1;
    }

    const double fyRelErr = relError(result.forceY, expectedFy);
    const double fxAbsErr = std::abs(result.forceX);
    const double tzAbsErr = std::abs(result.torqueZ);

    std::cout << "ProbeOutput: Fx=" << result.forceX << " Fy=" << result.forceY
              << " Tz=" << result.torqueZ << " Fy_relErr=" << fyRelErr
              << " Fx_absErr=" << fxAbsErr << " Tz_absErr=" << tzAbsErr << '\n';

    return 0;
}
