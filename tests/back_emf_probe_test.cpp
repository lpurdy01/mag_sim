#include "motorsim/ingest.hpp"
#include "motorsim/probes.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

int main() {
    using namespace motorsim;
    namespace fs = std::filesystem;

    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/back_emf_probe_test.json").lexically_normal();

    ScenarioSpec spec;
    try {
        spec = loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load back-EMF scenario: " << ex.what() << '\n';
        return 1;
    }

    if (spec.outputs.backEmfProbes.size() != 2) {
        std::cerr << "Expected two back_emf_probe outputs\n";
        return 1;
    }

    const auto& polygonProbe = spec.outputs.backEmfProbes[0];
    if (polygonProbe.component != FluxComponent::Bx) {
        std::cerr << "Polygon probe component mismatch\n";
        return 1;
    }
    if (polygonProbe.shape != ScenarioSpec::Outputs::BackEmfProbe::RegionShape::Polygon) {
        std::cerr << "Polygon probe shape mismatch\n";
        return 1;
    }
    if (polygonProbe.frameIndices.size() != 2 || polygonProbe.frameIndices[0] != 0 ||
        polygonProbe.frameIndices[1] != 1) {
        std::cerr << "Polygon probe frames not parsed correctly\n";
        return 1;
    }

    const auto& rectProbe = spec.outputs.backEmfProbes[1];
    if (rectProbe.component != FluxComponent::Bmag) {
        std::cerr << "Rect probe default component expected to be Bmag\n";
        return 1;
    }
    if (!rectProbe.frameIndices.empty()) {
        std::cerr << "Rect probe should default to all frames when frames array omitted\n";
        return 1;
    }
    if (rectProbe.shape != ScenarioSpec::Outputs::BackEmfProbe::RegionShape::Rect) {
        std::cerr << "Rect probe shape mismatch\n";
        return 1;
    }

    Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    const std::size_t count = grid.nx * grid.ny;
    grid.Bx.assign(count, 1.0);
    grid.By.assign(count, 0.0);

    double polygonFlux = 0.0;
    double rectFlux = 0.0;
    try {
        polygonFlux = integrate_polygon_flux_component(grid, spec.originX, spec.originY, spec.dx, spec.dy,
                                                        polygonProbe.xs, polygonProbe.ys, polygonProbe.minX,
                                                        polygonProbe.maxX, polygonProbe.minY, polygonProbe.maxY,
                                                        polygonProbe.component);
        rectFlux = integrate_rect_flux_component(grid, spec.originX, spec.originY, spec.dx, spec.dy,
                                                  rectProbe.minX, rectProbe.maxX, rectProbe.minY, rectProbe.maxY,
                                                  rectProbe.component);
    } catch (const std::exception& ex) {
        std::cerr << "Flux integration failed: " << ex.what() << '\n';
        return 1;
    }

    const double expectedFlux = 4.0 * spec.dx * spec.dy;  // four cells selected by both regions
    const auto approxEqual = [](double a, double b) {
        return std::abs(a - b) <= 1e-12;
    };
    if (!approxEqual(polygonFlux, expectedFlux)) {
        std::cerr << "Unexpected polygon flux: " << polygonFlux << " expected " << expectedFlux << '\n';
        return 1;
    }
    if (!approxEqual(rectFlux, 2.0 * spec.dx * spec.dy)) {
        // 2 cells * area * Bmag (Bmag = 1.0 in this synthetic field)
        std::cerr << "Unexpected rect flux: " << rectFlux << '\n';
        return 1;
    }

    std::vector<std::size_t> frames{0, 1};
    std::vector<double> times{spec.timeline[0].time, spec.timeline[1].time};
    std::vector<double> fluxes{0.0, polygonFlux};
    std::vector<BackEmfSample> samples;
    try {
        samples = compute_back_emf_series(frames, times, fluxes);
    } catch (const std::exception& ex) {
        std::cerr << "compute_back_emf_series failed: " << ex.what() << '\n';
        return 1;
    }
    if (samples.size() != 1) {
        std::cerr << "Expected a single back-EMF sample\n";
        return 1;
    }
    const BackEmfSample& sample = samples.front();
    const double expectedEmf = -(polygonFlux - 0.0) / (times[1] - times[0]);
    if (!approxEqual(sample.emf, expectedEmf)) {
        std::cerr << "Unexpected EMF value: " << sample.emf << " expected " << expectedEmf << '\n';
        return 1;
    }

    return 0;
}
