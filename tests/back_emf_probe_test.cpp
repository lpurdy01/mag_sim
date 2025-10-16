#include "motorsim/ingest.hpp"
#include "motorsim/probes.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

int main() {
    using namespace motorsim;
    namespace fs = std::filesystem;

    constexpr double kPi = 3.14159265358979323846;

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
    const double rectExpectedFlux = 2.0 * spec.dx * spec.dy;
    if (!approxEqual(rectFlux, rectExpectedFlux)) {
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

    const double phiAmplitude = 2.5e-3;
    const double electricalHz = 60.0;
    const double omega = 2.0 * kPi * electricalHz;
    const std::size_t samplesPerCycle = 48;
    const std::size_t frameCount = samplesPerCycle + 1;
    const double dt = 1.0 / electricalHz / samplesPerCycle;

    std::vector<std::size_t> seriesFrames(frameCount);
    std::vector<double> seriesTimes(frameCount);
    std::vector<double> fluxPhaseA(frameCount);
    std::vector<double> fluxPhaseB(frameCount);
    std::vector<double> fluxPhaseC(frameCount);
    for (std::size_t n = 0; n < frameCount; ++n) {
        const double t = static_cast<double>(n) * dt;
        seriesFrames[n] = n;
        seriesTimes[n] = t;
        fluxPhaseA[n] = phiAmplitude * std::sin(omega * t);
        fluxPhaseB[n] = phiAmplitude * std::sin(omega * t - 2.0 * kPi / 3.0);
        fluxPhaseC[n] = phiAmplitude * std::sin(omega * t + 2.0 * kPi / 3.0);
    }

    const double expectedAmplitude = phiAmplitude * omega;

    const auto checkSinusoid = [&](const std::vector<double>& flux, double phaseShift) {
        std::vector<BackEmfSample> emfSeries = compute_back_emf_series(seriesFrames, seriesTimes, flux);
        if (emfSeries.size() != frameCount - 1) {
            std::cerr << "Expected " << (frameCount - 1) << " EMF samples, got " << emfSeries.size() << '\n';
            return false;
        }
        double maxAbsEmf = 0.0;
        const double tolerance = 1.5e-3 * std::max(expectedAmplitude, 1.0);
        for (std::size_t idx = 0; idx < emfSeries.size(); ++idx) {
            const auto& interval = emfSeries[idx];
            const double tMid = 0.5 * (interval.timeStart + interval.timeEnd);
            const double expected = -phiAmplitude * omega * std::cos(omega * tMid + phaseShift);
            if (std::abs(interval.emf - expected) > tolerance) {
                std::cerr << "Sinusoidal EMF mismatch at index " << idx << ": got " << interval.emf
                          << " expected " << expected << '\n';
                return false;
            }
            maxAbsEmf = std::max(maxAbsEmf, std::abs(interval.emf));
        }
        if (std::abs(maxAbsEmf - expectedAmplitude) / expectedAmplitude > 0.01) {
            std::cerr << "EMF amplitude mismatch: got " << maxAbsEmf << " expected " << expectedAmplitude << '\n';
            return false;
        }
        return true;
    };

    if (!checkSinusoid(fluxPhaseA, 0.0)) {
        return 1;
    }
    if (!checkSinusoid(fluxPhaseB, -2.0 * kPi / 3.0)) {
        return 1;
    }
    if (!checkSinusoid(fluxPhaseC, 2.0 * kPi / 3.0)) {
        return 1;
    }

    const auto relError = [](double value, double reference) {
        const double denom = std::max(std::abs(reference), 1e-12);
        return std::abs(value - reference) / denom;
    };

    const double polygonFluxRelErr = relError(polygonFlux, expectedFlux);
    const double rectFluxRelErr = relError(rectFlux, rectExpectedFlux);
    const double emfRelErr = relError(sample.emf, expectedEmf);

    std::cout << "BackEmfTest: polygon_flux=" << polygonFlux << " rect_flux=" << rectFlux
              << " emf=" << sample.emf << " polygon_flux_relErr=" << polygonFluxRelErr
              << " rect_flux_relErr=" << rectFluxRelErr << " emf_relErr=" << emfRelErr
              << " sinusoid_amp=" << expectedAmplitude << '\n';

    return 0;
}
