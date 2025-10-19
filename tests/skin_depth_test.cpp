#include "motorsim/ingest.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <cmath>
#include <filesystem>
#include <iostream>
#include <vector>

namespace {
constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
}  // namespace

int main() {
    namespace fs = std::filesystem;

    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/skin_depth_test.json").lexically_normal();

    motorsim::ScenarioSpec spec;
    try {
        spec = motorsim::loadScenarioFromJson(scenarioPath.string());
    } catch (const std::exception& ex) {
        std::cerr << "Failed to load skin-depth scenario: " << ex.what() << '\n';
        return 1;
    }

    std::vector<motorsim::ScenarioFrame> frames;
    try {
        frames = motorsim::expandScenarioTimeline(spec);
    } catch (const std::exception& ex) {
        std::cerr << "Failed to expand scenario timeline: " << ex.what() << '\n';
        return 1;
    }
    if (frames.size() != 1) {
        std::cerr << "Expected exactly one frame in skin-depth scenario, got " << frames.size() << '\n';
        return 1;
    }

    motorsim::ScenarioFrame& frame = frames.front();
    motorsim::Grid2D grid(frame.spec.nx, frame.spec.ny, frame.spec.dx, frame.spec.dy);
    try {
        motorsim::rasterizeScenarioToGrid(frame.spec, grid);
    } catch (const std::exception& ex) {
        std::cerr << "Rasterisation failure: " << ex.what() << '\n';
        return 1;
    }

    const auto& solverSettings = frame.spec.solverSettings;
    double omega = 0.0;
    if (solverSettings.harmonicOmegaSpecified) {
        omega = solverSettings.harmonicOmega;
    } else if (solverSettings.harmonicFrequencySpecified) {
        omega = solverSettings.harmonicFrequencyHz * kTwoPi;
    }
    if (!(omega > 0.0)) {
        std::cerr << "Harmonic frequency missing from scenario" << '\n';
        return 1;
    }

    motorsim::SolveOptions options{};
    options.kind = motorsim::SolverKind::Harmonic;
    options.maxIters = 4000;
    options.tol = 1e-6;
    options.progressEverySec = 0.0;
    options.snapshotEveryIters = 0;
    options.harmonicOmega = omega;

    const motorsim::SolveResult result = motorsim::solveAz(grid, options, nullptr, nullptr);
    if (!result.converged) {
        std::cerr << "Harmonic solve did not converge: iters=" << result.iters
                  << " relResidual=" << result.relResidual << '\n';
        return 1;
    }

    std::vector<double> BxReal;
    std::vector<double> BxImag;
    std::vector<double> ByReal;
    std::vector<double> ByImag;
    motorsim::computeBHarmonic(grid, grid.Az, grid.AzImag, BxReal, BxImag, ByReal, ByImag);

    const double interfaceX = 0.0;
    const double yMin = -0.02;
    const double yMax = 0.02;
    std::vector<double> distances;
    std::vector<double> logMagnitudes;

    for (std::size_t i = 0; i < grid.nx; ++i) {
        const double x = frame.spec.originX + static_cast<double>(i) * frame.spec.dx;
        if (x <= interfaceX + grid.dx) {
            continue;
        }
        double accumMag = 0.0;
        std::size_t count = 0;
        for (std::size_t j = 0; j < grid.ny; ++j) {
            const double y = frame.spec.originY + static_cast<double>(j) * frame.spec.dy;
            if (y < yMin || y > yMax) {
                continue;
            }
            const std::size_t idx = grid.idx(i, j);
            const double bxR = BxReal[idx];
            const double bxI = BxImag[idx];
            const double byR = ByReal[idx];
            const double byI = ByImag[idx];
            const double mag = std::sqrt(bxR * bxR + bxI * bxI + byR * byR + byI * byI);
            accumMag += mag;
            ++count;
        }
        if (count == 0) {
            continue;
        }
        const double avgMag = accumMag / static_cast<double>(count);
        if (!(avgMag > 0.0)) {
            continue;
        }
        distances.push_back(x - interfaceX);
        logMagnitudes.push_back(std::log(avgMag));
    }

    if (distances.size() < 6) {
        std::cerr << "Insufficient samples inside conductor: " << distances.size() << '\n';
        return 1;
    }

    double sumX = 0.0;
    double sumY = 0.0;
    double sumXX = 0.0;
    double sumXY = 0.0;
    for (std::size_t idx = 0; idx < distances.size(); ++idx) {
        const double x = distances[idx];
        const double y = logMagnitudes[idx];
        sumX += x;
        sumY += y;
        sumXX += x * x;
        sumXY += x * y;
    }
    const double n = static_cast<double>(distances.size());
    const double denom = n * sumXX - sumX * sumX;
    if (!(std::abs(denom) > 1e-12)) {
        std::cerr << "Regression fit ill-conditioned" << '\n';
        return 1;
    }
    const double slope = (n * sumXY - sumX * sumY) / denom;

    constexpr double sigma = 5.96e7;
    const double mu = motorsim::MU0;
    const double skinDepth = std::sqrt(2.0 / (omega * mu * sigma));
    const double expectedSlope = -1.0 / skinDepth;
    const double relativeError = std::abs(slope - expectedSlope) / std::abs(expectedSlope);

    if (!(relativeError < 0.15)) {
        std::cerr << "Skin-depth slope error too high: " << relativeError << '\n';
        return 1;
    }
    if (!(slope < 0.0)) {
        std::cerr << "Recovered slope should be negative but was " << slope << '\n';
        return 1;
    }

    return 0;
}
