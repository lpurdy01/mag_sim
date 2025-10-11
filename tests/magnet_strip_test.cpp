#include "motorsim/ingest.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace {
std::pair<double, double> B_from_line(double x, double y, double xs, double ys, double current, double mu) {
    const double dx = x - xs;
    const double dy = y - ys;
    const double r2 = dx * dx + dy * dy;
    if (r2 <= 0.0) {
        return {0.0, 0.0};
    }
    const double coeff = (mu * current) / (2.0 * M_PI * r2);
    return {-coeff * dy, coeff * dx};
}

std::pair<double, double> reference_field(double x, double y, double magnetMinX, double magnetMaxX,
                                          double magnetMinY, double magnetMaxY, double magnetisation) {
    constexpr int segments = 4000;
    const double dy = (magnetMaxY - magnetMinY) / static_cast<double>(segments);
    const double leftX = magnetMinX;
    const double rightX = magnetMaxX;
    double bx_total = 0.0;
    double by_total = 0.0;
    for (int k = 0; k < segments; ++k) {
        const double yk = magnetMinY + (k + 0.5) * dy;
        const auto left = B_from_line(x, y, leftX, yk, magnetisation * dy, motorsim::MU0);
        const auto right = B_from_line(x, y, rightX, yk, -magnetisation * dy, motorsim::MU0);
        bx_total += left.first + right.first;
        by_total += left.second + right.second;
    }
    return {bx_total, by_total};
}

}  // namespace

int main() {
    using namespace motorsim;

    namespace fs = std::filesystem;
    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/magnet_strip_test.json").lexically_normal();
    ScenarioSpec spec = loadScenarioFromJson(scenarioPath.string());

    if (spec.magnetRegions.size() != 1) {
        std::cerr << "magnet_strip_test: scenario must contain exactly one magnet region\n";
        return 1;
    }
    const ScenarioSpec::MagnetRegion magnet = spec.magnetRegions.front();

    SolveOptions options{};
    options.maxIters = 25000;
    options.tol = 5e-6;
    options.omega = 1.6;

    const double sampleY = 0.0;
    const std::vector<double> sampleX = {-0.08, -0.06, -0.04, 0.0, 0.04, 0.06, 0.08};
    const std::size_t j_center =
        static_cast<std::size_t>(std::llround((sampleY - spec.originY) / spec.dy));
    if (j_center >= spec.ny) {
        std::cerr << "magnet_strip_test: sample row out of range\n";
        return 1;
    }

    double baselineResidual = 0.0;
    double baselineRelErr = 0.0;
    double baselineHy = 0.0;
    bool baselineRecorded = false;

    const std::array<std::pair<SolverKind, const char*>, 2> solverCases = {
        std::make_pair(SolverKind::SOR, "SOR"), std::make_pair(SolverKind::CG, "CG")};

    for (const auto& [solverKind, label] : solverCases) {
        options.kind = solverKind;
        Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
        rasterizeScenarioToGrid(spec, grid);

        const SolveResult report = solveAz(grid, options);
        if (!report.converged || !(report.relResidual < options.tol)) {
            std::cerr << "magnet_strip_test: " << label << " solver did not converge to tolerance\n";
            return 1;
        }

        computeB(grid);
        computeH(grid);

        double diff2 = 0.0;
        double ref2 = 0.0;
        for (double x : sampleX) {
            const std::size_t i = static_cast<std::size_t>(std::llround((x - spec.originX) / spec.dx));
            if (i >= spec.nx) {
                std::cerr << "magnet_strip_test: sample column out of range\n";
                return 1;
            }
            const std::size_t idx = grid.idx(i, j_center);
            const auto ref =
                reference_field(x, sampleY, magnet.min_x, magnet.max_x, magnet.min_y, magnet.max_y, magnet.My);
            const double by_num = grid.By[idx];
            diff2 += (by_num - ref.second) * (by_num - ref.second);
            ref2 += ref.second * ref.second;
        }

        const double relErr = (ref2 > 0.0) ? std::sqrt(diff2 / ref2) : 0.0;
        if (!(relErr < 0.15)) {
            std::cerr << "magnet_strip_test (" << label << "): field RMS relErr=" << relErr
                      << " exceeds tolerance\n";
            return 1;
        }

        const std::size_t i_center = static_cast<std::size_t>(std::llround((-spec.originX) / spec.dx));
        const std::size_t center_idx = grid.idx(i_center, j_center);
        const double Hy_center = grid.Hy[center_idx];
        if (!(std::abs(Hy_center) < 3.0e5)) {
            std::cerr << "magnet_strip_test (" << label
                      << "): Hy inside magnet expected near zero, got " << Hy_center << '\n';
            return 1;
        }

        std::cout << "magnet_strip_test (" << label << "): relErr=" << relErr
                  << ", Hy_center=" << Hy_center << ", iters=" << report.iters << '\n';

        if (!baselineRecorded) {
            baselineResidual = report.relResidual;
            baselineRelErr = relErr;
            baselineHy = std::abs(Hy_center);
            baselineRecorded = true;
        } else {
            const double residualBound = baselineResidual * 1.1 + 1e-12;
            if (report.relResidual > residualBound) {
                std::cerr << "magnet_strip_test: CG residual " << report.relResidual
                          << " exceeds parity bound " << residualBound << '\n';
                return 1;
            }
            if (relErr > baselineRelErr + 0.02) {
                std::cerr << "magnet_strip_test: CG relErr " << relErr
                          << " exceeds parity slack relative to SOR " << baselineRelErr << '\n';
                return 1;
            }
            if (std::abs(Hy_center) > baselineHy + 5.0e4) {
                std::cerr << "magnet_strip_test: CG |Hy| " << std::abs(Hy_center)
                          << " exceeds parity slack relative to SOR " << baselineHy << '\n';
                return 1;
            }
        }
    }

    return 0;
}
