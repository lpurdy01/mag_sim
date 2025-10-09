#include "motorsim/ingest.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <cmath>
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

    const double Lx = 0.30;
    const double Ly = 0.20;
    const std::size_t nx = 241;
    const std::size_t ny = 201;
    const double dx = Lx / static_cast<double>(nx - 1);
    const double dy = Ly / static_cast<double>(ny - 1);

    ScenarioSpec spec{};
    spec.version = "0.2";
    spec.Lx = Lx;
    spec.Ly = Ly;
    spec.nx = nx;
    spec.ny = ny;
    spec.dx = dx;
    spec.dy = dy;
    spec.originX = -0.5 * Lx;
    spec.originY = -0.5 * Ly;
    spec.mu_r_background = 1.0;

    ScenarioSpec::MagnetRegion magnet{};
    magnet.shape = ScenarioSpec::MagnetRegion::Shape::Polygon;
    magnet.Mx = 0.0;
    magnet.My = 6.0e5;
    magnet.xs = {-0.02, -0.02, 0.02, 0.02};
    magnet.ys = {-0.04, 0.04, 0.04, -0.04};
    magnet.min_x = -0.02;
    magnet.max_x = 0.02;
    magnet.min_y = -0.04;
    magnet.max_y = 0.04;
    spec.magnetRegions.push_back(magnet);

    Grid2D grid(nx, ny, dx, dy);
    rasterizeScenarioToGrid(spec, grid);

    SolveOptions options{};
    options.maxIters = 25000;
    options.tol = 5e-6;
    options.omega = 1.6;

    const SolveReport report = solveAz_GS_SOR(grid, options);
    if (!report.converged || !(report.relResidual < options.tol)) {
        std::cerr << "magnet_strip_test: solver did not converge to tolerance\n";
        return 1;
    }

    computeB(grid);
    computeH(grid);

    const double sampleY = 0.0;
    const std::vector<double> sampleX = {-0.08, -0.06, -0.04, 0.0, 0.04, 0.06, 0.08};
    const std::size_t j_center = static_cast<std::size_t>(std::llround((sampleY - spec.originY) / spec.dy));
    if (j_center >= ny) {
        std::cerr << "magnet_strip_test: sample row out of range\n";
        return 1;
    }

    double diff2 = 0.0;
    double ref2 = 0.0;
    for (double x : sampleX) {
        const std::size_t i = static_cast<std::size_t>(std::llround((x - spec.originX) / spec.dx));
        if (i >= nx) {
            std::cerr << "magnet_strip_test: sample column out of range\n";
            return 1;
        }
        const std::size_t idx = grid.idx(i, j_center);
        const auto ref = reference_field(x, sampleY, magnet.min_x, magnet.max_x, magnet.min_y, magnet.max_y, magnet.My);
        const double by_num = grid.By[idx];
        diff2 += (by_num - ref.second) * (by_num - ref.second);
        ref2 += ref.second * ref.second;
    }

    const double relErr = (ref2 > 0.0) ? std::sqrt(diff2 / ref2) : 0.0;
    if (!(relErr < 0.15)) {
        std::cerr << "magnet_strip_test: field RMS relErr=" << relErr << " exceeds tolerance\n";
        return 1;
    }

    const std::size_t i_center = static_cast<std::size_t>(std::llround((-spec.originX) / spec.dx));
    const std::size_t center_idx = grid.idx(i_center, j_center);
    const double Hy_center = grid.Hy[center_idx];
    if (!(std::abs(Hy_center) < 3.0e5)) {
        std::cerr << "magnet_strip_test: Hy inside magnet expected near zero, got " << Hy_center << '\n';
        return 1;
    }

    std::cout << "magnet_strip_test: relErr=" << relErr << ", Hy_center=" << Hy_center << '\n';
    return 0;
}
