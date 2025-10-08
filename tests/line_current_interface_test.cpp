#include "motorsim/ingest.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <utility>

namespace {
std::pair<double, double> B_from_line(double x, double y, double xs, double ys,
                                      double current, double mu) {
    const double dx = x - xs;
    const double dy = y - ys;
    const double r2 = dx * dx + dy * dy;
    if (r2 <= 0.0) {
        return {0.0, 0.0};
    }
    const double coeff = (mu * current) / (2.0 * M_PI * r2);
    return {-coeff * dy, coeff * dx};
}

}  // namespace

int main() {
    using namespace motorsim;

    const double Lx = 0.30;
    const double Ly = 0.20;
    const std::size_t nx = 241;
    const std::size_t ny = 161;
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

    const double mu_r_left = 1.0;
    const double mu_r_right = 1000.0;
    spec.mu_r_background = mu_r_left;
    ScenarioSpec::HalfspaceRegion left_region{};
    left_region.normal_x = 1.0;
    left_region.normal_y = 0.0;
    left_region.offset = 0.0;
    left_region.mu_r = mu_r_left;
    left_region.inv_mu = 1.0 / (MU0 * mu_r_left);
    spec.halfspaces.push_back(left_region);  // x < 0
    ScenarioSpec::RegionMask left_mask{};
    left_mask.kind = ScenarioSpec::RegionMask::Kind::Halfspace;
    left_mask.index = 0;
    spec.regionMasks.push_back(left_mask);

    ScenarioSpec::HalfspaceRegion right_region{};
    right_region.normal_x = -1.0;
    right_region.normal_y = 0.0;
    right_region.offset = 0.0;
    right_region.mu_r = mu_r_right;
    right_region.inv_mu = 1.0 / (MU0 * mu_r_right);
    spec.halfspaces.push_back(right_region);  // x > 0
    ScenarioSpec::RegionMask right_mask{};
    right_mask.kind = ScenarioSpec::RegionMask::Kind::Halfspace;
    right_mask.index = 1;
    spec.regionMasks.push_back(right_mask);

    const double current = 10.0;
    const double wire_x = -0.06;
    const double wire_y = 0.0;
    const double wire_radius = 3.0 * dx;
    spec.wires.push_back({wire_x, wire_y, wire_radius, current});

    Grid2D grid(nx, ny, dx, dy);
    rasterizeScenarioToGrid(spec, grid);

    SolveOptions options{};
    options.maxIters = 10000;
    options.tol = 1e-6;
    options.omega = 1.7;

    const SolveReport report = solveAz_GS_SOR(grid, options);
    if (!report.converged) {
        std::cerr << "line_current_interface_test: solver failed to converge\n";
        return 1;
    }
    if (!(report.relResidual < options.tol)) {
        std::cerr << "line_current_interface_test: residual " << report.relResidual
                  << " exceeds tolerance " << options.tol << "\n";
        return 1;
    }

    computeB(grid);

    const double mu1 = MU0 * mu_r_left;
    const double mu2 = MU0 * mu_r_right;
    const double rho = (mu1 - mu2) / (mu1 + mu2);
    const double tau = (2.0 * mu2) / (mu1 + mu2);

    const double sample_offset = 0.04;
    const std::size_t i_left = static_cast<std::size_t>(std::llround(((-sample_offset) - spec.originX) / spec.dx));
    const std::size_t i_right = static_cast<std::size_t>(std::llround(((sample_offset) - spec.originX) / spec.dx));

    if (i_left >= nx || i_right >= nx) {
        std::cerr << "line_current_interface_test: probe indices out of range\n";
        return 1;
    }

    const auto accumulate_error = [&](std::size_t column_index, bool left_side) {
        const double x = spec.originX + static_cast<double>(column_index) * spec.dx;
        double err2 = 0.0;
        double ref2 = 0.0;
        std::size_t samples = 0;
        for (std::size_t j = 2; j + 2 < ny; ++j) {
            const double y = spec.originY + static_cast<double>(j) * spec.dy;
            if (std::abs(y) > 0.08) {
                continue;
            }
            const std::size_t idx = grid.idx(column_index, j);
            const double bx_num = grid.Bx[idx];
            const double by_num = grid.By[idx];

            double bx_ana = 0.0;
            double by_ana = 0.0;
            if (left_side) {
                const auto real_field = B_from_line(x, y, wire_x, wire_y, current, mu1);
                const auto image_field = B_from_line(x, y, -wire_x, wire_y, rho * current, mu1);
                bx_ana = real_field.first + image_field.first;
                by_ana = real_field.second + image_field.second;
            } else {
                const auto transmitted = B_from_line(x, y, wire_x, wire_y, tau * current, mu2);
                bx_ana = transmitted.first;
                by_ana = transmitted.second;
            }

            const double bnum_mag = std::hypot(bx_num, by_num);
            const double bana_mag = std::hypot(bx_ana, by_ana);
            if (bana_mag < 1e-9) {
                continue;
            }
            const double diff = bnum_mag - bana_mag;
            err2 += diff * diff;
            ref2 += bana_mag * bana_mag;
            ++samples;
        }
        if (samples == 0 || ref2 <= 0.0) {
            return std::make_pair(std::numeric_limits<double>::quiet_NaN(), samples);
        }
        return std::make_pair(std::sqrt(err2 / ref2), samples);
    };

    const auto left_result = accumulate_error(i_left, true);
    const auto right_result = accumulate_error(i_right, false);

    if (!(left_result.second > ny / 4)) {
        std::cerr << "line_current_interface_test: insufficient left samples\n";
        return 1;
    }
    if (!(right_result.second > ny / 4)) {
        std::cerr << "line_current_interface_test: insufficient right samples\n";
        return 1;
    }

    if (!(left_result.first < 0.40)) {
        std::cerr << "line_current_interface_test: left relErr=" << left_result.first << "\n";
        return 1;
    }
    std::cout << "line_current_interface_test: relResidual=" << report.relResidual
              << ", leftRelErr=" << left_result.first
              << ", rightRelErr=" << right_result.first << "\n";
    return 0;
}
