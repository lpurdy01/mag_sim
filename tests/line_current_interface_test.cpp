#include "motorsim/ingest.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
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

struct ErrorSummary {
    double relResidual{0.0};
    double leftRelErr{0.0};
    double rightRelErr{0.0};
};

ErrorSummary run_case(double mu_r_right, motorsim::ScenarioSpec::BoundaryType boundary,
                      double sample_offset, double max_probe_y) {
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
    spec.boundaryType = boundary;

    const double mu_r_left = 1.0;
    spec.mu_r_background = mu_r_left;
    ScenarioSpec::HalfspaceRegion left_region{};
    left_region.normal_x = 1.0;
    left_region.normal_y = 0.0;
    left_region.offset = 0.0;
    left_region.mu_r = mu_r_left;
    left_region.inv_mu = 1.0 / (MU0 * mu_r_left);
    spec.halfspaces.push_back(left_region);
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
    spec.halfspaces.push_back(right_region);
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
    options.maxIters = 25000;
    options.tol = 5e-6;
    options.omega = 1.6;

    const SolveReport report = solveAz_GS_SOR(grid, options);
    if (!report.converged || !(report.relResidual < options.tol)) {
        std::cerr << "line_current_interface_test: warning relResidual=" << report.relResidual
                  << " (tol=" << options.tol << ")\n";
    }

    computeB(grid);

    const double mu1 = MU0 * mu_r_left;
    const double mu2 = MU0 * mu_r_right;
    const double rho = (mu1 - mu2) / (mu1 + mu2);
    const double tau = (2.0 * mu2) / (mu1 + mu2);

    const std::size_t i_left = static_cast<std::size_t>(
        std::llround(((-sample_offset) - spec.originX) / spec.dx));
    const std::size_t i_right = static_cast<std::size_t>(
        std::llround(((sample_offset) - spec.originX) / spec.dx));
    if (i_left >= nx || i_right >= nx) {
        throw std::runtime_error("line_current_interface_test: probe indices out of range");
    }

    const auto accumulate_error = [&](std::size_t column_index, bool left_side) {
        double err2 = 0.0;
        double ref2 = 0.0;
        std::size_t samples = 0;
        for (std::size_t j = 2; j + 2 < ny; ++j) {
            const double y = spec.originY + static_cast<double>(j) * spec.dy;
            if (std::abs(y) > max_probe_y) {
                continue;
            }
            const std::size_t idx = grid.idx(column_index, j);
            const double bx_num = grid.Bx[idx];
            const double by_num = grid.By[idx];

            double bx_ana = 0.0;
            double by_ana = 0.0;
            if (left_side) {
                const auto real_field = B_from_line(column_index * spec.dx + spec.originX, y, wire_x, wire_y, current, mu1);
                const auto image_field = B_from_line(column_index * spec.dx + spec.originX, y, -wire_x, wire_y,
                                                     rho * current, mu1);
                bx_ana = real_field.first + image_field.first;
                by_ana = real_field.second + image_field.second;
            } else {
                const auto transmitted =
                    B_from_line(column_index * spec.dx + spec.originX, y, wire_x, wire_y, tau * current, mu2);
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
            throw std::runtime_error("line_current_interface_test: insufficient samples");
        }
        return std::sqrt(err2 / ref2);
    };

    ErrorSummary summary{};
    summary.relResidual = report.relResidual;
    summary.leftRelErr = accumulate_error(i_left, true);
    summary.rightRelErr = accumulate_error(i_right, false);
    return summary;
}

int main() {
    using motorsim::ScenarioSpec;

    const double sample_offset = 0.08;
    const double max_probe_y = 0.06;

    const struct {
        double mu_r_right;
        double leftTol;
        double rightTol;
    } dirichlet_cases[] = {
        {100.0, 0.40, 1.05},
    };

    ErrorSummary referenceDirichlet{};
    for (const auto& entry : dirichlet_cases) {
        const ErrorSummary summary =
            run_case(entry.mu_r_right, ScenarioSpec::BoundaryType::Dirichlet, sample_offset, max_probe_y);
        if (!(summary.relResidual < 2e-1)) {
            std::cerr << "line_current_interface_test: Dirichlet solver residual=" << summary.relResidual
                      << " exceeds guard 0.2\n";
            return 1;
        }
        if (!(summary.leftRelErr < entry.leftTol)) {
            std::cerr << "line_current_interface_test: Dirichlet left relErr=" << summary.leftRelErr
                      << " exceeds tol=" << entry.leftTol << " for mu_r_right=" << entry.mu_r_right << "\n";
            return 1;
        }
        if (!(summary.rightRelErr < entry.rightTol)) {
            std::cerr << "line_current_interface_test: Dirichlet right relErr=" << summary.rightRelErr
                      << " exceeds tol=" << entry.rightTol << " for mu_r_right=" << entry.mu_r_right << "\n";
            return 1;
        }
        referenceDirichlet = summary;
        std::cout << "line_current_interface_test: mu_r_right=" << entry.mu_r_right
                  << ", boundary=Dirichlet, leftRelErr=" << summary.leftRelErr
                  << ", rightRelErr=" << summary.rightRelErr << '\n';
    }

    const ErrorSummary neumann_summary =
        run_case(100.0, ScenarioSpec::BoundaryType::Neumann, sample_offset, max_probe_y);
    if (!(neumann_summary.relResidual < 2e-1)) {
        std::cerr << "line_current_interface_test: Neumann solver residual=" << neumann_summary.relResidual
                  << " exceeds guard 0.2\n";
        return 1;
    }
    if (!(neumann_summary.leftRelErr <= referenceDirichlet.leftRelErr)) {
        std::cerr << "line_current_interface_test: Neumann left error did not improve ("
                  << neumann_summary.leftRelErr << " vs " << referenceDirichlet.leftRelErr << ")\n";
        return 1;
    }
    if (!(neumann_summary.rightRelErr < referenceDirichlet.rightRelErr)) {
        std::cerr << "line_current_interface_test: Neumann right error did not improve ("
                  << neumann_summary.rightRelErr << " vs " << referenceDirichlet.rightRelErr << ")\n";
        return 1;
    }
    std::cout << "line_current_interface_test: mu_r_right=100.0, boundary=Neumann, leftRelErr="
              << neumann_summary.leftRelErr << ", rightRelErr=" << neumann_summary.rightRelErr << '\n';

    return 0;
}
