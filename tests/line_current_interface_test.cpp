#include "motorsim/ingest.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <array>
#include <cmath>
#include <filesystem>
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

ErrorSummary run_case(const motorsim::ScenarioSpec& baseSpec, double mu_r_right,
                      motorsim::ScenarioSpec::BoundaryType boundary, double sample_offset,
                      double max_probe_y, motorsim::SolverKind solverKind) {
    using namespace motorsim;

    ScenarioSpec spec = baseSpec;
    spec.boundaryType = boundary;

    if (spec.wires.size() != 1) {
        throw std::runtime_error("line_current_interface_test: expected exactly one wire in scenario");
    }

    const ScenarioSpec::Wire wire = spec.wires.front();

    std::size_t leftIndex = std::numeric_limits<std::size_t>::max();
    std::size_t rightIndex = std::numeric_limits<std::size_t>::max();
    for (std::size_t idx = 0; idx < spec.halfspaces.size(); ++idx) {
        const auto& hs = spec.halfspaces[idx];
        if (hs.normal_x > 0.5 && leftIndex == std::numeric_limits<std::size_t>::max()) {
            leftIndex = idx;
        } else if (hs.normal_x < -0.5 && rightIndex == std::numeric_limits<std::size_t>::max()) {
            rightIndex = idx;
        }
    }

    if (leftIndex == std::numeric_limits<std::size_t>::max() ||
        rightIndex == std::numeric_limits<std::size_t>::max()) {
        throw std::runtime_error("line_current_interface_test: scenario is missing expected halfspaces");
    }

    const double mu_r_left = spec.halfspaces[leftIndex].mu_r;
    spec.halfspaces[rightIndex].mu_r = mu_r_right;
    spec.halfspaces[rightIndex].inv_mu = 1.0 / (MU0 * mu_r_right);

    Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    rasterizeScenarioToGrid(spec, grid);

    SolveOptions options{};
    options.maxIters = 25000;
    options.tol = 5e-6;
    options.omega = 1.6;
    options.kind = solverKind;

    const SolveResult report = solveAz(grid, options);
    if (!report.converged || !(report.relResidual < options.tol)) {
        std::cerr << "line_current_interface_test: warning relResidual=" << report.relResidual
                  << " (tol=" << options.tol << ")\n";
    }

    computeB(grid);

    const double mu1 = MU0 * mu_r_left;
    const double mu2 = MU0 * mu_r_right;
    const double rho = (mu2 - mu1) / (mu1 + mu2);
    const double tau = (2.0 * mu1) / (mu1 + mu2);

    const std::size_t i_left = static_cast<std::size_t>(
        std::llround(((-sample_offset) - spec.originX) / spec.dx));
    const std::size_t i_right = static_cast<std::size_t>(
        std::llround(((sample_offset) - spec.originX) / spec.dx));
    if (i_left >= spec.nx || i_right >= spec.nx) {
        throw std::runtime_error("line_current_interface_test: probe indices out of range");
    }

    const auto accumulate_error = [&](std::size_t column_index, bool left_side) {
        double err2 = 0.0;
        double ref2 = 0.0;
        std::size_t samples = 0;
        for (std::size_t j = 2; j + 2 < spec.ny; ++j) {
            const double y = spec.originY + static_cast<double>(j) * spec.dy;
            if (std::abs(y) > max_probe_y) {
                continue;
            }
            const std::size_t idx = grid.idx(column_index, j);
            const double bx_num = grid.Bx[idx];
            const double by_num = grid.By[idx];

            double bx_ana = 0.0;
            double by_ana = 0.0;
            const double sample_x = spec.originX + static_cast<double>(column_index) * spec.dx;
            if (left_side) {
                const auto real_field = B_from_line(sample_x, y, wire.x, wire.y, wire.current, mu1);
                const auto image_field =
                    B_from_line(sample_x, y, -wire.x, wire.y, rho * wire.current, mu1);
                bx_ana = real_field.first + image_field.first;
                by_ana = real_field.second + image_field.second;
            } else {
                const auto transmitted =
                    B_from_line(sample_x, y, wire.x, wire.y, tau * wire.current, mu2);
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
    namespace fs = std::filesystem;

    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/line_current_interface_test.json").lexically_normal();

    motorsim::ScenarioSpec baseSpec = motorsim::loadScenarioFromJson(scenarioPath.string());

    const double sample_offset = 0.08;
    const double max_probe_y = 0.06;

    const struct {
        double mu_r_right;
        double leftTol;
        double rightTol;
    } dirichlet_cases[] = {
        {100.0, 0.20, 0.70},
    };

    ErrorSummary referenceDirichletSor{};
    for (const auto& entry : dirichlet_cases) {
        const ErrorSummary sorSummary =
            run_case(baseSpec, entry.mu_r_right, ScenarioSpec::BoundaryType::Dirichlet, sample_offset,
                     max_probe_y, motorsim::SolverKind::SOR);
        const ErrorSummary cgSummary =
            run_case(baseSpec, entry.mu_r_right, ScenarioSpec::BoundaryType::Dirichlet, sample_offset,
                     max_probe_y, motorsim::SolverKind::CG);

        if (!(sorSummary.relResidual < 2e-1)) {
            std::cerr << "line_current_interface_test: Dirichlet solver residual=" << sorSummary.relResidual
                      << " exceeds guard 0.2\n";
            return 1;
        }
        if (!(sorSummary.leftRelErr < entry.leftTol)) {
            std::cerr << "line_current_interface_test: Dirichlet left relErr=" << sorSummary.leftRelErr
                      << " exceeds tol=" << entry.leftTol << " for mu_r_right=" << entry.mu_r_right << "\n";
            return 1;
        }
        if (!(sorSummary.rightRelErr < entry.rightTol)) {
            std::cerr << "line_current_interface_test: Dirichlet right relErr=" << sorSummary.rightRelErr
                      << " exceeds tol=" << entry.rightTol << " for mu_r_right=" << entry.mu_r_right << "\n";
            return 1;
        }
        if (!(cgSummary.relResidual <= sorSummary.relResidual * 1.1 + 1e-12)) {
            std::cerr << "line_current_interface_test: Dirichlet CG residual=" << cgSummary.relResidual
                      << " exceeds parity bound relative to SOR " << sorSummary.relResidual << "\n";
            return 1;
        }
        if (!(cgSummary.leftRelErr <= sorSummary.leftRelErr + 0.02)) {
            std::cerr << "line_current_interface_test: Dirichlet CG left relErr=" << cgSummary.leftRelErr
                      << " exceeds parity slack relative to SOR " << sorSummary.leftRelErr << "\n";
            return 1;
        }
        if (!(cgSummary.rightRelErr <= sorSummary.rightRelErr + 0.02)) {
            std::cerr << "line_current_interface_test: Dirichlet CG right relErr=" << cgSummary.rightRelErr
                      << " exceeds parity slack relative to SOR " << sorSummary.rightRelErr << "\n";
            return 1;
        }

        referenceDirichletSor = sorSummary;
        std::cout << "line_current_interface_test: mu_r_right=" << entry.mu_r_right
                  << ", boundary=Dirichlet (SOR) leftRelErr=" << sorSummary.leftRelErr
                  << ", rightRelErr=" << sorSummary.rightRelErr << '\n';
        std::cout << "line_current_interface_test: mu_r_right=" << entry.mu_r_right
                  << ", boundary=Dirichlet (CG) leftRelErr=" << cgSummary.leftRelErr
                  << ", rightRelErr=" << cgSummary.rightRelErr << '\n';
    }

    const ErrorSummary neumannSor =
        run_case(baseSpec, 100.0, ScenarioSpec::BoundaryType::Neumann, sample_offset, max_probe_y,
                 motorsim::SolverKind::SOR);
    const ErrorSummary neumannCg =
        run_case(baseSpec, 100.0, ScenarioSpec::BoundaryType::Neumann, sample_offset, max_probe_y,
                 motorsim::SolverKind::CG);
    std::cout << "line_current_interface_test: mu_r_right=100.0, boundary=Neumann (SOR) leftRelErr="
              << neumannSor.leftRelErr << ", rightRelErr=" << neumannSor.rightRelErr << '\n';
    std::cout << "line_current_interface_test: mu_r_right=100.0, boundary=Neumann (CG) leftRelErr="
              << neumannCg.leftRelErr << ", rightRelErr=" << neumannCg.rightRelErr << '\n';

    if (!(neumannSor.relResidual < 2e-1)) {
        std::cerr << "line_current_interface_test: Neumann solver residual=" << neumannSor.relResidual
                  << " exceeds guard 0.2\n";
        return 1;
    }
    if (!(neumannSor.leftRelErr < 0.25)) {
        std::cerr << "line_current_interface_test: Neumann left error " << neumannSor.leftRelErr
                  << " exceeds tolerance 0.25\n";
        return 1;
    }
    if (!(neumannSor.rightRelErr < referenceDirichletSor.rightRelErr)) {
        std::cerr << "line_current_interface_test: Neumann right error did not improve ("
                  << neumannSor.rightRelErr << " vs " << referenceDirichletSor.rightRelErr << ")\n";
        return 1;
    }
    if (!(neumannCg.relResidual <= neumannSor.relResidual * 1.1 + 1e-12)) {
        std::cerr << "line_current_interface_test: Neumann CG residual=" << neumannCg.relResidual
                  << " exceeds parity bound relative to SOR " << neumannSor.relResidual << "\n";
        return 1;
    }
    if (!(neumannCg.leftRelErr <= neumannSor.leftRelErr + 0.02)) {
        std::cerr << "line_current_interface_test: Neumann CG left relErr=" << neumannCg.leftRelErr
                  << " exceeds parity slack relative to SOR " << neumannSor.leftRelErr << "\n";
        return 1;
    }
    if (!(neumannCg.rightRelErr <= neumannSor.rightRelErr + 0.02)) {
        std::cerr << "line_current_interface_test: Neumann CG right relErr=" << neumannCg.rightRelErr
                  << " exceeds parity slack relative to SOR " << neumannSor.rightRelErr << "\n";
        return 1;
    }

    return 0;
}
