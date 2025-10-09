// filename: solver.cpp
// part of 2D Electromagnetic Motor Simulator
// MIT License

#include "motorsim/solver.hpp"

#include "motorsim/types.hpp"

#include <cmath>
#include <cstddef>
#include <iostream>

namespace motorsim {
namespace {

[[nodiscard]] double harmonicAverage(double a, double b) {
    const double denom = a + b;
    if (denom == 0.0) {
        return 0.0;
    }
    return 2.0 * a * b / denom;
}

}  // namespace

SolveReport solveAz_GS_SOR(Grid2D& grid, const SolveOptions& options) {
    SolveReport report{};
    if (grid.nx < 3 || grid.ny < 3) {
        return report;
    }

    const double invDx2 = 1.0 / (grid.dx * grid.dx);
    const double invDy2 = 1.0 / (grid.dy * grid.dy);

    double sumJ2 = 0.0;
    for (double value : grid.Jz) {
        sumJ2 += value * value;
    }
    const double denom = std::sqrt(sumJ2 + 1e-30);

    for (std::size_t iter = 0; iter < options.maxIters; ++iter) {
        double residualNorm2 = 0.0;

        for (std::size_t j = 1; j + 1 < grid.ny; ++j) {
            for (std::size_t i = 1; i + 1 < grid.nx; ++i) {
                const std::size_t p = grid.idx(i, j);
                const std::size_t e = grid.idx(i + 1, j);
                const std::size_t w = grid.idx(i - 1, j);
                const std::size_t n = grid.idx(i, j + 1);
                const std::size_t s = grid.idx(i, j - 1);

                const double nuP = grid.invMu[p];
                const double nuE = harmonicAverage(nuP, grid.invMu[e]);
                const double nuW = harmonicAverage(nuP, grid.invMu[w]);
                const double nuN = harmonicAverage(nuP, grid.invMu[n]);
                const double nuS = harmonicAverage(nuP, grid.invMu[s]);

                // 5-point stencil denominator assembled from harmonic-averaged face coefficients
                const double diag = (nuE + nuW) * invDx2 + (nuN + nuS) * invDy2;
                if (diag == 0.0) {
                    continue;
                }

                // Discrete Poisson update with variable coefficients and source term Jz
                const double numerator = (nuE * grid.Az[e] + nuW * grid.Az[w]) * invDx2 +
                                         (nuN * grid.Az[n] + nuS * grid.Az[s]) * invDy2 +
                                         grid.Jz[p];
                const double candidate = numerator / diag;
                const double updated = (1.0 - options.omega) * grid.Az[p] + options.omega * candidate;
                grid.Az[p] = updated;

                if (grid.boundaryCondition == Grid2D::BoundaryKind::Neumann) {
                    if (i == 1) {
                        grid.Az[w] = updated;
                    }
                    if (i + 2 == grid.nx) {
                        grid.Az[e] = updated;
                    }
                    if (j == 1) {
                        grid.Az[s] = updated;
                    }
                    if (j + 2 == grid.ny) {
                        grid.Az[n] = updated;
                    }
                }

                // Residual: divergence(nu grad A) + J should approach zero
                const double flux = (nuE * (grid.Az[e] - updated) - nuW * (updated - grid.Az[w])) * invDx2 +
                                    (nuN * (grid.Az[n] - updated) - nuS * (updated - grid.Az[s])) * invDy2;
                const double cellResidual = flux + grid.Jz[p];
                residualNorm2 += cellResidual * cellResidual;
            }
        }

        const double relResidual = std::sqrt(residualNorm2) / denom;
        report.relResidual = relResidual;
        report.iters = iter + 1;
        if (options.verbose) {
            std::cout << "Iteration " << report.iters << ": relResidual=" << relResidual << '\n';
        }
        if (relResidual < options.tol) {
            report.converged = true;
            break;
        }
    }

    return report;
}

void computeB(Grid2D& grid) {
    const std::size_t count = grid.nx * grid.ny;
    if (grid.Bx.size() != count) {
        grid.Bx.assign(count, 0.0);
    }
    if (grid.By.size() != count) {
        grid.By.assign(count, 0.0);
    }

    const double inv2dx = 1.0 / (2.0 * grid.dx);
    const double inv2dy = 1.0 / (2.0 * grid.dy);

    for (std::size_t j = 0; j < grid.ny; ++j) {
        for (std::size_t i = 0; i < grid.nx; ++i) {
            const std::size_t p = grid.idx(i, j);

            double bx = 0.0;
            if (j == 0) {
                bx = (grid.Az[grid.idx(i, j + 1)] - grid.Az[p]) / grid.dy;
            } else if (j + 1 == grid.ny) {
                bx = (grid.Az[p] - grid.Az[grid.idx(i, j - 1)]) / grid.dy;
            } else {
                bx = (grid.Az[grid.idx(i, j + 1)] - grid.Az[grid.idx(i, j - 1)]) * inv2dy;
            }

            double by = 0.0;
            if (i == 0) {
                by = -(grid.Az[grid.idx(i + 1, j)] - grid.Az[p]) / grid.dx;
            } else if (i + 1 == grid.nx) {
                by = -(grid.Az[p] - grid.Az[grid.idx(i - 1, j)]) / grid.dx;
            } else {
                by = -(grid.Az[grid.idx(i + 1, j)] - grid.Az[grid.idx(i - 1, j)]) * inv2dx;
            }

            grid.Bx[p] = bx;
            grid.By[p] = by;
        }
    }
}

void computeH(Grid2D& grid) {
    const std::size_t count = grid.nx * grid.ny;
    if (grid.Hx.size() != count) {
        grid.Hx.assign(count, 0.0);
    }
    if (grid.Hy.size() != count) {
        grid.Hy.assign(count, 0.0);
    }

    for (std::size_t idx = 0; idx < count; ++idx) {
        double hx = grid.invMu[idx] * grid.Bx[idx];
        double hy = grid.invMu[idx] * grid.By[idx];
        if (grid.Mx.size() == count && grid.My.size() == count) {
            const double invMuCell = grid.invMu[idx];
            if (invMuCell > 0.0) {
                const double mu_r = 1.0 / (MU0 * invMuCell);
                hx -= grid.Mx[idx] / mu_r;
                hy -= grid.My[idx] / mu_r;
            } else {
                hx -= grid.Mx[idx];
                hy -= grid.My[idx];
            }
        }
        grid.Hx[idx] = hx;
        grid.Hy[idx] = hy;
    }
}

SolveReport solveAz_CG(Grid2D& grid, const SolveOptions& options) {
    (void)grid;
    (void)options;
    SolveReport report{};
    report.converged = false;
    report.relResidual = 1.0;
    return report;
}

}  // namespace motorsim
