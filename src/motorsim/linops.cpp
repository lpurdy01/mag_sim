#include "motorsim/linops.hpp"

#include <algorithm>
#include <cmath>

namespace motorsim::linops {
namespace {

#if defined(MOTORSIM_SIMD_HINTS)
#    if defined(__clang__)
#        define MOTORSIM_LOOP_HINT _Pragma("clang loop vectorize(enable) interleave(enable)")
#    elif defined(__GNUC__)
#        define MOTORSIM_LOOP_HINT _Pragma("GCC ivdep")
#    elif defined(_MSC_VER)
#        define MOTORSIM_LOOP_HINT __pragma(loop(ivdep))
#    else
#        define MOTORSIM_LOOP_HINT
#    endif
#else
#    define MOTORSIM_LOOP_HINT
#endif

[[nodiscard]] double harmonicAverage(double a, double b) {
    const double denom = a + b;
    if (denom == 0.0) {
        return 0.0;
    }
    return 2.0 * a * b / denom;
}

}  // namespace

void zero(std::size_t n, double* x) {
    if (n == 0 || x == nullptr) {
        return;
    }
    MOTORSIM_LOOP_HINT
    for (std::size_t idx = 0; idx < n; ++idx) {
        x[idx] = 0.0;
    }
}

void applyA(const Grid2D& grid, const double* x, double* y) {
    const std::size_t n = grid.nx * grid.ny;
    if (n == 0 || x == nullptr || y == nullptr) {
        return;
    }

    zero(n, y);

    if (grid.nx < 3 || grid.ny < 3) {
        return;
    }

    const double invDx2 = 1.0 / (grid.dx * grid.dx);
    const double invDy2 = 1.0 / (grid.dy * grid.dy);

    for (std::size_t j = 1; j + 1 < grid.ny; ++j) {
        for (std::size_t i = 1; i + 1 < grid.nx; ++i) {
            const std::size_t p = grid.idx(i, j);
            const std::size_t e = grid.idx(i + 1, j);
            const std::size_t w = grid.idx(i - 1, j);
            const std::size_t nIdx = grid.idx(i, j + 1);
            const std::size_t s = grid.idx(i, j - 1);

            const double nuP = grid.invMu[p];
            const double nuE = harmonicAverage(nuP, grid.invMu[e]);
            const double nuW = harmonicAverage(nuP, grid.invMu[w]);
            const double nuN = harmonicAverage(nuP, grid.invMu[nIdx]);
            const double nuS = harmonicAverage(nuP, grid.invMu[s]);

            double azCenter = x[p];
            double azEast = x[e];
            double azWest = x[w];
            double azNorth = x[nIdx];
            double azSouth = x[s];

            if (grid.boundaryCondition == Grid2D::BoundaryKind::Neumann) {
                if (i == 1) {
                    azWest = azCenter;
                }
                if (i + 2 == grid.nx) {
                    azEast = azCenter;
                }
                if (j == 1) {
                    azSouth = azCenter;
                }
                if (j + 2 == grid.ny) {
                    azNorth = azCenter;
                }
            }

            const double flux = (nuE * (azEast - azCenter) - nuW * (azCenter - azWest)) * invDx2 +
                                (nuN * (azNorth - azCenter) - nuS * (azCenter - azSouth)) * invDy2;
            y[p] = -flux;
        }
    }

    if (grid.boundaryCondition == Grid2D::BoundaryKind::Dirichlet) {
        for (std::size_t i = 0; i < grid.nx; ++i) {
            y[grid.idx(i, 0)] = x[grid.idx(i, 0)];
            y[grid.idx(i, grid.ny - 1)] = x[grid.idx(i, grid.ny - 1)];
        }
        for (std::size_t j = 0; j < grid.ny; ++j) {
            y[grid.idx(0, j)] = x[grid.idx(0, j)];
            y[grid.idx(grid.nx - 1, j)] = x[grid.idx(grid.nx - 1, j)];
        }
    }
}

void computeResidual(const Grid2D& grid, const double* x, double* residual) {
    const std::size_t n = grid.nx * grid.ny;
    if (n == 0 || residual == nullptr) {
        return;
    }
    applyA(grid, x, residual);

    if (grid.Jz.size() < n) {
        return;
    }

    const double* rhs = grid.Jz.data();
    MOTORSIM_LOOP_HINT
    for (std::size_t idx = 0; idx < n; ++idx) {
        residual[idx] -= rhs[idx];
    }
}

double dot(std::size_t n, const double* a, const double* b) {
    if (n == 0 || a == nullptr || b == nullptr) {
        return 0.0;
    }
    double accum = 0.0;
    MOTORSIM_LOOP_HINT
    for (std::size_t idx = 0; idx < n; ++idx) {
        accum += a[idx] * b[idx];
    }
    return accum;
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    const std::size_t n = std::min(a.size(), b.size());
    return dot(n, a.data(), b.data());
}

double squaredNorm(std::size_t n, const double* x) {
    if (n == 0 || x == nullptr) {
        return 0.0;
    }
    double accum = 0.0;
    MOTORSIM_LOOP_HINT
    for (std::size_t idx = 0; idx < n; ++idx) {
        accum += x[idx] * x[idx];
    }
    return accum;
}

double squaredNorm(const std::vector<double>& x) {
    return squaredNorm(x.size(), x.data());
}

double norm(std::size_t n, const double* x) {
    return std::sqrt(squaredNorm(n, x));
}

double norm(const std::vector<double>& x) {
    return std::sqrt(squaredNorm(x));
}

double residualNorm(const Grid2D& grid, const double* x, std::vector<double>& scratch) {
    const std::size_t n = grid.nx * grid.ny;
    scratch.resize(n);
    computeResidual(grid, x, scratch.data());
    return norm(n, scratch.data());
}

void axpy(std::size_t n, double alpha, const double* x, double* y) {
    if (n == 0 || x == nullptr || y == nullptr) {
        return;
    }
    MOTORSIM_LOOP_HINT
    for (std::size_t idx = 0; idx < n; ++idx) {
        y[idx] += alpha * x[idx];
    }
}

void scal(std::size_t n, double alpha, double* x) {
    if (n == 0 || x == nullptr) {
        return;
    }
    MOTORSIM_LOOP_HINT
    for (std::size_t idx = 0; idx < n; ++idx) {
        x[idx] *= alpha;
    }
}

}  // namespace motorsim::linops
