// filename: solver.hpp
// part of 2D Electromagnetic Motor Simulator
// MIT License

#pragma once

#include "grid.hpp"

#include <chrono>
#include <functional>

namespace motorsim {

/**
 * @brief Snapshot of the iterative solver's progress.
 */
struct SolverProgress {
    std::size_t iteration{0};
    double relResidual{0.0};
};

/**
 * @brief Solver options controlling the Gauss–Seidel / SOR iteration.
 */
struct SolveOptions {
    std::size_t maxIters{2000};
    double tol{1e-6};
    double omega{1.7};
    bool verbose{false};
    std::chrono::steady_clock::duration progressInterval{std::chrono::steady_clock::duration::zero()};
    std::function<void(const SolverProgress&)> progressCallback{};
};

/**
 * @brief Reports the outcome of the magnetostatic solve.
 */
struct SolveReport {
    std::size_t iters{0};
    double relResidual{0.0};
    bool converged{false};
};

SolveReport solveAz_GS_SOR(Grid2D& grid, const SolveOptions& options);

/**
 * @brief Compute magnetic flux density components from the solved vector potential.
 * @param grid Grid storing the magnetic vector potential; filled with Bx and By on return.
 */
void computeB(Grid2D& grid);

/**
 * @brief Compute magnetic field intensity components using the stored flux density
 *        and magnetisation.
 * @param grid Grid with populated Bx/By and magnetisation vectors.
 */
void computeH(Grid2D& grid);

/**
 * @brief Placeholder for a future conjugate gradient solver.
 * @note Currently unimplemented; returns a non-converged report.
 */
SolveReport solveAz_CG(Grid2D& grid, const SolveOptions& options);

}  // namespace motorsim
