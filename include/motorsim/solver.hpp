// filename: solver.hpp
// part of 2D Electromagnetic Motor Simulator
// MIT License

#pragma once

#include "grid.hpp"

#include <chrono>
#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace motorsim {

/**
 * @brief Snapshot of the iterative solver's progress.
 */
enum class SolverKind { SOR, CG, Harmonic };

struct InitialGuess {
    const std::vector<double>* Az0{nullptr};
    const std::vector<double>* AzImag0{nullptr};
};

struct ProgressSample {
    std::size_t iter{0};
    double relResidual{0.0};
    double elapsedSeconds{0.0};
};

struct ProgressSink {
    virtual ~ProgressSink() = default;
    virtual bool onProgress(const ProgressSample& sample) = 0;
    virtual std::optional<std::string> requestFieldDump(const ProgressSample& sample) {
        (void)sample;
        return std::nullopt;
    }
};

struct SolveOptions {
    SolverKind kind{SolverKind::SOR};
    std::size_t maxIters{20000};
    double tol{1e-6};
    double omega{1.7};
    std::size_t snapshotEveryIters{0};
    double progressEverySec{2.0};
    bool warmStart{false};
    bool useProlongation{false};
    std::size_t coarseNx{0};
    std::size_t coarseNy{0};
    bool verbose{false};
    double harmonicOmega{0.0};
};

struct SolveResult {
    std::size_t iters{0};
    double relResidual{0.0};
    bool converged{false};
};

SolveResult solveAz(Grid2D& grid,
                    const SolveOptions& options,
                    const InitialGuess* initialGuess = nullptr,
                    ProgressSink* progress = nullptr);

/**
 * @brief Compute the harmonic magnetic flux density components from complex vector potential parts.
 */
void computeBHarmonic(const Grid2D& grid,
                      const std::vector<double>& AzReal,
                      const std::vector<double>& AzImag,
                      std::vector<double>& BxReal,
                      std::vector<double>& BxImag,
                      std::vector<double>& ByReal,
                      std::vector<double>& ByImag);

void computeHHarmonic(const Grid2D& grid,
                      const std::vector<double>& BxReal,
                      const std::vector<double>& BxImag,
                      const std::vector<double>& ByReal,
                      const std::vector<double>& ByImag,
                      std::vector<double>& HxReal,
                      std::vector<double>& HxImag,
                      std::vector<double>& HyReal,
                      std::vector<double>& HyImag);

SolveResult solveAz_GS_SOR(Grid2D& grid,
                           const SolveOptions& options,
                           const InitialGuess* initialGuess = nullptr,
                           ProgressSink* progress = nullptr);

SolveResult solveAz_CG(Grid2D& grid,
                       const SolveOptions& options,
                       const InitialGuess* initialGuess = nullptr,
                       ProgressSink* progress = nullptr);

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

}  // namespace motorsim
