#include "motorsim/solver.hpp"

#include "motorsim/grid.hpp"
#include "motorsim/io_csv.hpp"
#include "motorsim/linops.hpp"
#include "motorsim/types.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace motorsim {
namespace {

using Clock = std::chrono::steady_clock;

[[nodiscard]] double harmonicAverage(double a, double b) {
    const double denom = a + b;
    if (denom == 0.0) {
        return 0.0;
    }
    return 2.0 * a * b / denom;
}

void enforceNeumannBoundaries(const Grid2D& grid, std::vector<double>& field) {
    if (grid.boundaryCondition != Grid2D::BoundaryKind::Neumann) {
        return;
    }
    if (grid.nx < 2 || grid.ny < 2) {
        return;
    }

    for (std::size_t j = 0; j < grid.ny; ++j) {
        const std::size_t leftInterior = grid.idx(1, j);
        const std::size_t rightInterior = grid.idx(grid.nx - 2, j);
        field[grid.idx(0, j)] = field[leftInterior];
        field[grid.idx(grid.nx - 1, j)] = field[rightInterior];
    }
    for (std::size_t i = 0; i < grid.nx; ++i) {
        const std::size_t bottomInterior = grid.idx(i, 1);
        const std::size_t topInterior = grid.idx(i, grid.ny - 2);
        field[grid.idx(i, 0)] = field[bottomInterior];
        field[grid.idx(i, grid.ny - 1)] = field[topInterior];
    }
}

void zeroDirichletBoundary(const Grid2D& grid, std::vector<double>& field) {
    if (grid.boundaryCondition != Grid2D::BoundaryKind::Dirichlet) {
        return;
    }
    if (grid.nx == 0 || grid.ny == 0) {
        return;
    }
    for (std::size_t i = 0; i < grid.nx; ++i) {
        field[grid.idx(i, 0)] = 0.0;
        field[grid.idx(i, grid.ny - 1)] = 0.0;
    }
    for (std::size_t j = 0; j < grid.ny; ++j) {
        field[grid.idx(0, j)] = 0.0;
        field[grid.idx(grid.nx - 1, j)] = 0.0;
    }
}

void removeMeanIfNeumann(const Grid2D& grid, std::vector<double>& field) {
    if (grid.boundaryCondition != Grid2D::BoundaryKind::Neumann) {
        return;
    }
    const std::size_t n = grid.nx * grid.ny;
    if (n == 0) {
        return;
    }
    double sum = 0.0;
    for (std::size_t idx = 0; idx < n; ++idx) {
        sum += field[idx];
    }
    const double mean = sum / static_cast<double>(n);
    for (std::size_t idx = 0; idx < n; ++idx) {
        field[idx] -= mean;
    }
}

void enforceNeumannGauge(const Grid2D& grid, std::vector<double>& field) {
    if (grid.boundaryCondition != Grid2D::BoundaryKind::Neumann) {
        return;
    }
    if (grid.nx == 0 || grid.ny == 0) {
        return;
    }
    field[grid.idx(0, 0)] = 0.0;
}

void computeJacobiDiagonal(const Grid2D& grid, std::vector<double>& diag) {
    const std::size_t n = grid.nx * grid.ny;
    diag.assign(n, 1.0);

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
            const double diagValue = (nuE + nuW) * invDx2 + (nuN + nuS) * invDy2;
            diag[p] = diagValue <= 0.0 ? 1.0 : diagValue;
        }
    }
}

bool emitProgress(ProgressSink* sink,
                  const SolveOptions& options,
                  std::size_t iter,
                  double relResidual,
                  const Clock::time_point& start,
                  Clock::time_point& lastEmission,
                  bool force) {
    if (sink == nullptr) {
        return true;
    }

    const auto now = Clock::now();
    const double elapsed = std::chrono::duration<double>(now - start).count();
    const bool intervalSatisfied = force || options.progressEverySec <= 0.0 ||
                                   std::chrono::duration<double>(now - lastEmission).count() >= options.progressEverySec;
    if (!intervalSatisfied) {
        return true;
    }

    ProgressSample sample{};
    sample.iter = iter;
    sample.relResidual = relResidual;
    sample.elapsedSeconds = elapsed;
    lastEmission = now;
    return sink->onProgress(sample);
}

void maybeWriteSnapshot(ProgressSink* sink,
                        const SolveOptions& options,
                        const Grid2D& grid,
                        std::size_t iter,
                        double relResidual,
                        const Clock::time_point& start) {
    if (sink == nullptr || options.snapshotEveryIters == 0 || iter == 0 ||
        (iter % options.snapshotEveryIters) != 0) {
        return;
    }

    const auto now = Clock::now();
    ProgressSample sample{};
    sample.iter = iter;
    sample.relResidual = relResidual;
    sample.elapsedSeconds = std::chrono::duration<double>(now - start).count();

    const std::optional<std::string> requestedPath = sink->requestFieldDump(sample);
    if (!requestedPath || requestedPath->empty()) {
        return;
    }

    std::filesystem::path path(*requestedPath);
    std::error_code ec;
    std::filesystem::create_directories(path.parent_path(), ec);

    const std::size_t strideX = std::max<std::size_t>(std::size_t(1), grid.nx / 64);
    const std::size_t strideY = std::max<std::size_t>(std::size_t(1), grid.ny / 64);

    std::vector<double> xs;
    std::vector<double> ys;
    std::vector<double> azSamples;
    xs.reserve((grid.nx / strideX + 1) * (grid.ny / strideY + 1));
    ys.reserve(xs.capacity());
    azSamples.reserve(xs.capacity());

    for (std::size_t j = 0; j < grid.ny; j += strideY) {
        for (std::size_t i = 0; i < grid.nx; i += strideX) {
            const std::size_t idx = grid.idx(i, j);
            xs.push_back(static_cast<double>(i) * grid.dx);
            ys.push_back(static_cast<double>(j) * grid.dy);
            azSamples.push_back(grid.Az[idx]);
        }
    }

    try {
        const std::vector<FieldColumnView> columns{{FieldColumnView{"Az", &azSamples}}};
        write_csv_field_map(path.string(), xs, ys, columns);
    } catch (const std::exception& ex) {
        std::cerr << "Warning: failed to write snapshot '" << path << "': " << ex.what() << "\n";
    }
}

SolveResult solveSORInternal(Grid2D& grid,
                             const SolveOptions& options,
                             ProgressSink* progressSink) {
    SolveResult report{};
    if (grid.nx < 3 || grid.ny < 3) {
        return report;
    }

    const double denom = std::sqrt(linops::squaredNorm(grid.Jz) + 1e-30);
    const double invDx2 = 1.0 / (grid.dx * grid.dx);
    const double invDy2 = 1.0 / (grid.dy * grid.dy);

    const auto start = Clock::now();
    auto lastProgress = start;

    for (std::size_t iter = 0; iter < options.maxIters; ++iter) {
        double residualNorm2 = 0.0;

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

                const double diag = (nuE + nuW) * invDx2 + (nuN + nuS) * invDy2;
                if (diag == 0.0) {
                    continue;
                }

                const double numerator = (nuE * grid.Az[e] + nuW * grid.Az[w]) * invDx2 +
                                         (nuN * grid.Az[nIdx] + nuS * grid.Az[s]) * invDy2 + grid.Jz[p];
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
                        grid.Az[nIdx] = updated;
                    }
                }

                const double flux = (nuE * (grid.Az[e] - updated) - nuW * (updated - grid.Az[w])) * invDx2 +
                                    (nuN * (grid.Az[nIdx] - updated) - nuS * (updated - grid.Az[s])) * invDy2;
                const double cellResidual = flux + grid.Jz[p];
                residualNorm2 += cellResidual * cellResidual;
            }
        }

        const double relResidual = std::sqrt(residualNorm2) / denom;
        report.relResidual = relResidual;
        report.iters = iter + 1;

        if (options.verbose) {
            std::cout << "Iteration " << report.iters << ": relResidual=" << relResidual << "\n";
        }

        const bool force = options.snapshotEveryIters > 0 && (report.iters % options.snapshotEveryIters) == 0;
        if (!emitProgress(progressSink, options, report.iters, relResidual, start, lastProgress, force)) {
            return report;
        }
        maybeWriteSnapshot(progressSink, options, grid, report.iters, relResidual, start);

        if (relResidual < options.tol) {
            report.converged = true;
            break;
        }
    }

    return report;
}

SolveResult solveCGInternal(Grid2D& grid,
                            const SolveOptions& options,
                            ProgressSink* progressSink) {
    SolveResult report{};
    const std::size_t n = grid.nx * grid.ny;
    if (grid.nx < 3 || grid.ny < 3 || n == 0) {
        return report;
    }

    enforceNeumannBoundaries(grid, grid.Az);
    removeMeanIfNeumann(grid, grid.Az);
    enforceNeumannGauge(grid, grid.Az);

    std::vector<double> r(n, 0.0);
    std::vector<double> z(n, 0.0);
    std::vector<double> p(n, 0.0);
    std::vector<double> Ap(n, 0.0);
    std::vector<double> diag(n, 1.0);

    computeJacobiDiagonal(grid, diag);

    std::vector<double> residualScratch;
    const double rhsNorm = std::sqrt(linops::squaredNorm(grid.Jz) + 1e-30);

    auto rebuildResidual = [&]() {
        residualScratch.resize(n);
        linops::computeResidual(grid, grid.Az.data(), residualScratch.data());
        removeMeanIfNeumann(grid, residualScratch);
        enforceNeumannGauge(grid, residualScratch);
        for (std::size_t idx = 0; idx < n; ++idx) {
            r[idx] = -residualScratch[idx];
        }
        zeroDirichletBoundary(grid, r);
        removeMeanIfNeumann(grid, r);
        enforceNeumannBoundaries(grid, r);
        enforceNeumannGauge(grid, r);
        return linops::norm(residualScratch) / rhsNorm;
    };

    double relResidual = rebuildResidual();
    if (options.warmStart && relResidual > 1.0) {
        std::fill(grid.Az.begin(), grid.Az.end(), 0.0);
        enforceNeumannBoundaries(grid, grid.Az);
        removeMeanIfNeumann(grid, grid.Az);
        enforceNeumannGauge(grid, grid.Az);
        relResidual = rebuildResidual();
    }
    double bestRelResidual = relResidual;
    std::size_t stagnationCounter = 0;
    constexpr std::size_t stagnationWindow = 200;
    bool warnedStagnation = false;

    const auto start = Clock::now();
    auto lastProgress = start;

    report.relResidual = relResidual;
    if (relResidual < options.tol) {
        report.converged = true;
        report.iters = 0;
        return report;
    }

    if (!emitProgress(progressSink, options, 0, relResidual, start, lastProgress, true)) {
        return report;
    }

    double rzPrev = linops::dot(r, z);  // initialise to zero via initial z
    bool firstIteration = true;

    for (std::size_t iter = 0; iter < options.maxIters; ++iter) {
        for (std::size_t idx = 0; idx < n; ++idx) {
            const double d = diag[idx];
            z[idx] = (d > 0.0) ? r[idx] / d : r[idx];
        }
        zeroDirichletBoundary(grid, z);
        removeMeanIfNeumann(grid, z);
        enforceNeumannBoundaries(grid, z);
        enforceNeumannGauge(grid, z);

        const double rz = linops::dot(r, z);
        if (!std::isfinite(rz)) {
            std::cerr << "CG encountered non-finite rz" << std::endl;
            break;
        }

        if (firstIteration) {
            p = z;
            firstIteration = false;
        } else {
            const double beta = (rzPrev != 0.0) ? (rz / rzPrev) : 0.0;
            for (std::size_t idx = 0; idx < n; ++idx) {
                p[idx] = z[idx] + beta * p[idx];
            }
        }
        zeroDirichletBoundary(grid, p);
        removeMeanIfNeumann(grid, p);
        enforceNeumannBoundaries(grid, p);
        enforceNeumannGauge(grid, p);

        linops::applyA(grid, p.data(), Ap.data());
        removeMeanIfNeumann(grid, Ap);
        enforceNeumannBoundaries(grid, Ap);
        enforceNeumannGauge(grid, Ap);
        double denom = linops::dot(p, Ap);
        if (denom <= 0.0 || !std::isfinite(denom)) {
            std::cerr << "CG denominator non-positive at iter " << iter << " (denom=" << denom
                      << ", rz=" << rz << ")" << std::endl;
            break;
        }

        const double alpha = rz / denom;
        for (std::size_t idx = 0; idx < n; ++idx) {
            grid.Az[idx] += alpha * p[idx];
            r[idx] -= alpha * Ap[idx];
        }
        zeroDirichletBoundary(grid, grid.Az);
        enforceNeumannBoundaries(grid, grid.Az);
        removeMeanIfNeumann(grid, grid.Az);
        removeMeanIfNeumann(grid, r);
        enforceNeumannBoundaries(grid, r);
        enforceNeumannGauge(grid, grid.Az);
        enforceNeumannGauge(grid, r);

        rzPrev = rz;
        relResidual = linops::norm(r) / rhsNorm;
        report.iters = iter + 1;
        report.relResidual = relResidual;

        if (options.verbose) {
            std::cout << "CG iteration " << report.iters << ": relResidual=" << relResidual << "\n";
        }

        const bool force = options.snapshotEveryIters > 0 && (report.iters % options.snapshotEveryIters) == 0;
        if (!emitProgress(progressSink, options, report.iters, relResidual, start, lastProgress, force)) {
            return report;
        }
        maybeWriteSnapshot(progressSink, options, grid, report.iters, relResidual, start);

        if (relResidual < options.tol) {
            report.converged = true;
            break;
        }

        if (relResidual + 1e-18 < bestRelResidual) {
            bestRelResidual = relResidual;
            stagnationCounter = 0;
        } else {
            ++stagnationCounter;
            if (stagnationCounter >= stagnationWindow && !warnedStagnation) {
                std::cerr << "Warning: CG residual stagnation detected; consider using --warm-start or --use-prolongation,"
                          << " or rerun with --solver sor for diagnostics." << std::endl;
                warnedStagnation = true;
            }
        }
    }

    return report;
}

Grid2D buildCoarseGrid(const Grid2D& fine, const SolveOptions& options) {
    Grid2D coarse;
    if (fine.nx < 3 || fine.ny < 3) {
        return coarse;
    }

    const std::size_t targetNx = options.coarseNx > 0 ? options.coarseNx : std::max<std::size_t>(3, (fine.nx + 1) / 2);
    const std::size_t targetNy = options.coarseNy > 0 ? options.coarseNy : std::max<std::size_t>(3, (fine.ny + 1) / 2);

    coarse.resize(targetNx, targetNy, 1.0, 1.0);
    coarse.boundaryCondition = fine.boundaryCondition;

    const double fineWidth = (fine.nx > 1) ? (fine.dx * static_cast<double>(fine.nx - 1)) : fine.dx;
    const double fineHeight = (fine.ny > 1) ? (fine.dy * static_cast<double>(fine.ny - 1)) : fine.dy;
    coarse.dx = (coarse.nx > 1) ? fineWidth / static_cast<double>(coarse.nx - 1) : fineWidth;
    coarse.dy = (coarse.ny > 1) ? fineHeight / static_cast<double>(coarse.ny - 1) : fineHeight;

    const auto sampleField = [&](const std::vector<double>& field, double u, double v) {
        const double x = std::clamp(u, 0.0, static_cast<double>(fine.nx - 1));
        const double y = std::clamp(v, 0.0, static_cast<double>(fine.ny - 1));
        const std::size_t i0 = static_cast<std::size_t>(std::floor(x));
        const std::size_t j0 = static_cast<std::size_t>(std::floor(y));
        const std::size_t i1 = std::min(fine.nx - 1, i0 + 1);
        const std::size_t j1 = std::min(fine.ny - 1, j0 + 1);
        const double tx = x - static_cast<double>(i0);
        const double ty = y - static_cast<double>(j0);

        const double v00 = field[fine.idx(i0, j0)];
        const double v10 = field[fine.idx(i1, j0)];
        const double v01 = field[fine.idx(i0, j1)];
        const double v11 = field[fine.idx(i1, j1)];

        const double v0 = v00 * (1.0 - tx) + v10 * tx;
        const double v1 = v01 * (1.0 - tx) + v11 * tx;
        return v0 * (1.0 - ty) + v1 * ty;
    };

    for (std::size_t j = 0; j < coarse.ny; ++j) {
        for (std::size_t i = 0; i < coarse.nx; ++i) {
            const double u = (coarse.nx > 1)
                                 ? static_cast<double>(i) * static_cast<double>(fine.nx - 1) /
                                       static_cast<double>(coarse.nx - 1)
                                 : 0.0;
            const double v = (coarse.ny > 1)
                                 ? static_cast<double>(j) * static_cast<double>(fine.ny - 1) /
                                       static_cast<double>(coarse.ny - 1)
                                 : 0.0;
            const std::size_t idx = coarse.idx(i, j);
            coarse.Jz[idx] = sampleField(fine.Jz, u, v);
            coarse.invMu[idx] = sampleField(fine.invMu, u, v);
        }
    }

    return coarse;
}

}  // namespace

SolveResult solveAz(Grid2D& grid,
                    const SolveOptions& options,
                    const InitialGuess* initialGuess,
                    ProgressSink* progress) {
    SolveOptions adjusted = options;

    if (!adjusted.warmStart && (!initialGuess || initialGuess->Az0 == nullptr) && !adjusted.useProlongation) {
        std::fill(grid.Az.begin(), grid.Az.end(), 0.0);
    }

    if (initialGuess != nullptr && initialGuess->Az0 != nullptr && initialGuess->Az0->size() == grid.Az.size()) {
        grid.Az = *initialGuess->Az0;
        enforceNeumannBoundaries(grid, grid.Az);
        removeMeanIfNeumann(grid, grid.Az);
        enforceNeumannGauge(grid, grid.Az);
        if (adjusted.kind == SolverKind::CG && adjusted.warmStart) {
            SolveOptions smoothOptions = adjusted;
            smoothOptions.kind = SolverKind::SOR;
            smoothOptions.maxIters = std::min<std::size_t>(smoothOptions.maxIters, std::size_t(50));
            smoothOptions.tol = std::max(smoothOptions.tol, 1e-3);
            smoothOptions.progressEverySec = 0.0;
            smoothOptions.snapshotEveryIters = 0;
            smoothOptions.verbose = false;
            solveSORInternal(grid, smoothOptions, nullptr);
            enforceNeumannBoundaries(grid, grid.Az);
            removeMeanIfNeumann(grid, grid.Az);
            enforceNeumannGauge(grid, grid.Az);
        }
    }

    Grid2D coarse;
    if (adjusted.useProlongation) {
        coarse = buildCoarseGrid(grid, adjusted);
        if (coarse.nx >= 3 && coarse.ny >= 3) {
            SolveOptions coarseOptions = adjusted;
            coarseOptions.kind = SolverKind::SOR;
            coarseOptions.useProlongation = false;
            coarseOptions.progressEverySec = 0.0;
            coarseOptions.snapshotEveryIters = 0;
            solveSORInternal(coarse, coarseOptions, nullptr);
            prolongateAzBilinear(coarse, grid);
            enforceNeumannBoundaries(grid, grid.Az);
            removeMeanIfNeumann(grid, grid.Az);
            enforceNeumannGauge(grid, grid.Az);
        }
    }

    if (adjusted.kind == SolverKind::CG && grid.boundaryCondition == Grid2D::BoundaryKind::Neumann) {
        if (options.verbose) {
            std::cout << "CG fallback to SOR for Neumann boundary condition" << std::endl;
        }
        SolveOptions sorOptions = adjusted;
        sorOptions.kind = SolverKind::SOR;
        return solveSORInternal(grid, sorOptions, progress);
    }

    switch (adjusted.kind) {
        case SolverKind::SOR:
            return solveSORInternal(grid, adjusted, progress);
        case SolverKind::CG:
            return solveCGInternal(grid, adjusted, progress);
    }
    return {};
}

SolveResult solveAz_GS_SOR(Grid2D& grid,
                           const SolveOptions& options,
                           const InitialGuess* initialGuess,
                           ProgressSink* progress) {
    SolveOptions sorOptions = options;
    sorOptions.kind = SolverKind::SOR;
    return solveAz(grid, sorOptions, initialGuess, progress);
}

SolveResult solveAz_CG(Grid2D& grid,
                       const SolveOptions& options,
                       const InitialGuess* initialGuess,
                       ProgressSink* progress) {
    SolveOptions cgOptions = options;
    cgOptions.kind = SolverKind::CG;
    return solveAz(grid, cgOptions, initialGuess, progress);
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

}  // namespace motorsim
