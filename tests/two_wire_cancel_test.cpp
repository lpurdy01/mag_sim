#include "motorsim/ingest.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <numeric>
#include <utility>

int main() {
    using namespace motorsim;

    namespace fs = std::filesystem;
    const fs::path scenarioPath =
        (fs::path(__FILE__).parent_path() / "../inputs/tests/two_wire_cancel_test.json").lexically_normal();
    ScenarioSpec spec = loadScenarioFromJson(scenarioPath.string());

    if (spec.wires.size() != 2) {
        std::cerr << "two_wire_cancel_test: expected two wires in scenario\n";
        return 1;
    }

    constexpr double pi = 3.14159265358979323846;
    const auto analyticField = [&](double x, double y) {
        double bx = 0.0;
        double by = 0.0;
        for (const auto& wire : spec.wires) {
            const double dx = x - wire.x;
            const double dy = y - wire.y;
            const double r2 = dx * dx + dy * dy;
            if (r2 <= 0.0) {
                return std::pair<double, double>{std::numeric_limits<double>::quiet_NaN(),
                                                 std::numeric_limits<double>::quiet_NaN()};
            }
            const double coeff = motorsim::MU0 * wire.current / (2.0 * pi * r2);
            bx += coeff * (-dy);
            by += coeff * dx;
        }
        return std::pair<double, double>{bx, by};
    };

    const auto nearWire = [&](double x, double y) {
        for (const auto& wire : spec.wires) {
            const double r = std::hypot(x - wire.x, y - wire.y);
            if (r <= 1.2 * wire.radius) {
                return true;
            }
        }
        return false;
    };

    const std::size_t iMid = spec.nx / 2;
    const std::size_t jMid = spec.ny / 2;
    const double xMid = spec.originX + static_cast<double>(iMid) * spec.dx;
    const double cellArea = spec.dx * spec.dy;

    const auto evaluateField = [&](const Grid2D& grid, double& verticalRelErr, double& horizontalRelErr,
                                   double& avgBxMid) -> bool {
        double verticalRelErrSum = 0.0;
        double verticalBxAbsSum = 0.0;
        std::size_t verticalSamples = 0;

        for (std::size_t j = 1; j + 1 < spec.ny; ++j) {
            const double y = spec.originY + static_cast<double>(j) * spec.dy;
            if (std::abs(y) > 0.05) {
                continue;
            }
            if (nearWire(xMid, y)) {
                continue;
            }
            const auto [bxAnalytic, byAnalytic] = analyticField(xMid, y);
            if (!std::isfinite(bxAnalytic) || !std::isfinite(byAnalytic)) {
                continue;
            }
            const std::size_t idx = grid.idx(iMid, j);
            const double bxNum = grid.Bx[idx];
            const double byNum = grid.By[idx];
            const double analyticMag = std::hypot(bxAnalytic, byAnalytic);
            if (analyticMag < 1e-9) {
                continue;
            }
            const double diffMag = std::hypot(bxNum - bxAnalytic, byNum - byAnalytic);
            verticalRelErrSum += diffMag / analyticMag;
            verticalBxAbsSum += std::abs(bxNum);
            ++verticalSamples;
        }

        if (verticalSamples < spec.ny / 4) {
            std::cerr << "two_wire_cancel_test: insufficient vertical samples (" << verticalSamples << ")\n";
            return false;
        }

        verticalRelErr = verticalRelErrSum / static_cast<double>(verticalSamples);
        avgBxMid = verticalBxAbsSum / static_cast<double>(verticalSamples);

        double horizontalRelErrSum = 0.0;
        std::size_t horizontalSamples = 0;
        const double yRow = spec.originY + static_cast<double>(jMid) * spec.dy;
        for (std::size_t i = 1; i + 1 < spec.nx; ++i) {
            const double x = spec.originX + static_cast<double>(i) * spec.dx;
            if (std::abs(x) > 0.05) {
                continue;
            }
            if (nearWire(x, yRow)) {
                continue;
            }
            const auto [bxAnalytic, byAnalytic] = analyticField(x, yRow);
            if (!std::isfinite(bxAnalytic) || !std::isfinite(byAnalytic)) {
                continue;
            }
            const std::size_t idx = grid.idx(i, jMid);
            const double bxNum = grid.Bx[idx];
            const double byNum = grid.By[idx];
            const double analyticMag = std::hypot(bxAnalytic, byAnalytic);
            if (analyticMag < 1e-9) {
                continue;
            }
            const double diffMag = std::hypot(bxNum - bxAnalytic, byNum - byAnalytic);
            horizontalRelErrSum += diffMag / analyticMag;
            ++horizontalSamples;
        }

        if (horizontalSamples < spec.nx / 4) {
            std::cerr << "two_wire_cancel_test: insufficient horizontal samples (" << horizontalSamples << ")\n";
            return false;
        }

        horizontalRelErr = horizontalRelErrSum / static_cast<double>(horizontalSamples);

        for (const auto& wire : spec.wires) {
            double integrated = 0.0;
            for (std::size_t j = 0; j < spec.ny; ++j) {
                const double y = spec.originY + static_cast<double>(j) * spec.dy;
                for (std::size_t i = 0; i < spec.nx; ++i) {
                    const double x = spec.originX + static_cast<double>(i) * spec.dx;
                    const double r = std::hypot(x - wire.x, y - wire.y);
                    if (r <= wire.radius) {
                        integrated += grid.Jz[grid.idx(i, j)] * cellArea;
                    }
                }
            }
            const double tolerance = 0.15 * std::abs(wire.current) + 1e-9;
            if (std::abs(integrated - wire.current) > tolerance) {
                std::cerr << "two_wire_cancel_test: wire current mismatch. expected " << wire.current
                          << ", got " << integrated << "\n";
                return false;
            }
        }

        return true;
    };

    SolveOptions options{};
    options.maxIters = 5000;
    options.tol = 1e-6;
    options.omega = 1.7;

    double baselineResidual = 0.0;
    double baselineVertical = 0.0;
    double baselineHorizontal = 0.0;
    double baselineBxMid = 0.0;
    bool baselineRecorded = false;

    const std::array<std::pair<SolverKind, const char*>, 2> solverCases = {
        std::make_pair(SolverKind::SOR, "SOR"), std::make_pair(SolverKind::CG, "CG")};

    for (const auto& [solverKind, label] : solverCases) {
        options.kind = solverKind;
        Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
        rasterizeScenarioToGrid(spec, grid);

        const SolveResult report = solveAz(grid, options);
        if (!report.converged) {
            std::cerr << label << " solver failed to converge\n";
            return 1;
        }
        if (report.relResidual > options.tol) {
            std::cerr << label << " residual " << report.relResidual
                      << " exceeds tolerance " << options.tol << "\n";
            return 1;
        }

        computeB(grid);

        double verticalErr = 0.0;
        double horizontalErr = 0.0;
        double avgBx = 0.0;
        if (!evaluateField(grid, verticalErr, horizontalErr, avgBx)) {
            return 1;
        }

        if (verticalErr >= 0.25) {
            std::cerr << label << " vertical relative error too high (" << verticalErr << ")\n";
            return 1;
        }
        if (horizontalErr >= 0.22) {
            std::cerr << label << " horizontal relative error too high (" << horizontalErr << ")\n";
            return 1;
        }
        if (avgBx >= 2.0e-5) {
            std::cerr << label << " avg |Bx| on midline not near zero (" << avgBx << ")\n";
            return 1;
        }

        std::cout << label << " iters=" << report.iters << " verticalErr=" << verticalErr
                  << " horizontalErr=" << horizontalErr << " avg|Bx|=" << avgBx << "\n";

        if (!baselineRecorded) {
            baselineResidual = report.relResidual;
            baselineVertical = verticalErr;
            baselineHorizontal = horizontalErr;
            baselineBxMid = avgBx;
            baselineRecorded = true;
        } else {
            const double residualBound = baselineResidual * 1.1 + 1e-12;
            if (report.relResidual > residualBound) {
                std::cerr << label << " residual " << report.relResidual
                          << " exceeds parity bound " << residualBound << "\n";
                return 1;
            }
            if (verticalErr > baselineVertical + 0.02) {
                std::cerr << label << " vertical error " << verticalErr
                          << " exceeds allowed slack relative to baseline " << baselineVertical
                          << "\n";
                return 1;
            }
            if (horizontalErr > baselineHorizontal + 0.02) {
                std::cerr << label << " horizontal error " << horizontalErr
                          << " exceeds allowed slack relative to baseline " << baselineHorizontal
                          << "\n";
                return 1;
            }
            if (avgBx > baselineBxMid + 5.0e-6) {
                std::cerr << label << " avg |Bx| " << avgBx
                          << " exceeds allowed slack relative to baseline " << baselineBxMid
                          << "\n";
                return 1;
            }
        }
    }

    return 0;
}
