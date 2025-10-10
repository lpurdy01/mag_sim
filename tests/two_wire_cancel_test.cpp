#include "motorsim/ingest.hpp"
#include "motorsim/solver.hpp"
#include "motorsim/types.hpp"

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

    Grid2D grid(spec.nx, spec.ny, spec.dx, spec.dy);
    rasterizeScenarioToGrid(spec, grid);

    SolveOptions options{};
    options.maxIters = 5000;
    options.tol = 1e-6;
    options.omega = 1.7;

    const SolveReport report = solveAz_GS_SOR(grid, options);
    if (!report.converged) {
        std::cerr << "two_wire_cancel_test: solver failed to converge\n";
        return 1;
    }
    if (report.relResidual > options.tol) {
        std::cerr << "two_wire_cancel_test: residual " << report.relResidual
                  << " exceeds tolerance " << options.tol << "\n";
        return 1;
    }

    computeB(grid);

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

    double verticalRelErrSum = 0.0;
    double verticalBxAbsSum = 0.0;
    std::size_t verticalSamples = 0;

    for (std::size_t j = 1; j + 1 < spec.ny; ++j) {
        const double y = spec.originY + static_cast<double>(j) * spec.dy;
        if (std::abs(y) > 0.05) {  // avoid boundary influence near y = ±0.1 m
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
        return 1;
    }

    const double avgVerticalRelErr = verticalRelErrSum / static_cast<double>(verticalSamples);
    if (!(avgVerticalRelErr < 0.25)) {
        std::cerr << "two_wire_cancel_test: vertical relative error too high (" << avgVerticalRelErr
                  << ")\n";
        return 1;
    }

    const double avgBxMid = verticalBxAbsSum / static_cast<double>(verticalSamples);
    if (!(avgBxMid < 2.0e-5)) {
        std::cerr << "two_wire_cancel_test: Bx on vertical midline not near zero (avg |Bx|=" << avgBxMid << ")\n";
        return 1;
    }

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
        return 1;
    }

    const double avgHorizontalRelErr = horizontalRelErrSum / static_cast<double>(horizontalSamples);
    if (!(avgHorizontalRelErr < 0.22)) {
        std::cerr << "two_wire_cancel_test: horizontal relative error too high (" << avgHorizontalRelErr
                  << ")\n";
        return 1;
    }

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
            return 1;
        }
    }

    std::cout << "two_wire_cancel_test: verticalRelErr=" << avgVerticalRelErr
              << ", horizontalRelErr=" << avgHorizontalRelErr << ", iters=" << report.iters << "\n";
    return 0;
}
