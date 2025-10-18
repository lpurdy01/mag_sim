#include "motorsim/probes.hpp"

#include "motorsim/types.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace motorsim {
namespace {

struct BilinearSample {
    double bx{0.0};
    double by{0.0};
};

BilinearSample sample_bilinear(const Grid2D& grid, double originX, double originY, double x, double y) {
    if (grid.Bx.empty() || grid.By.empty()) {
        throw std::runtime_error("Grid magnetic field arrays are empty; call computeB() first");
    }

    const double fx = (x - originX) / grid.dx;
    const double fy = (y - originY) / grid.dy;
    if (!(fx >= 0.0 && fy >= 0.0 && fx <= static_cast<double>(grid.nx - 1) &&
          fy <= static_cast<double>(grid.ny - 1))) {
        throw std::runtime_error("Probe loop samples lie outside the simulation domain");
    }

    const std::size_t maxI = grid.nx - 2;
    const std::size_t maxJ = grid.ny - 2;
    const std::size_t i0 = static_cast<std::size_t>(std::min(std::floor(fx), static_cast<double>(maxI)));
    const std::size_t j0 = static_cast<std::size_t>(std::min(std::floor(fy), static_cast<double>(maxJ)));
    const std::size_t i1 = i0 + 1;
    const std::size_t j1 = j0 + 1;

    const double tx = std::clamp(fx - static_cast<double>(i0), 0.0, 1.0);
    const double ty = std::clamp(fy - static_cast<double>(j0), 0.0, 1.0);

    const std::size_t idx00 = grid.idx(i0, j0);
    const std::size_t idx10 = grid.idx(i1, j0);
    const std::size_t idx01 = grid.idx(i0, j1);
    const std::size_t idx11 = grid.idx(i1, j1);

    const auto interpolate = [&](const std::vector<double>& field) {
        const double v00 = field[idx00];
        const double v10 = field[idx10];
        const double v01 = field[idx01];
        const double v11 = field[idx11];
        const double v0 = (1.0 - tx) * v00 + tx * v10;
        const double v1 = (1.0 - tx) * v01 + tx * v11;
        return (1.0 - ty) * v0 + ty * v1;
    };

    BilinearSample sample{};
    sample.bx = interpolate(grid.Bx);
    sample.by = interpolate(grid.By);
    return sample;
}

double polygon_signed_area(const std::vector<double>& xs, const std::vector<double>& ys) {
    const std::size_t count = xs.size();
    double area = 0.0;
    for (std::size_t i = 0; i < count; ++i) {
        const std::size_t j = (i + 1) % count;
        area += xs[i] * ys[j] - xs[j] * ys[i];
    }
    return 0.5 * area;
}

bool point_in_polygon(double x, double y, const std::vector<double>& xs, const std::vector<double>& ys) {
    const std::size_t count = xs.size();
    bool inside = false;
    for (std::size_t i = 0, j = count - 1; i < count; j = i++) {
        const double xi = xs[i];
        const double yi = ys[i];
        const double xj = xs[j];
        const double yj = ys[j];
        const double denom = yj - yi;
        if (std::abs(denom) < 1e-15) {
            continue;
        }
        const bool intersect = ((yi > y) != (yj > y)) &&
                               (x < (xj - xi) * (y - yi) / denom + xi);
        if (intersect) {
            inside = !inside;
        }
    }
    return inside;
}

double sample_component(const Grid2D& grid, std::size_t idx, FluxComponent component) {
    switch (component) {
        case FluxComponent::Bx:
            return grid.Bx[idx];
        case FluxComponent::By:
            return grid.By[idx];
        case FluxComponent::Bmag:
        default:
            return std::hypot(grid.Bx[idx], grid.By[idx]);
    }
}

double average_cell_component(const Grid2D& grid, std::size_t i, std::size_t j, FluxComponent component) {
    const std::size_t i1 = i + 1;
    const std::size_t j1 = j + 1;
    const std::size_t idx00 = grid.idx(i, j);
    const std::size_t idx10 = grid.idx(i1, j);
    const std::size_t idx01 = grid.idx(i, j1);
    const std::size_t idx11 = grid.idx(i1, j1);

    const double v00 = sample_component(grid, idx00, component);
    const double v10 = sample_component(grid, idx10, component);
    const double v01 = sample_component(grid, idx01, component);
    const double v11 = sample_component(grid, idx11, component);
    return 0.25 * (v00 + v10 + v01 + v11);
}

int clamp_start_index(double minCoord, double origin, double spacing, std::size_t limit) {
    const double continuous = (minCoord - origin) / spacing - 0.5;
    const double raw = std::ceil(continuous);
    const int idx = static_cast<int>(std::max(0.0, raw));
    return std::min(idx, static_cast<int>(limit));
}

int clamp_end_index(double maxCoord, double origin, double spacing, std::size_t limit) {
    const double continuous = (maxCoord - origin) / spacing - 0.5;
    const double raw = std::floor(continuous);
    const int idx = static_cast<int>(std::max(0.0, raw));
    return std::min(idx, static_cast<int>(limit));
}

}  // namespace

StressTensorResult evaluate_maxwell_stress_probe(const Grid2D& grid, double originX, double originY,
                                                 double dx, double dy, const std::vector<double>& xs,
                                                 const std::vector<double>& ys) {
    if (xs.size() != ys.size()) {
        throw std::runtime_error("Probe loop coordinate arrays must have matching lengths");
    }
    if (xs.size() < 3) {
        throw std::runtime_error("Probe loop requires at least three vertices");
    }
    if (dx <= 0.0 || dy <= 0.0) {
        throw std::runtime_error("Grid spacing must be positive for probe evaluation");
    }

    const std::size_t vertexCount = xs.size();
    const double area = polygon_signed_area(xs, ys);
    if (std::abs(area) < 1e-18) {
        throw std::runtime_error("Probe loop has near-zero area");
    }
    const double orientationSign = area >= 0.0 ? 1.0 : -1.0;

    const double sampleSpacing = std::min(dx, dy);
    StressTensorResult accum{};

    for (std::size_t i = 0; i < vertexCount; ++i) {
        const std::size_t j = (i + 1) % vertexCount;
        const double x0 = xs[i];
        const double y0 = ys[i];
        const double x1 = xs[j];
        const double y1 = ys[j];
        const double segDx = x1 - x0;
        const double segDy = y1 - y0;
        const double segLen = std::hypot(segDx, segDy);
        if (segLen < 1e-12) {
            continue;
        }

        const double invLen = 1.0 / segLen;
        const double tx = segDx * invLen;
        const double ty = segDy * invLen;

        double nx = orientationSign >= 0.0 ? ty : -ty;
        double ny = orientationSign >= 0.0 ? -tx : tx;
        const double normalLen = std::hypot(nx, ny);
        if (normalLen < 1e-18) {
            continue;
        }
        nx /= normalLen;
        ny /= normalLen;

        const int steps = std::max(1, static_cast<int>(std::ceil(segLen / sampleSpacing)));
        const double ds = segLen / static_cast<double>(steps);

        for (int step = 0; step < steps; ++step) {
            const double sMid = (static_cast<double>(step) + 0.5) * ds;
            const double px = x0 + tx * sMid;
            const double py = y0 + ty * sMid;

            const BilinearSample sample = sample_bilinear(grid, originX, originY, px, py);
            const double bx = sample.bx;
            const double by = sample.by;
            const double b2 = bx * bx + by * by;

            const double stress_xx = (bx * bx - 0.5 * b2) / MU0;
            const double stress_xy = (bx * by) / MU0;
            const double stress_yx = stress_xy;
            const double stress_yy = (by * by - 0.5 * b2) / MU0;

            const double tractionX = stress_xx * nx + stress_xy * ny;
            const double tractionY = stress_yx * nx + stress_yy * ny;

            const double dFx = tractionX * ds;
            const double dFy = tractionY * ds;

            accum.forceX += dFx;
            accum.forceY += dFy;
            accum.torqueZ += px * dFy - py * dFx;
        }
    }

    return accum;
}

double integrate_polygon_flux_component(const Grid2D& grid, double originX, double originY, double dx, double dy,
                                        const std::vector<double>& xs, const std::vector<double>& ys, double minX,
                                        double maxX, double minY, double maxY, FluxComponent component) {
    if (grid.Bx.empty() || grid.By.empty()) {
        throw std::runtime_error("Grid magnetic field arrays are empty; call computeB() first");
    }
    if (xs.size() != ys.size() || xs.size() < 3) {
        throw std::runtime_error("Polygon flux region requires matching coordinate arrays with at least three vertices");
    }
    if (!(dx > 0.0) || !(dy > 0.0)) {
        throw std::runtime_error("Grid spacing must be positive for flux integration");
    }

    const std::size_t maxCellI = grid.nx >= 2 ? grid.nx - 2 : 0;
    const std::size_t maxCellJ = grid.ny >= 2 ? grid.ny - 2 : 0;
    const int iStart = clamp_start_index(minX, originX, dx, maxCellI);
    const int iEnd = clamp_end_index(maxX, originX, dx, maxCellI);
    const int jStart = clamp_start_index(minY, originY, dy, maxCellJ);
    const int jEnd = clamp_end_index(maxY, originY, dy, maxCellJ);
    if (iEnd < iStart || jEnd < jStart) {
        return 0.0;
    }

    double flux = 0.0;
    for (int j = jStart; j <= jEnd; ++j) {
        const double yCenter = originY + (static_cast<double>(j) + 0.5) * dy;
        for (int i = iStart; i <= iEnd; ++i) {
            const double xCenter = originX + (static_cast<double>(i) + 0.5) * dx;
            if (!point_in_polygon(xCenter, yCenter, xs, ys)) {
                continue;
            }
            const double value = average_cell_component(grid, static_cast<std::size_t>(i),
                                                        static_cast<std::size_t>(j), component);
            flux += value * dx * dy;
        }
    }

    return flux;
}

double integrate_rect_flux_component(const Grid2D& grid, double originX, double originY, double dx, double dy,
                                     double minX, double maxX, double minY, double maxY, FluxComponent component) {
    if (grid.Bx.empty() || grid.By.empty()) {
        throw std::runtime_error("Grid magnetic field arrays are empty; call computeB() first");
    }
    if (!(dx > 0.0) || !(dy > 0.0)) {
        throw std::runtime_error("Grid spacing must be positive for flux integration");
    }
    if (!(maxX > minX) || !(maxY > minY)) {
        throw std::runtime_error("Rectangle flux region must have positive width and height");
    }

    const std::size_t maxCellI = grid.nx >= 2 ? grid.nx - 2 : 0;
    const std::size_t maxCellJ = grid.ny >= 2 ? grid.ny - 2 : 0;
    const int iStart = clamp_start_index(minX, originX, dx, maxCellI);
    const int iEnd = clamp_end_index(maxX, originX, dx, maxCellI);
    const int jStart = clamp_start_index(minY, originY, dy, maxCellJ);
    const int jEnd = clamp_end_index(maxY, originY, dy, maxCellJ);
    if (iEnd < iStart || jEnd < jStart) {
        return 0.0;
    }

    double flux = 0.0;
    for (int j = jStart; j <= jEnd; ++j) {
        const double yCenter = originY + (static_cast<double>(j) + 0.5) * dy;
        if (yCenter < minY || yCenter > maxY) {
            continue;
        }
        for (int i = iStart; i <= iEnd; ++i) {
            const double xCenter = originX + (static_cast<double>(i) + 0.5) * dx;
            if (xCenter < minX || xCenter > maxX) {
                continue;
            }
            const double value = average_cell_component(grid, static_cast<std::size_t>(i),
                                                        static_cast<std::size_t>(j), component);
            flux += value * dx * dy;
        }
    }

    return flux;
}

double compute_magnetic_coenergy(const Grid2D& grid, double dx, double dy) {
    if (!(dx > 0.0) || !(dy > 0.0)) {
        throw std::runtime_error("Grid spacing must be positive for co-energy integration");
    }

    const std::size_t count = grid.nx * grid.ny;
    if (grid.Bx.size() != count || grid.By.size() != count) {
        throw std::runtime_error("Magnetic field components unavailable; call computeB() first");
    }
    if (grid.Hx.size() != count || grid.Hy.size() != count) {
        throw std::runtime_error("Magnetic field intensity unavailable; call computeH() before co-energy integration");
    }

    const double cellArea = dx * dy;
    double energy = 0.0;
    const bool hasMagnetization = grid.Mx.size() == count && grid.My.size() == count;
    for (std::size_t idx = 0; idx < count; ++idx) {
        const double hx = grid.Hx[idx];
        const double hy = grid.Hy[idx];
        double mx = 0.0;
        double my = 0.0;
        if (hasMagnetization) {
            mx = grid.Mx[idx];
            my = grid.My[idx];
        }
        energy += grid.Bx[idx] * (hx + mx) + grid.By[idx] * (hy + my);
    }
    return 0.5 * energy * cellArea;
}

std::vector<BackEmfSample> compute_back_emf_series(const std::vector<std::size_t>& frameIndices,
                                                   const std::vector<double>& times,
                                                   const std::vector<double>& fluxes) {
    if (frameIndices.size() != times.size() || times.size() != fluxes.size()) {
        throw std::runtime_error("Back-EMF series requires matching frame, time, and flux arrays");
    }
    if (frameIndices.size() < 2) {
        throw std::runtime_error("Back-EMF series requires at least two samples");
    }

    std::vector<BackEmfSample> samples;
    samples.reserve(frameIndices.size() - 1);
    for (std::size_t idx = 0; idx + 1 < frameIndices.size(); ++idx) {
        const double dt = times[idx + 1] - times[idx];
        if (std::abs(dt) < 1e-18) {
            throw std::runtime_error("Back-EMF frames must have distinct timestamps");
        }
        const double emf = -(fluxes[idx + 1] - fluxes[idx]) / dt;
        BackEmfSample sample{};
        sample.frameStart = frameIndices[idx];
        sample.frameEnd = frameIndices[idx + 1];
        sample.timeStart = times[idx];
        sample.timeEnd = times[idx + 1];
        sample.fluxStart = fluxes[idx];
        sample.fluxEnd = fluxes[idx + 1];
        sample.emf = emf;
        samples.push_back(sample);
    }

    return samples;
}

}  // namespace motorsim
