#pragma once

#include <cstddef>
#include <vector>

#include "motorsim/grid.hpp"
#include "motorsim/types.hpp"

namespace motorsim {

struct StressTensorResult {
    double forceX{0.0};
    double forceY{0.0};
    double torqueZ{0.0};
};

StressTensorResult evaluate_maxwell_stress_probe(const Grid2D& grid, double originX, double originY,
                                                 double dx, double dy, const std::vector<double>& xs,
                                                 const std::vector<double>& ys);

struct BackEmfSample {
    std::size_t frameStart{0};
    std::size_t frameEnd{0};
    double timeStart{0.0};
    double timeEnd{0.0};
    double fluxStart{0.0};
    double fluxEnd{0.0};
    double emf{0.0};
};

double integrate_polygon_flux_component(const Grid2D& grid, double originX, double originY, double dx, double dy,
                                        const std::vector<double>& xs, const std::vector<double>& ys, double minX,
                                        double maxX, double minY, double maxY, FluxComponent component);

double integrate_rect_flux_component(const Grid2D& grid, double originX, double originY, double dx, double dy,
                                     double minX, double maxX, double minY, double maxY, FluxComponent component);

double compute_magnetic_coenergy(const Grid2D& grid, double dx, double dy);

std::vector<BackEmfSample> compute_back_emf_series(const std::vector<std::size_t>& frameIndices,
                                                   const std::vector<double>& times,
                                                   const std::vector<double>& fluxes);

}  // namespace motorsim
