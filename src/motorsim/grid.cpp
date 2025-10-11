#include "motorsim/grid.hpp"

#include <algorithm>
#include <cmath>

namespace motorsim {

void prolongateAzBilinear(const Grid2D& coarse, Grid2D& fine) {
    if (coarse.Az.empty() || fine.Az.empty()) {
        return;
    }
    if (coarse.nx < 2 || coarse.ny < 2) {
        fine.Az.assign(fine.Az.size(), 0.0);
        return;
    }

    const auto sample = [&](double u, double v) {
        const double x = std::clamp(u, 0.0, static_cast<double>(coarse.nx - 1));
        const double y = std::clamp(v, 0.0, static_cast<double>(coarse.ny - 1));
        const std::size_t i0 = static_cast<std::size_t>(std::floor(x));
        const std::size_t j0 = static_cast<std::size_t>(std::floor(y));
        const std::size_t i1 = std::min(coarse.nx - 1, i0 + 1);
        const std::size_t j1 = std::min(coarse.ny - 1, j0 + 1);
        const double tx = x - static_cast<double>(i0);
        const double ty = y - static_cast<double>(j0);

        const double v00 = coarse.Az[coarse.idx(i0, j0)];
        const double v10 = coarse.Az[coarse.idx(i1, j0)];
        const double v01 = coarse.Az[coarse.idx(i0, j1)];
        const double v11 = coarse.Az[coarse.idx(i1, j1)];

        const double vx0 = v00 * (1.0 - tx) + v10 * tx;
        const double vx1 = v01 * (1.0 - tx) + v11 * tx;
        return vx0 * (1.0 - ty) + vx1 * ty;
    };

    for (std::size_t j = 0; j < fine.ny; ++j) {
        for (std::size_t i = 0; i < fine.nx; ++i) {
            const double u = (fine.nx > 1)
                                 ? static_cast<double>(i) * static_cast<double>(coarse.nx - 1) /
                                       static_cast<double>(fine.nx - 1)
                                 : 0.0;
            const double v = (fine.ny > 1)
                                 ? static_cast<double>(j) * static_cast<double>(coarse.ny - 1) /
                                       static_cast<double>(fine.ny - 1)
                                 : 0.0;
            fine.Az[fine.idx(i, j)] = sample(u, v);
        }
    }
}

}  // namespace motorsim
