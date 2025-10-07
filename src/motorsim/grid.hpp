#pragma once

#include <cstddef>
#include <vector>

namespace motorsim {

struct Grid2D {
    std::size_t nx{0}, ny{0};
    double dx{1.0}, dy{1.0};
    std::vector<double> data;

    Grid2D(std::size_t nx_, std::size_t ny_, double dx_ = 1.0, double dy_ = 1.0)
        : nx(nx_), ny(ny_), dx(dx_), dy(dy_), data(nx_ * ny_, 0.0) {}

    inline double &at(std::size_t i, std::size_t j) { return data[j * nx + i]; }
    inline const double &at(std::size_t i, std::size_t j) const { return data[j * nx + i]; }
};

}  // namespace motorsim
