// filename: grid.hpp
// part of 2D Electromagnetic Motor Simulator
// MIT License

#pragma once

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace motorsim {

/**
 * @brief Structured uniform grid container for 2D magnetostatic simulations.
 */
struct Grid2D {
    enum class BoundaryKind { Dirichlet, Neumann };

    std::size_t nx{0};
    std::size_t ny{0};
    double dx{1.0};
    double dy{1.0};
    BoundaryKind boundaryCondition{BoundaryKind::Dirichlet};
    std::vector<double> Az;
    std::vector<double> Jz;
    std::vector<double> invMu;
    std::vector<double> Mx;
    std::vector<double> My;
    std::vector<double> Bx;
    std::vector<double> By;
    std::vector<double> Hx;
    std::vector<double> Hy;

    Grid2D() = default;

    Grid2D(std::size_t nxIn, std::size_t nyIn, double dxIn, double dyIn)
        : nx(nxIn), ny(nyIn), dx(dxIn), dy(dyIn),
          Az(nxIn * nyIn, 0.0), Jz(nxIn * nyIn, 0.0), invMu(nxIn * nyIn, 0.0),
          Mx(nxIn * nyIn, 0.0), My(nxIn * nyIn, 0.0) {}

    [[nodiscard]] inline std::size_t idx(std::size_t i, std::size_t j) const {
        return j * nx + i;
    }

    [[nodiscard]] inline bool inBounds(std::size_t i, std::size_t j) const {
        return i < nx && j < ny;
    }

    [[nodiscard]] inline double& az(std::size_t i, std::size_t j) {
        if (!inBounds(i, j)) {
            throw std::out_of_range("Grid2D::az index out of range");
        }
        return Az[idx(i, j)];
    }

    [[nodiscard]] inline const double& az(std::size_t i, std::size_t j) const {
        if (!inBounds(i, j)) {
            throw std::out_of_range("Grid2D::az index out of range");
        }
        return Az[idx(i, j)];
    }

    void resize(std::size_t nxIn, std::size_t nyIn, double dxIn, double dyIn) {
        nx = nxIn;
        ny = nyIn;
        dx = dxIn;
        dy = dyIn;
        const std::size_t count = nx * ny;
        Az.assign(count, 0.0);
        Jz.assign(count, 0.0);
        invMu.assign(count, 0.0);
        Mx.assign(count, 0.0);
        My.assign(count, 0.0);
        Bx.clear();
        By.clear();
        Hx.clear();
        Hy.clear();
    }
};

}  // namespace motorsim
