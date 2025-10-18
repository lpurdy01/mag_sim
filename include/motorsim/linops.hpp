#pragma once

#include "grid.hpp"

#include <cstddef>
#include <vector>

namespace motorsim::linops {

/**
 * @brief Apply the discrete magnetostatic operator A to a vector.
 *
 * Computes y = A x using the current grid configuration.
 */
void applyA(const Grid2D& grid, const double* x, double* y);

/**
 * @brief Compute the residual r = A x - J for the current grid state.
 */
void computeResidual(const Grid2D& grid, const double* x, double* residual);

/**
 * @brief Compute the Euclidean dot product of two vectors.
 */
double dot(std::size_t n, const double* a, const double* b);

double dot(const std::vector<double>& a, const std::vector<double>& b);

/**
 * @brief Compute the squared Euclidean norm of a vector.
 */
double squaredNorm(std::size_t n, const double* x);

double squaredNorm(const std::vector<double>& x);

/**
 * @brief Compute the Euclidean norm of a vector.
 */
double norm(std::size_t n, const double* x);

double norm(const std::vector<double>& x);

/**
 * @brief Compute the residual norm ||A x - J||_2.
 *
 * Scratch storage is reused across calls to avoid repeated allocations.
 */
double residualNorm(const Grid2D& grid, const double* x, std::vector<double>& scratch);

/**
 * @brief Perform the axpy operation y = y + alpha * x.
 */
void axpy(std::size_t n, double alpha, const double* x, double* y);

/**
 * @brief Scale a vector in-place: x = alpha * x.
 */
void scal(std::size_t n, double alpha, double* x);

/**
 * @brief Set all entries of a vector to zero.
 */
void zero(std::size_t n, double* x);

}  // namespace motorsim::linops
