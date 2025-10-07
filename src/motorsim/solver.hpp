#pragma once

#include "grid.hpp"

namespace motorsim {

struct SolveOptions {
    std::size_t max_iters{1000};
    double tol{1e-6};
};

inline void solve_placeholder(Grid2D& /*Az*/, const SolveOptions& /*opts*/) {
    // no-op for bootstrap; replace with Gauss-Seidel/SOR later
}

}  // namespace motorsim
