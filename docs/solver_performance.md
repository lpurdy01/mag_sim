# Solver Performance & Complexity Notes

This note captures rule-of-thumb complexity estimates for the Gauss–Seidel/SOR
solver and records timing measurements gathered with `tools/solver_benchmark`.
The figures are intended to help contributors choose realistic grid sizes when
working inside lightweight environments such as GitHub Codespaces.

## 1. Complexity heuristics

* Each SOR sweep visits every interior cell once, so the cost per iteration is
  \(\mathcal{O}(n_x n_y)\).
* The number of iterations required to hit a fixed tolerance grows roughly with
  the grid dimension because Gauss–Seidel damps high-frequency error quickly but
  converges slowly on low-frequency modes. For uniform grids, iteration counts on
  the wire test scale approximately with \(\mathcal{O}(n_x)\).
* Wall-clock time is therefore well-approximated by
  \[
  t \approx \alpha \; n_x n_y N_\text{iter},
  \]
  where \(\alpha\) is the time per cell-iteration. On the reference hardware
  below \(\alpha \approx 2.2 \times 10^{-8}\) s (≈ 22 ms per million
  cell-iterations).

When planning new tests, estimate the runtime via the reported
"ms per million cell-iterations" (Section 2). For example, a 321×321 grid that
requires 4500 iterations will involve ≈464 million cell-iterations → about
10.2 seconds under the Codespaces baseline.

## 2. Reference timings (GitHub Codespaces 2 vCPU)

Hardware/software snapshot:

* GitHub Codespaces basic instance (2 vCPU / 4 GB RAM).
* Ubuntu 22.04 container with GCC 13.3, CMake 3.22.
* Solver compiled in `Release` mode.

| Grid (nx×ny) | Omega | Iterations | Avg time [ms] | ms / 1e6 cell-its | Throughput [Melements/s] | Status |
| ------------ | ----- | ---------- | ------------- | ----------------- | ------------------------ | ------ |
| 129×129      | 1.7   | 3760       | 1362          | 21.8              | 45.9                     | ✓ (relResidual ≈ 1.0×10⁻⁶) |
| 201×201      | 1.9   | 2872       | 2539          | 21.9              | 45.7                     | ✓ (relResidual ≈ 1.0×10⁻⁶) |
| 257×257      | 1.9   | 4654       | 6673          | 21.7              | 46.1                     | ✓ (relResidual ≈ 1.0×10⁻⁶) |

The solver sustains ~46 million cell-iterations per second on this machine,
which stays remarkably constant as grids scale up. The dominant lever on runtime
is the iteration count, which is reduced by using more aggressive relaxation
factors (up to \(\omega \approx 1.9\) before instability).

## 3. Benchmark workflow

1. Build the tool: `cmake --build build --target solver_benchmark -j`.
2. Run with custom parameters. Example:

   ```bash
   ./build/solver_benchmark --nx 201 --ny 201 --omega 1.9 --tol 1e-6 --repeats 5
   ```

3. Optional: append to a CSV log for historical comparison using
   `--csv outputs/solver_benchmarks.csv`.

The binary resets the solution between repeats so the reported timings are not
influenced by warm caches from previous runs. Include the hardware description
whenever posting new numbers to this table.

## 4. Recommendations for automated tests

* Keep CI-oriented tests at ≤ 201² cells with \(\omega\) tuned for rapid
  convergence.
* Increase tolerances (e.g., `tol = 5e-6`) for exploratory notebooks to trade
  accuracy for turnaround.
* Use the CSV logging mode when experimenting with algorithmic tweaks so
  performance regressions are easy to detect.

