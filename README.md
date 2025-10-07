# 2D Electromagnetic Motor Simulator

Minimal C++17 playground for magnetostatic experiments. The repo ships with a
uniform-grid solver, analytic validation test, and lightweight tooling that is
friendly to GitHub Codespaces and VS Code Remote workflows.

## Build & Test (CLI)

```bash
scripts/setup_env.sh           # optional helper for fresh environments
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
ctest --test-dir build --output-on-failure
```

`build/motor_sim` is a placeholder executable. The automated regression lives
in `analytic_wire_test` and mirrors the analytic \(\mu_0 I / (2\pi r)\) wire case.

## VS Code & Codespaces quickstart

The repository already contains `.vscode/` settings tuned for the CMake Tools
and C/C++ extensions.

* **Tasks** (`Terminal → Run Task`):
  * `cmake-configure` — configure the build directory in Debug mode.
  * `cmake-build` — build all targets, depends on configure.
  * `ctest` — execute the full test suite with failure output.
  * `run-analytic-wire-test` — compile and run the validation test binary.
  * `run-solver-benchmark` — compile and launch the benchmark helper with a
    lean 129×129 grid.
* **Debugging** (`Run and Debug` sidebar): launch configs exist for the main
  application, the analytic test, and the solver benchmark.

See `docs/dev_environment.md` for a complete walkthrough including advice on
Codespaces resource usage and customization tips for other IDEs.

## Solver benchmarking & performance guidance

Use `build/solver_benchmark` to profile convergence speed for different grid
resolutions and tolerances:

```bash
./build/solver_benchmark --nx 129 --ny 129 --max-iters 5000 --tol 1e-6
```

The tool prints throughput metrics and can append results to CSV for historical
comparison. Practical rules of thumb and sample measurements are documented in
`docs/solver_performance.md`.

## Documentation

* `docs/math_and_solver.md` — physics background, discretisation, and test plan.
* `docs/dev_environment.md` — IDE setup, VS Code workflow, and debugging tips.
* `docs/solver_performance.md` — complexity discussion and benchmark reference.

Python helpers for plotting the analytic wire validation live under
`python/`. Run `python/visualize_wire.py` after the regression test to compare
the numerical centreline magnitude against the analytic profile for a finite
radius wire (the drop to zero at the core is expected for a solid conductor).
For a top-down view and to inspect the azimuthal direction of **B**, use
`python/visualize_wire_field.py`, which overlays a quiver plot on the 2D field
map emitted by the test.
