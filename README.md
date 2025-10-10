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

`build/motor_sim` loads JSON scenarios that describe regions, currents, and
optional permanent magnets. The automated regression suite lives in the
`tests/` directory and includes the analytic \(\mu_0 I / (2\pi r)\) wire case,
the planar permeability interface comparison, a magnet strip scenario, and a
two-wire cancellation sanity check.

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
* `CONTRIBUTING.md` — development workflow expectations, including the
  `scripts/run_ci_checks.sh` harness that mirrors GitHub Actions locally.

Python helpers for plotting the analytic wire validation live under
`python/`. Run `python/visualize_wire.py` after the regression test to compare
the numerical centreline magnitude against the analytic profile for a finite
radius wire (the drop to zero at the core is expected for a solid conductor).
For a top-down view or to inspect heterogeneous regions, use
`python/visualize_scenario_field.py`. The script accepts the scenario JSON and a
field-map CSV (`--field-map`) and provides switches such as `--draw-boundaries`,
`--streamlines`, and `--color-scale log` for CI-friendly renders.

## Sample scenarios

* `inputs/line_current_interface.json` — two-material validation that exercises
  the planar permeability interface analytic case described in
  `docs/math_and_solver.md`.
* `inputs/iron_ring_demo.json` — heterogeneous permeability demo with a polygon
  iron ring and six alternating-current wires in the bore. Solve it via
  `./build/motor_sim --scenario inputs/iron_ring_demo.json --solve` and render a
  PNG with `python/visualize_scenario_field.py --scenario inputs/iron_ring_demo.json \
  --field-map outputs/iron_ring_field.csv --save outputs/iron_ring_field.png`.
  (CI runs capture the same render as an artifact, so the repository stays free
  of committed binaries.)
