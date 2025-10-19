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

Always finish a feature branch by running `scripts/run_ci_checks.sh`. The helper
mirrors `.github/workflows/ci.yml`, solving the demo scenarios, exercising the
Python visualisers, and collecting the same artefacts that GitHub Actions
uploads. When the script passes locally you can be confident the hosted CI will
behave the same way.

`build/motor_sim` loads JSON scenarios that describe regions, currents, and
optional permanent magnets. The automated regression suite lives in the
`tests/` directory and spans the analytic \(\mu_0 I / (2\pi r)\) wire case, the
planar permeability interface comparison, a magnet strip scenario, torque and
back-EMF probes, rotor ripple timelines, a frequency-domain skin-depth check,
and a transient magnetic-diffusion fixture. Each scenario ships with a compact
JSON under `inputs/tests/` so CI-friendly runs stay quick while maintaining
coverage across the solver features.

Runtime flags worth highlighting:

* `--solver {sor|cg|harmonic}` chooses between Gauss–Seidel, PCG, and the
  frequency-domain formulation.
* `--pc {none|jacobi|ssor}` selects the CG preconditioner. Jacobi is a good
  default on lightly conductive grids; SSOR speeds up highly anisotropic cases.
* `--progress-history <csv>` together with `--progress-every 0` captures
  per-iteration residuals for later plotting or CI artefacts.

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
Codespaces resource usage and customization tips for other IDEs. The same guide
outlines the CI mirroring workflow in more detail.

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
`--streamlines`, and `--color-scale log` for CI-friendly renders. For geometry-
centric animations, use `python/generate_rotor_animation.py` to combine a
scenario, mechanical trace, and optional circuit currents into coloured slot
animations; the CI workflow renders GIF/PNG pairs for the synchronous, DC, and
induction demos.

## Sample scenarios

* `inputs/line_current_interface.json` — two-material validation that exercises
  the planar permeability interface analytic case described in
  `docs/math_and_solver.md`.
* `inputs/iron_ring_demo.json` — heterogeneous permeability demo with a polygon
  iron ring and six alternating-current wires in the bore. Solve it via
  `./build/motor_sim --scenario inputs/iron_ring_demo.json --solve --outputs all`
  and render a
  PNG with `python/visualize_scenario_field.py --scenario inputs/iron_ring_demo.json \
  --field-map outputs/iron_ring_field.csv --save outputs/iron_ring_field.png`.
  (CI runs capture the same render as an artifact, so the repository stays free
  of committed binaries.)

## Time-series demos

The synchronous motor walkthroughs live under `python/` and emit JSON scenarios
that `motor_sim` can solve frame-by-frame. Both generators expose compact
`ci` profiles for regression and larger `hires` presets for offline studies.

* **Rotating stator field** –
  `python/gen_three_phase_stator.py --profile ci --out inputs/three_phase_stator_ci.json`
  reproduces the original rotating-field showcase without a rotor. Run it with
  `./build/motor_sim --scenario inputs/three_phase_stator_ci.json --solve --parallel-frames \
  --vtk-series outputs/three_phase_ci.pvd --tol 5e-6 --max-iters 40000` and inspect
  the bore-angle CSV or ParaView `.pvd` series. `docs/three_phase_stator.md`
  collects further tips and animation commands.
* **PM motor spin-up** –
  `python/gen_three_phase_pm_motor.py --profile ci --mode spinup --out inputs/three_phase_pm_motor_spinup_ci.json`
  adds the permanent-magnet rotor, lumped RL circuits, and the RK4 mechanical
  integrator. Slot polygons carry explicit turn counts (60 per slot) and a 0.55
  copper fill fraction so the stator’s ampere-turn budget matches the magnet
  linkage. The generator carves the magnet out of the rotor iron and assigns it
  a near-air permeability (μᵣ≈1.05) so the trimmed 1×10⁵ A/m magnetisation
  produces bore flux on the same order as the stator coils. The CI profile keeps
  polygon tessellation lean (tens of vertices per circle) and rounds coordinates
  to four decimals so the stored scenario stays readable while still honouring
  slot symmetry. Solving the generated JSON with `--vtk-series` writes
  `pm_motor_spinup_frame_###.vti`, torque CSVs, and a
  `pm_motor_spinup_mechanical.csv` history that `python/check_pm_spinup.py`
  validates. See `docs/three_phase_pm_motor.md` for scenario parameters and
  workflow guidance.
* **DC motor spin-up** –
  `python/gen_dc_motor.py --profile ci --mode spinup --out inputs/dc_motor_spinup_ci.json`
  introduces a commutated armature whose coil links reference rotor angle driven
  segments. The generator balances stator field and armature strength by pairing
  matched ampere-turn budgets (220-turn field poles at 18 V / 1.6 Ω, 110-turn
  armature slots at 12 V / 0.75 Ω) and exposes copper fill percentages so the
  deposited current density stays realistic. During timeline solves the circuit
  layer flips coil orientation as `dc_rotor` sweeps through ±90° so the torque
  sign remains positive. Outputs mirror the other demos—VTK field frames,
  `dc_motor_torque.csv`, outlines, and a `dc_motor_mechanical.csv` trace that the
  shared spin-up checker can validate with `--rotor dc_rotor`. See
  `docs/dc_commutated_motor.md` for commutator schema details and tuning tips.
* **Induction motor spin-up** –
  `python/gen_three_phase_induction_motor.py --profile ci --mode spinup --out inputs/three_phase_induction_motor_spinup_ci.json`
  swaps the permanent magnet for a conductive-bar cage, enables the transient
  Crank–Nicolson solve, and lets the RK4 mechanical loop react to eddy-current
  torque. The helper emits VTK series, outline polydata, and a mechanical trace
  (`induction_motor_mechanical.csv`) that the shared
  `python/check_pm_spinup.py --rotor induction_rotor` validator can inspect. See
  `docs/three_phase_induction_motor.md` for the full walkthrough and CI fixture
  details.
