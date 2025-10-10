# Agent Onboarding Notes

Start by reviewing the existing documentation before making changes:

- **Project overview:** `README.md`
- **Contribution rules and coding standards:** `CONTRIBUTING.md`
- **Physics and solver details:** `docs/math_and_solver.md`
- **Developer workflow tips (CLI, VS Code, Codespaces):** `docs/dev_environment.md`
- **Performance guidance and benchmarking workflow:** `docs/solver_performance.md`
- **Visualization helper:** `python/visualize_wire.py`

Check for additional nested `AGENTS.md` files in subdirectories when you work within them.

## Scenario authoring quick tip
- For heterogeneous permeability cases, build a two-halfspace scenario (left/right split)
  using the new `halfspace` region type introduced in spec v0.2.
- Complex cross-sections (e.g. iron rings with inner bores) can be described via
  layered `polygon` regions—later entries override earlier materials, so an
  outer iron polygon followed by an inner air polygon carves out a bore.

## Authoring a Scenario (v0.1)
- Use `python/scenario_api.py` to build scenarios in Python and emit JSON. The
  helper mirrors the C++ `ScenarioSpec` structure and writes to `inputs/`.
- Declare any desired exports via `FieldMapOutput`/`LineProbeOutput` so the
  solver knows which CSV artefacts to generate.
- Run `./build/motor_sim --scenario inputs/two_wire_cancel.json --solve` to
  rasterise, solve, and emit the requested outputs. `--list-outputs` shows the
  available IDs and `--outputs midline` restricts the run to a subset (use
  `none` to skip emission). The legacy `--write-midline` flag is still present
  for quick one-off dumps.
- Reserved fields such as `"timeline"` are ignored by the current ingestor, so
  you can sketch future extensions without breaking compatibility.
- Use `python/visualize_scenario_field.py` to plot the exported field maps with
  wire overlays (defaults to the bundled two-wire cancellation scenario). Pass
  `--color-scale log` for a logarithmic magnitude map when the dynamic range is
  large.

## Visualisation and CI artefacts
- Do not commit binary artefacts (PNGs, etc.) to the repository. The CI
  workflow renders figures at runtime and uploads them as GitHub Action
  artefacts via `actions/upload-artifact`. Expect the `solver-artifacts`
  bundle to contain field-map renders for the two-wire, interface, iron-ring,
  and magnet-strip scenarios, the analytic wire line plot, all exported CSVs,
  and a `test_accuracy_report.txt` summary captured from the regression
  executables.
- Before pushing CI-affecting changes, run `scripts/run_ci_checks.sh` to mirror
  the workflow locally; it exercises the same solves, visualisations, and
  artifact collection steps so argument regressions (e.g. flag typos) are
  caught early.
- The field visualiser exposes extra switches to aid debugging: use
  `--draw-boundaries` to overlay material and magnet outlines, `--streamlines`
  to trace the flow, and `--overlay-analytic-interface` to compare against the
  method-of-images solution for planar interface scenarios. Pair with
  `--hide-vectors` if the quiver clutter is distracting.
- The iron-ring demo requires a slightly looser nonlinear solve than the other
  regressions; the CI workflow invokes `motor_sim` with
  `--tol 5e-6 --max-iters 40000` for that scenario so the artefact rendering
  succeeds without manual intervention.
