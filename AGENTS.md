# Agent Onboarding Notes

Start by reviewing the existing documentation before making changes:

- **Project overview:** `README.md`
- **Contribution rules and coding standards:** `CONTRIBUTING.md`
- **Physics and solver details:** `docs/math_and_solver.md`
- **Developer workflow tips (CLI, VS Code, Codespaces):** `docs/dev_environment.md`
- **Performance guidance and benchmarking workflow:** `docs/solver_performance.md`
- **Visualization helper:** `python/visualize_wire.py`

Check for additional nested `AGENTS.md` files in subdirectories when you work within them.

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
  wire overlays (defaults to the bundled two-wire cancellation scenario).
