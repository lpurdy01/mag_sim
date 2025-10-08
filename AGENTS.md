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
- The CLI entry point `./build/motor_sim --scenario inputs/two_wire_cancel.json --solve --write-midline`
  loads the JSON, rasterises it onto the solver grid, runs the magnetostatic
  solve, and writes an optional midline CSV to `outputs/` for quick plotting.
- Reserved fields such as `"timeline"` are ignored by the current ingestor, so
  you can sketch future extensions without breaking compatibility.
