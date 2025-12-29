# Flask GUI Agent Notes

- This package hosts the Flask-based prototype web UI for `mag_sim`.
- Keep browser dependencies minimal: prefer server-sent events (SSE) and
  lightweight Bootstrap styling already included in the templates. Avoid adding
  large frontend frameworks unless explicitly requested.
- Reuse the existing `SimulationManager` helper when wiring new routes. Ensure
  updates remain thread-safe and continue to enforce the single-simulation
  constraint.
- Stage 1 of the GUI overhaul introduced geometry previews, live field-map
  visualisation, the DXF workspace, and the `/visualization.png` customisation
  endpoint. The UI now includes a CAD-style top bar, workflow sidebar, and a
  project-aware workspace summary. Keep the File-menu actions
  (`/project/reset`, `/project/scenario`, `/project/export_dxf`) and
  per-session project state intact when iterating on features.
- The results panel now exposes playback controls powered by SSE frame events.
  Preserve the `sequence_index`, `frame_index`, and `frames` metadata so the
  animation gallery can stay in sync with the solver log, and extend the JS
  helpers (`registerFrame`, `advanceFrame`) when adding new visual streams.
- DXF handling lives in `python/gui/dxf_utils.py`; reuse its helpers when
  parsing uploads, rendering previews, or exporting geometry. Layer selections
  are stored in `PROJECTS[pid]['dxf_files']`—mutations should be mirrored in
  tests to preserve coverage.
- Sample DXFs for the induction workflow live under
  `docs/assets/dxf/induction_demo/`; keep them aligned with
  `inputs/induction_gui_demo.json` and update the workflow doc when you add new
  layers or assets.
- Follow-up work should expand tests when touching related code—fixtures clear
  `PROJECTS`, so new project-oriented behaviour needs explicit coverage.
- End-to-end tests live under `tests/e2e/` and rely on `data-testid` hooks in
  the templates. Preserve or extend those attributes when adjusting the markup
  so Playwright selectors remain stable.
- When modifying behaviour, update the user guide (`docs/user-guide/gui_flask.md`)
  and extend `tests/test_gui_flask.py` so new code paths stay covered.
- Runtime artefacts such as uploads and generated logs should stay out of the
  repository. Use or extend the ignore rules in this directory if new folders
  are introduced.
- The helper scripts `scripts/setup_gui_env.sh` and
  `scripts/maintain_gui_env.sh` keep local environments reproducible—run them
  before hacking on the GUI or when tidying accumulated artefacts.
