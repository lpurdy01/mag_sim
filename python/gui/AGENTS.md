# Flask GUI Agent Notes

- This package hosts the Flask-based prototype web UI for `mag_sim`.
- Keep browser dependencies minimal: prefer server-sent events (SSE) and
  lightweight Bootstrap styling already included in the templates. Avoid adding
  large frontend frameworks unless explicitly requested.
- Reuse the existing `SimulationManager` helper when wiring new routes. Ensure
  updates remain thread-safe and continue to enforce the single-simulation
  constraint.
- Stage 1 of the GUI overhaul introduced geometry previews, live field-map
  visualisation, and the `/visualization.png` customisation endpoint. Any
  follow-up work should keep those panels functional and expand tests when
  touching related code.
- When modifying behaviour, update the user guide (`docs/user-guide/gui_flask.md`)
  and extend `tests/test_gui_flask.py` so new code paths stay covered.
- Runtime artefacts such as uploads and generated logs should stay out of the
  repository. Use or extend the ignore rules in this directory if new folders
  are introduced.
