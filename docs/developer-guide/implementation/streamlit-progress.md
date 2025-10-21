# Streamlit GUI Implementation Progress

This log captures the major milestones while building the Streamlit-based GUI.
It complements the solver progress log and helps future contributors understand
what remains before feature parity with the CLI tools.

## Checklist

- [x] Add Streamlit, `numpy`, `matplotlib`, and `ezdxf` to the development
      environment (`.devcontainer/Dockerfile`, `scripts/setup_env.sh`, and
      `requirements-gui.txt`).
- [x] Scaffold the `python/gui` package with `app_streamlit.py`, helper tests,
      and agent notes.
- [x] Implement sidebar controls (file upload, solver configuration, scenario
      download) and main-panel placeholders (geometry preview, status, logs,
      progress bar, and output downloads).
- [x] Wire the Streamlit UI to `motor_sim` via `subprocess`, including live log
      streaming and a Stop button backed by session state.
- [x] Provide a headless `--test-run` mode for smoke testing and CI integration.
- [x] Document usage in `docs/user-guide/gui_streamlit.md` and reference the
      alternative Flask branch.
- [x] Add unit tests covering scenario staging, metadata generation, process log
      streaming, and the test hook.
- [ ] Visual geometry previews and DXF-to-scenario conversion pipeline.
- [ ] Richer progress reporting (parse residual percentages rather than counting
      log lines).
- [ ] In-app visualisation of solver outputs (e.g., field maps via matplotlib).

## Notes

- DXF uploads currently short-circuit after staging; once the geometry pipeline
  is ready, replace the placeholder message with conversion to scenario JSON and
  reuse the existing run loop.
- The Stop button terminates `motor_sim` by sending `terminate()` and relaying a
  log message. If the simulator adds a graceful shutdown flag, prefer that over
  termination to avoid partial output files.
- Output download links are built from the scenario JSON. Ensure new scenario
  schemas continue to populate `outputs[*].path` so the GUI can locate files.
- Keep long-running operations (mesh generation, solves) within the streaming
  loop to preserve responsiveness—avoid blocking the Streamlit script without
  emitting log messages.
