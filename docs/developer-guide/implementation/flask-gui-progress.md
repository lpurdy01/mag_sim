# Flask GUI Implementation Progress

This log tracks the development of the Flask-based web prototype for `mag_sim`.
It complements the existing implementation status and progress notes by focusing
on the GUI-specific milestones.

## Completed work

- Provisioned the devcontainer with Flask, Flask-CORS, NumPy, Matplotlib, and
  `ezdxf` so the GUI can reuse plotting or DXF parsing helpers when they become
  available.
- Established a Flask application (`python/gui/app_flask.py`) with:
  - A background `SimulationManager` that keeps the `motor_sim` subprocess
    isolated from request threads and streams stdout via server-sent events.
  - Upload handling for scenario JSON files with validation and command-line
    flag mapping for solver, tolerance, max iterations, and optional outputs.
  - Download endpoints for the uploaded scenario and live log artefacts.
- Added Bootstrap-backed templates that expose upload controls, progress bar,
  log viewer, and download list with minimal JavaScript (EventSource + DOM
  updates only).
- Wrote unit tests that exercise upload flow, stop signalling, and the SSE
  endpoint payloads without invoking the real solver binary.
- Documented usage (`docs/user-guide/gui_flask.md`) and linked the new page into
  the MkDocs navigation.
- **Stage 1 UI overhaul**:
  - Split the main page into dedicated geometry, control, progress, downloads,
    and visualisation panels with placeholders for future DXF/animation work.
  - Introduced a CAD-style top bar and workflow sidebar so navigation feels
    closer to desktop tooling while leaving hooks for future menus.
  - Added geometry preview rendering on demand plus live in-process and final
    field-map visualisation streamed via SSE.
  - Introduced configurable plotting controls (vector scaling, density,
    overlays, logarithmic modes) that regenerate Matplotlib renders on request.
  - Persisted scenario uploads as per-session projects, enabling preview and run
    actions without re-uploading; added File-menu routes for exporting or
    clearing the workspace.
  - Expanded automated coverage for preview flows, solver output detection,
    project reuse, and the visualisation endpoint to keep regression protection
    high.
  - Added helper scripts (`scripts/setup_gui_env.sh` and
    `scripts/maintain_gui_env.sh`) so contributors can bootstrap dependencies,
    prune artefacts, and monitor Python package freshness locally.

## In-flight / next steps

- DXF ingestion is still stubbed. Once the geometry-to-scenario converter is
  ready we can replace the placeholder guard in `/upload` with real processing.
- The progress parser currently looks for percentage tokens in stdout; refining
  this to understand the solver's exact log format will improve the progress bar
  fidelity.
- Surface more solver artefacts (line probes, streamlines) alongside field maps
  once those exporters are exercised in real scenarios.
- Harden long-running process management for multiple users (per-session queues
  instead of global state) if the GUI graduates beyond a single-user demo.
- Add screenshot-driven regression tests for the Bootstrap layout when we wire
  up a lightweight browser harness.
