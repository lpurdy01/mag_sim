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

## In-flight / next steps

- DXF ingestion is still stubbed. Once the geometry-to-scenario converter is
  ready we can replace the placeholder guard in `/upload` with real processing.
- The progress parser currently looks for percentage tokens in stdout; refining
  this to understand the solver's exact log format will improve the progress bar
  fidelity.
- Consider persisting additional solver outputs (e.g., field map CSVs) to the
  downloads list once the CLI contract is finalised.
- Harden long-running process management for multiple users (per-session queues
  instead of global state) if the GUI graduates beyond a single-user demo.
