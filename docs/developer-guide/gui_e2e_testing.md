# GUI end-to-end testing

This guide explains how to exercise the Flask GUI in a real browser using
Playwright. The suite drives the CAD upload, scenario load, and run workflow so
regressions in the interactive flow are caught early.

## What is covered

The Playwright test in `tests/e2e/test_gui_e2e.py` replicates the induction
motor walkthrough:

1. Start a clean project.
2. Upload domain and rotor-bar DXFs, mark the layers as included, and save the
   mapping so the combined preview refreshes.
3. Upload `inputs/induction_gui_demo.json` in **Simulation setup**.
4. Click **Run simulation** and watch the progress/log cards populate.
5. Confirm the downloads list exposes the staged field map and the
   visualisation panel renders the frame.

The test uses the template `data-testid` hooks to keep selectors stable even if
the layout evolves.

## Environment setup

Run the helper script to install the browser automation dependencies and the
GUI stack:

```bash
./scripts/setup_gui_e2e_env.sh
```

This installs Flask, Matplotlib, ezdxf, pytest, `pytest-playwright`, and the
Playwright Chromium binary so CI or local runs have everything required.

## Running the suite

Launch the Playwright-driven test from the repository root:

```bash
pytest tests/e2e/test_gui_e2e.py
```

The fixture spins up the Flask app on a random local port, stubs the solver
with a fast in-process field-map emission, and exercises the UI in a headless
Chromium instance. Screenshots can be added to the test if debugging future
regressions.
