# Web GUI (Streamlit) Usage

The Streamlit interface provides a lightweight browser-based front-end for
running `motor_sim` without leaving the development environment. It is designed
for quick iteration: upload a scenario, tweak solver options, and watch the
progress log update live as the solver runs.

## Prerequisites

The development container and setup script install the required Python
dependencies (`streamlit`, `numpy`, `matplotlib`, and `ezdxf`). If you are
working outside the container, install them manually:

```bash
python3 -m pip install -r requirements-gui.txt
```

Compile the simulator before launching the GUI:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

The GUI expects the `build/motor_sim` binary to be present.

## Launching the app

From the repository root run:

```bash
streamlit run python/gui/app_streamlit.py
```

Streamlit listens on port `8501` by default. In VS Code Remote or GitHub
Codespaces forward that port to your browser. When the page loads you will see
an idle progress bar, a placeholder geometry preview, and an empty log panel.

## Sidebar controls

The sidebar contains everything required to stage a simulation:

- **File uploader** – Accepts JSON scenarios or raw DXF geometry. JSON files are
  passed straight to the solver. DXF upload support is staged for future work;
  the GUI records the file and surfaces a reminder that parsing is not yet
  available.
- **Solver configuration** – Choose between the SOR and CG solvers, set the
  convergence tolerance, maximum iterations, and the progress sampling cadence.
  These values are forwarded as CLI flags when the solver launches.
- **Scenario download** – Press *Download Scenario JSON* to export the current
  configuration (including solver options) as a JSON file. This is useful for
  versioning or running the same configuration from the command line.
- **Region selection note** – A short reminder marks the spot where future
  multi-region controls will appear once DXF ingestion lands.

## Main panel walkthrough

1. **Geometry preview** – Indicates whether a scenario JSON or DXF file is in
   use. Visualisation hooks are stubbed for now and will be populated once the
   geometry pipeline is integrated.
2. **Run/Stop buttons** – Use *Run Simulation* to start a solve. The button is
   disabled while a simulation is running to prevent duplicate launches. The
   *Stop Simulation* button sets a session flag that terminates the solver
   process on the next log update.
3. **Status & progress** – A live progress bar increments with each log line and
   resets when runs complete or abort. Status text summarises the most recent
   action (starting, stopping, completion, or failure).
4. **Log panel** – Streams combined stdout/stderr from `motor_sim`. Messages are
   appended as the solver runs so you can watch residuals and diagnostics in
   real time.
5. **Outputs** – After a successful run the GUI scans the scenario's `outputs`
   section and renders download buttons for each file that exists on disk. This
   is ideal for quickly grabbing CSVs for plotting.

If the solver terminates early (either via the Stop button or due to an error)
all state flags reset and the progress bar returns to zero. Errors surface in
red above the log panel for quick diagnosis.

## Limitations and tips

- DXF support is intentionally minimal in this iteration. Uploads are preserved
  for later conversion but are not yet passed through to the solver.
- Ensure that the simulator has been built before launching the GUI. The app
  reports a helpful error if `build/motor_sim` is missing.
- Long-running simulations still stream logs, but Streamlit queues updates. If
  you need higher fidelity residual history, also enable `--progress-history`
  in your scenario and download the CSV afterwards.
- The download buttons read files on each render; large artefacts should remain
  outside of Git history in keeping with repository guidelines.

## Alternate implementation

A Flask-based prototype lives on the `feat/gui-flask` branch. Check out that
branch if you would like to compare approaches:

```bash
git fetch origin
git checkout feat/gui-flask
```

Switch back to the current Streamlit implementation with `git checkout main`
(or the branch you are developing on). The Streamlit UI provides the simplest
on-ramp for experimentation and is the recommended path going forward.
