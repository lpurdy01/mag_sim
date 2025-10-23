# Web GUI (Flask)

The Flask-based GUI provides a lightweight browser front-end for running
`mag_sim` scenarios. It mirrors the functionality explored in the Streamlit
prototype while staying close to the underlying command-line workflow.

## Requirements

- A built `motor_sim` executable in `./build/motor_sim`.
- Flask, Matplotlib, and NumPy available in the active Python environment.
  Execute `./scripts/setup_gui_env.sh` to install the trio (and create the GUI
  runtime folders) when a pre-provisioned devcontainer is not available.

## Starting the server

From the repository root:

```bash
python -m python.gui.app_flask
```

The development server listens on `http://127.0.0.1:5000` by default. When
running inside GitHub Codespaces or a similar environment, forward port 5000 to
access the interface from your browser.

Alternatively, set `FLASK_APP=python.gui.app_flask` and run `flask run` if you
prefer Flask's CLI wrapper.

## Using the interface

The CAD-inspired layout splits the page into a top navigation bar, a workflow
sidebar, and content panels for each stage. The **File** menu exposes
project-wide actions (**New project** clears the workspace, **Save scenario…**
downloads the current JSON), while the sidebar links jump to the key panels
without leaving the page.

1. Visit the **Simulation setup** panel to upload a scenario JSON file. If
   you're continuing work on an existing project, the hidden form state preserves
   the previously uploaded file—you only need to upload again when swapping to a
   different scenario. DXF uploads remain disabled until the converter ships.
2. Adjust solver parameters in the same panel. Solver (`cg` or `sor`),
   tolerance, maximum iterations, and an optional `--outputs` override all map
   directly to the CLI arguments used when launching `motor_sim`.
3. Use **Preview geometry** to render the domain outline without running the
   solver. The preview is cached with the project so you can revisit other
   panels without losing the image. Any new upload replaces the stored preview.
4. Click **Run simulation**. The run button disables itself, the **Stop** button
   becomes active, and the **Progress**/**Solver log** cards appear. The log is
   streamed via server‑sent events to avoid polling.
5. Watch the run unfold: percentage tokens update the progress bar, log entries
   append live, and any `field_map` outputs mentioned in stdout trigger a fresh
   render in the **Results visualisation** panel. Intermediate images replace
   the geometry preview automatically.
6. When the solver exits, the **Downloads** panel lists the uploaded scenario,
   the log file, and any field-map artefacts made available by the solver. The
   visualisation panel shows the final render and stays wired to the plotting
   controls so you can iterate on overlays without rerunning the simulation.

Only one solve is allowed at a time; the interface reports an error if you try
to launch a second run before the first finishes. Use **Stop** to request
termination—`terminate()` is sent to the subprocess and the log records the
request before the child exits.

### Environment upkeep

Generated artefacts accumulate under `python/gui/uploads` and
`python/gui/results`. Run `./scripts/maintain_gui_env.sh` periodically (set
`KEEP_DAYS` or pass a day count to keep more history) to prune stale files and
surface any outdated Python packages.

### Visualisation controls

The results card includes a lightweight form for regenerating the field-map
image without leaving the browser:

- **Vector scaling**: switch between linear arrows, logarithmic scaling, or hide
  vectors entirely.
- **Colour scale**: toggle between linear and logarithmic magnitude mapping.
- **Arrow skip**: control the quiver density (defaults to every fourth sample).
- **Log floor / Vector floor**: clamp low values when using logarithmic modes.
- **Region outlines / Streamlines**: enable or disable material/magnet borders
  and streamline overlays.

Choose the desired parameters and click **Update visualisation**. The server
re-renders the image with Matplotlib and returns a fresh PNG via the
`/visualization.png` endpoint, keeping the UI responsive with minimal
JavaScript.

## Known limitations

- DXF geometry ingestion is stubbed; use JSON scenarios generated via the
  existing Python helpers until the converter is available.
- Preview and visualisation rendering rely on Matplotlib and NumPy. These
  dependencies ship with the devcontainer, but a custom environment must supply
  the same stack.
- The progress parser uses a simple percentage heuristic. If the solver log does
  not emit percentage tokens the bar will stay at 0 % until the run completes.
- The current prototype keeps global state and targets single-user usage. A
  production deployment should move to per-session queues and authentication.

## Related work

A Streamlit-based GUI prototype lives on the `feat/gui-streamlit` branch. Both
paths explore similar workflows; use this Flask version when you need explicit
control over routing and templating, or switch branches to compare the
Streamlit-powered experience.
