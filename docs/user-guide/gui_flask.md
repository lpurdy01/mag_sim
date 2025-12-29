# Web GUI (Flask)

The Flask-based GUI provides a lightweight browser front-end for running
`mag_sim` scenarios. It mirrors the functionality explored in the Streamlit
prototype while staying close to the underlying command-line workflow.

## Requirements

- A built `motor_sim` executable in `./build/motor_sim`.
- The devcontainer (or local virtual environment) must include Flask and its
  lightweight companions. The repository's `.devcontainer/Dockerfile` installs
  Flask, Flask-CORS, NumPy, Matplotlib, and `ezdxf` so the GUI has everything it
  needs when launched inside Codespaces.

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

1. Upload a scenario JSON file via the **Scenario file** input. DXF uploads are
   acknowledged but remain disabled until the converter lands in a later stage.
2. Pick a solver (`cg` or `sor`), adjust the tolerance and maximum iteration
   count if required, and optionally pass a value for `--outputs`.
3. Use **Preview geometry** to render the domain outline without running the
   solver. The preview appears in the left-hand card and is regenerated on each
   upload.
4. Click **Run simulation**. The form disables itself, the **Stop** button
   becomes available, and progress/log cards appear underneath the controls.
5. While the solve runs the server streams stdout via SSE. The progress bar
   reacts to percentage tokens and the log window scrolls in real time. When the
   solver writes a field-map output the UI replaces the placeholder image with
   the rendered snapshot.
6. At completion a **Downloads** card lists the uploaded scenario, the log, and
   any field-map artefacts surfaced by the solver. The visualisation panel also
   shows the final render with a caption so you can inspect results immediately.

Only one solve is allowed at a time; the interface reports an error if you try
to launch a second run before the first finishes. Use **Stop** to request
termination—`terminate()` is sent to the subprocess and the log records the
request before the child exits.

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
