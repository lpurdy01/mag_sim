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

1. Upload a scenario JSON file via the **Scenario file** input. (DXF uploads are
   reserved for a future update and currently return a validation error.)
2. Pick a solver (`cg` or `sor`), adjust the tolerance and maximum iteration
   count if required, and optionally pass a value for `--outputs`.
3. Click **Run simulation**. The form disables itself, the **Stop** button
   becomes available, and a progress card appears.
4. The server streams the solver's stdout via Server-Sent Events (SSE). The
   progress bar reacts to any percentage tokens in the log, and the log panel
   fills with the live output.
5. When the run finishes, the log persists and a **Downloads** card appears with
   links to the uploaded scenario and the captured log file. The log is useful
   for sharing solver traces without leaving the browser.
6. Click **Run simulation** again to process another scenario. Only one solve is
   allowed at a time; the interface reports an error if you attempt to start a
   second run before the first completes.

Use the **Stop** button to request early termination. The server sends
`terminate()` to the subprocess and the final log entry marks the stop request.
Depending on solver state, expect a short delay while the child process exits.

## Known limitations

- DXF geometry ingestion is stubbed; use JSON scenarios generated via the
  existing Python helpers until the converter is available.
- The progress parser uses a simple percentage heuristic. If the solver log does
  not emit percentage tokens the bar will stay at 0 % until the run completes.
- The current prototype keeps global state and targets single-user usage. A
  production deployment should move to per-session queues and authentication.

## Related work

A Streamlit-based GUI prototype lives on the `feat/gui-streamlit` branch. Both
paths explore similar workflows; use this Flask version when you need explicit
control over routing and templating, or switch branches to compare the
Streamlit-powered experience.
