"""Flask application providing a lightweight web GUI for ``mag_sim``."""

from __future__ import annotations

import json
import queue
import subprocess
import threading
import time
from pathlib import Path
from typing import Any, Dict, List, Optional

from flask import (
    Flask,
    Response,
    abort,
    redirect,
    render_template,
    request,
    url_for,
    send_file,
)
from werkzeug.utils import secure_filename


ALLOWED_EXTENSIONS = {".json", ".dxf"}
DEFAULT_SOLVER = "cg"
DEFAULT_TOLERANCE = 1e-6
DEFAULT_MAX_ITERS = 20000
MAX_UPLOAD_SIZE_MB = 50


class SimulationManager:
    """Manage the lifecycle of a single background simulation process."""

    def __init__(self, flask_app: Flask) -> None:
        self._app = flask_app
        self._lock = threading.Lock()
        self._queue: "queue.Queue[Dict[str, Any]]" = queue.Queue()
        self._process: Optional[subprocess.Popen[str]] = None
        self._thread: Optional[threading.Thread] = None
        self._running = False
        self._stop_requested = False
        self._metadata: Dict[str, Any] = {}

    @property
    def is_running(self) -> bool:
        return self._running

    @property
    def queue(self) -> "queue.Queue[Dict[str, Any]]":
        return self._queue

    def start(
        self,
        command: List[str],
        *,
        scenario_path: Path,
        log_path: Path,
        extra_downloads: Optional[List[Dict[str, Any]]] = None,
    ) -> None:
        with self._lock:
            if self._running:
                raise RuntimeError("A simulation is already running.")
            self._drain_queue_locked()
            self._running = True
            self._stop_requested = False
            self._metadata = {
                "scenario_path": scenario_path,
                "log_path": log_path,
                "extra_downloads": extra_downloads or [],
            }
            self._queue.put({
                "started": True,
                "message": "Simulation launched.",
                "progress": 0,
            })
            thread = threading.Thread(
                target=self._run_process, args=(command,), daemon=True
            )
            self._thread = thread
            thread.start()

    def stop(self) -> bool:
        with self._lock:
            if self._process is None or not self._running:
                return False
            self._stop_requested = True
            self._queue.put({"message": "Termination requested…", "warning": True})
            try:
                self._process.terminate()
            except OSError:
                # Process already gone.
                pass
            return True

    def reset(self) -> None:
        """Test helper to clear internal state."""

        with self._lock:
            self._drain_queue_locked()
            self._process = None
            self._thread = None
            self._running = False
            self._stop_requested = False
            self._metadata = {}

    def _drain_queue_locked(self) -> None:
        while True:
            try:
                self._queue.get_nowait()
            except queue.Empty:
                break

    def _run_process(self, command: List[str]) -> None:
        log_path: Path = self._metadata["log_path"]
        log_path.parent.mkdir(parents=True, exist_ok=True)
        scenario_path: Path = self._metadata.get("scenario_path", log_path)

        success = False
        message = "Simulation did not start."

        try:
            with log_path.open("w", encoding="utf-8") as log_file:
                try:
                    process = subprocess.Popen(
                        command,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        text=True,
                        bufsize=1,
                    )
                except FileNotFoundError:
                    message = (
                        "Unable to start motor_sim; build the project before running the GUI."
                    )
                    self._queue.put({"message": message, "error": True})
                    return
                except Exception as exc:  # pragma: no cover - defensive guard
                    message = f"Failed to launch simulation: {exc}"
                    self._queue.put({"message": message, "error": True})
                    return

                with self._lock:
                    self._process = process

                assert process.stdout is not None
                for raw_line in process.stdout:
                    line = raw_line.rstrip()
                    log_file.write(line + "\n")
                    log_file.flush()
                    event: Dict[str, Any] = {"message": line}
                    progress = _extract_progress(line)
                    if progress is not None:
                        event["progress"] = progress
                    self._queue.put(event)

                return_code = process.wait()
                success = return_code == 0 and not self._stop_requested
                if self._stop_requested and return_code != 0:
                    # Some environments report 143 when terminate() is used; normalise.
                    success = False

                if self._stop_requested:
                    message = "Simulation stopped by user."
                elif success:
                    message = "Simulation complete."
                else:
                    message = f"Simulation exited with code {return_code}."
        finally:
            with self._lock:
                self._running = False
                self._process = None
            self._emit_completion(success=success, message=message, scenario_path=scenario_path)

    def _emit_completion(self, *, success: bool, message: str, scenario_path: Path) -> None:
        downloads = []
        log_path: Optional[Path] = self._metadata.get("log_path")
        if scenario_path and scenario_path.exists():
            downloads.append(
                {
                    "category": "scenario",
                    "filename": scenario_path.name,
                    "label": "Uploaded scenario",
                }
            )
        if log_path and Path(log_path).exists():
            downloads.append(
                {
                    "category": "log",
                    "filename": Path(log_path).name,
                    "label": "Simulation log",
                }
            )
        for extra in self._metadata.get("extra_downloads", []):
            path = Path(extra.get("path"))
            if not path.exists():
                continue
            downloads.append(
                {
                    "category": extra.get("category", "result"),
                    "filename": path.name,
                    "label": extra.get("label", path.name),
                }
            )

        event: Dict[str, Any] = {"complete": True, "message": message, "success": success}
        if success:
            event["progress"] = 100
        if self._stop_requested:
            event["stopped"] = True
        if downloads:
            event["downloads"] = downloads

        self._queue.put(event)
        self._stop_requested = False


def _extract_progress(line: str) -> Optional[int]:
    for token in line.replace("%", " % ").split():
        if token.endswith("%"):
            token = token[:-1]
        try:
            value = float(token)
        except ValueError:
            continue
        if 0 <= value <= 100:
            return int(value)
    return None


BASE_DIR = Path(__file__).resolve().parent
TEMPLATE_DIR = BASE_DIR / "templates"
STATIC_DIR = BASE_DIR / "static"

app = Flask(__name__, template_folder=str(TEMPLATE_DIR), static_folder=str(STATIC_DIR))
app.config.update(
    UPLOAD_FOLDER=str(BASE_DIR / "uploads"),
    RESULTS_FOLDER=str(BASE_DIR / "results"),
    MAX_CONTENT_LENGTH=MAX_UPLOAD_SIZE_MB * 1024 * 1024,
    JSON_SORT_KEYS=False,
)

Path(app.config["UPLOAD_FOLDER"]).mkdir(parents=True, exist_ok=True)
Path(app.config["RESULTS_FOLDER"]).mkdir(parents=True, exist_ok=True)

manager = SimulationManager(app)


def allowed_file(filename: str) -> bool:
    return Path(filename).suffix.lower() in ALLOWED_EXTENSIONS


@app.get("/")
def index() -> str:
    error = request.args.get("error")
    return render_template(
        "index.html",
        running=manager.is_running,
        error=error,
        default_solver=DEFAULT_SOLVER,
        default_tolerance=DEFAULT_TOLERANCE,
        default_max_iters=DEFAULT_MAX_ITERS,
    )


@app.post("/upload")
def upload() -> Response:
    if manager.is_running:
        return redirect(url_for("index", error="A simulation is already running."))

    geometry = request.files.get("geometry_file")
    if geometry is None or geometry.filename == "":
        return redirect(url_for("index", error="Select a scenario JSON before running."))

    if not allowed_file(geometry.filename):
        return redirect(url_for("index", error="Unsupported file type. Use JSON scenarios."))

    solver = (request.form.get("solver") or DEFAULT_SOLVER).lower()
    if solver not in {"cg", "sor"}:
        return redirect(url_for("index", error="Solver must be either CG or SOR."))

    tol_raw = request.form.get("tol", str(DEFAULT_TOLERANCE))
    try:
        tolerance = float(tol_raw)
    except ValueError:
        return redirect(url_for("index", error="Tolerance must be numeric."))

    max_iters_raw = request.form.get("max_iters", str(DEFAULT_MAX_ITERS))
    try:
        max_iters = int(max_iters_raw)
    except ValueError:
        return redirect(url_for("index", error="Max iterations must be an integer."))

    outputs_value = (request.form.get("outputs") or "").strip()

    suffix = Path(geometry.filename).suffix.lower()
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    safe_name = secure_filename(geometry.filename) or f"scenario_{timestamp}{suffix}"
    scenario_path = Path(app.config["UPLOAD_FOLDER"]) / f"{timestamp}_{safe_name}"

    if suffix == ".json":
        scenario_path.parent.mkdir(parents=True, exist_ok=True)
        geometry.save(scenario_path)
        try:
            json.loads(scenario_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            scenario_path.unlink(missing_ok=True)
            return redirect(url_for("index", error="Uploaded JSON is invalid."))
    else:
        # DXF support is not yet implemented for the Flask prototype.
        return redirect(url_for("index", error="DXF uploads are not supported yet."))

    command = [
        "./build/motor_sim",
        "--scenario",
        str(scenario_path),
        "--solver",
        solver,
        "--tol",
        str(tolerance),
        "--max-iters",
        str(max_iters),
        "--solve",
    ]

    if outputs_value:
        command.extend(["--outputs", outputs_value])

    log_path = Path(app.config["RESULTS_FOLDER"]) / f"simulation_{timestamp}.log"

    try:
        manager.start(
            command,
            scenario_path=scenario_path,
            log_path=log_path,
        )
    except RuntimeError as exc:
        return redirect(url_for("index", error=str(exc)))

    return redirect(url_for("index"))


@app.post("/stop")
def stop() -> Response:
    stopped = manager.stop()
    status_code = 202 if stopped else 200
    return ("", status_code)


@app.get("/progress")
def stream_progress() -> Response:
    def event_stream() -> Any:
        while True:
            try:
                payload = manager.queue.get(timeout=1.0)
            except queue.Empty:
                yield ": heartbeat\n\n"
                continue

            yield f"data: {json.dumps(payload)}\n\n"
            if payload.get("complete"):
                break

    return Response(event_stream(), mimetype="text/event-stream")


@app.get("/download")
def download() -> Response:
    category = request.args.get("category")
    filename = request.args.get("filename")
    if not category or not filename:
        abort(400)

    safe_name = secure_filename(filename)
    if not safe_name:
        abort(400)

    if category == "scenario":
        base_dir = Path(app.config["UPLOAD_FOLDER"])
    elif category in {"log", "result"}:
        base_dir = Path(app.config["RESULTS_FOLDER"])
    else:
        abort(404)

    file_path = base_dir / safe_name
    if not file_path.exists():
        abort(404)

    return send_file(file_path, as_attachment=True)


if __name__ == "__main__":  # pragma: no cover - manual execution helper
    app.run(host="0.0.0.0", port=5000, debug=True)
