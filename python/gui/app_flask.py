"""Flask application providing a lightweight web GUI for ``mag_sim``."""

from __future__ import annotations

import json
import queue
import re
import subprocess
import threading
import time
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

from flask import (
    Flask,
    Response,
    abort,
    redirect,
    render_template,
    request,
    send_file,
    session,
    url_for,
)
from werkzeug.utils import secure_filename

from python.gui.render import (
    DEFAULT_LOG_FLOOR,
    render_field_map_image,
    render_geometry_preview,
)


ALLOWED_EXTENSIONS = {".json", ".dxf"}
DEFAULT_SOLVER = "cg"
DEFAULT_TOLERANCE = 1e-6
DEFAULT_MAX_ITERS = 20000
MAX_UPLOAD_SIZE_MB = 50
DEFAULT_VECTOR_MODE = "linear"
DEFAULT_COLOR_SCALE = "linear"
DEFAULT_QUIVER_SKIP = 4
DEFAULT_STREAMLINES = False
FIELD_MAP_REGEX = re.compile(
    r"Frame\\s+(?P<frame>\\d+):\\s+wrote\\s+field_map\\s+'(?P<id>[^']+)'\\s+to\\s+(?P<path>.+)"
)


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
        self._last_result: Dict[str, Any] = {}

    @property
    def is_running(self) -> bool:
        return self._running

    @property
    def queue(self) -> "queue.Queue[Dict[str, Any]]":
        return self._queue

    def get_last_result(self) -> Dict[str, Any]:
        with self._lock:
            return dict(self._last_result)

    def start(
        self,
        command: List[str],
        *,
        scenario_path: Path,
        log_path: Path,
        extra_downloads: Optional[List[Dict[str, Any]]] = None,
        field_outputs: Optional[List[Dict[str, Any]]] = None,
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
                "field_outputs": field_outputs or [],
                "processed_field_maps": set(),
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
            self._last_result = {}

    def clear_last_result(self) -> None:
        """Discard the cached visualisation payload."""

        with self._lock:
            self._last_result = {}

    def _drain_queue_locked(self) -> None:
        while True:
            try:
                self._queue.get_nowait()
            except queue.Empty:
                break

    def _handle_field_map_event(self, *, frame: Optional[int], field_id: str, path_str: str) -> None:
        scenario_path: Optional[Path] = self._metadata.get("scenario_path")
        if scenario_path is None or not scenario_path.exists():
            return

        cleaned_path = path_str.strip().strip('"').strip("'")
        if not cleaned_path:
            return

        output_path = Path(cleaned_path)
        if not output_path.is_absolute():
            output_path = scenario_path.parent / output_path

        if not output_path.exists():
            return

        with self._lock:
            processed: Set[Path] = self._metadata.setdefault("processed_field_maps", set())
            if output_path in processed:
                return
            processed.add(output_path)

        results_dir = Path(self._app.config["RESULTS_FOLDER"])
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        frame_suffix = f"_frame{frame}" if frame is not None else ""
        filename = f"field_progress_{timestamp}{frame_suffix}.png"
        image_path = results_dir / filename

        try:
            render_field_map_image(
                scenario_path,
                output_path,
                output_path=image_path,
                vector_mode=DEFAULT_VECTOR_MODE,
                quiver_skip=DEFAULT_QUIVER_SKIP,
                color_scale=DEFAULT_COLOR_SCALE,
                draw_boundaries=True,
                streamlines=DEFAULT_STREAMLINES,
            )
        except Exception as exc:  # pragma: no cover - defensive guard
            self._queue.put({
                "message": f"Failed to render field map '{field_id}': {exc}",
                "warning": True,
            })
            return

        with self._lock:
            self._metadata["latest_field_map"] = output_path
            self._metadata["latest_image"] = image_path

        caption = f"Frame {frame} field map" if frame is not None else f"Field map '{field_id}'"
        self._queue.put(
            {
                "visualization": {
                    "image": f"/result-image/{image_path.name}",
                    "caption": caption,
                }
            }
        )

    def _candidate_field_maps(self) -> List[Path]:
        scenario_path: Optional[Path] = self._metadata.get("scenario_path")
        scenario_dir = scenario_path.parent if scenario_path else None
        candidates: List[Path] = []
        latest = self._metadata.get("latest_field_map")
        if isinstance(latest, Path):
            candidates.append(latest)
        for entry in self._metadata.get("field_outputs", []):
            path_value = entry.get("path")
            if not path_value:
                continue
            path = Path(path_value)
            if not path.is_absolute() and scenario_dir is not None:
                path = scenario_dir / path
            candidates.append(path)
        return candidates

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

                    match = FIELD_MAP_REGEX.match(line)
                    if match:
                        frame_value = match.group("frame")
                        frame = int(frame_value) if frame_value is not None else None
                        field_id = match.group("id")
                        path_str = match.group("path")
                        self._handle_field_map_event(
                            frame=frame,
                            field_id=field_id,
                            path_str=path_str,
                        )

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

        visualization_payload: Optional[Dict[str, Any]] = None
        if success:
            for candidate in self._candidate_field_maps():
                if not candidate.exists():
                    continue
                scenario_path_resolved = scenario_path if scenario_path.exists() else None
                results_dir = Path(self._app.config["RESULTS_FOLDER"])
                timestamp = time.strftime("%Y%m%d-%H%M%S")
                image_path = results_dir / f"field_result_{timestamp}.png"
                try:
                    render_field_map_image(
                        scenario_path_resolved or scenario_path,
                        candidate,
                        output_path=image_path,
                        vector_mode=DEFAULT_VECTOR_MODE,
                        quiver_skip=DEFAULT_QUIVER_SKIP,
                        color_scale=DEFAULT_COLOR_SCALE,
                        draw_boundaries=True,
                        streamlines=DEFAULT_STREAMLINES,
                    )
                except Exception as exc:  # pragma: no cover - defensive guard
                    self._queue.put({
                        "message": f"Failed to render final field map: {exc}",
                        "warning": True,
                    })
                    break

                downloads.append(
                    {
                        "category": "result",
                        "filename": image_path.name,
                        "label": "Field map image",
                    }
                )
                visualization_payload = {
                    "image": f"/result-image/{image_path.name}",
                    "caption": "Final field map",
                }
                with self._lock:
                    self._last_result = {
                        "scenario_path": str(scenario_path),
                        "field_map": str(candidate),
                        "image": str(image_path),
                        "caption": "Final field map",
                        "settings": {
                            "vector_mode": DEFAULT_VECTOR_MODE,
                            "color_scale": DEFAULT_COLOR_SCALE,
                            "quiver_skip": DEFAULT_QUIVER_SKIP,
                            "draw_boundaries": True,
                            "streamlines": DEFAULT_STREAMLINES,
                        },
                    }
                break
            else:
                with self._lock:
                    self._last_result = {"message": "No field map output generated."}
                visualization_payload = {"message": "No field map output generated."}
        else:
            with self._lock:
                self._last_result = {}

        event: Dict[str, Any] = {"complete": True, "message": message, "success": success}
        if success:
            event["progress"] = 100
        if self._stop_requested:
            event["stopped"] = True
        if downloads:
            event["downloads"] = downloads
        if visualization_payload:
            event["visualization"] = visualization_payload

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

app.secret_key = "mag-sim-gui"

manager = SimulationManager(app)

PROJECTS: Dict[str, Dict[str, Any]] = {}


def allowed_file(filename: str) -> bool:
    return Path(filename).suffix.lower() in ALLOWED_EXTENSIONS


def _build_form_state(overrides: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    state = {
        "solver": DEFAULT_SOLVER,
        "tol": f"{DEFAULT_TOLERANCE}",
        "max_iters": f"{DEFAULT_MAX_ITERS}",
        "outputs": "",
    }
    if overrides:
        for key, value in overrides.items():
            if key in state and value is not None:
                state[key] = value
    return state


def _register_project(
    scenario_path: Path, spec: Dict[str, Any], original_name: str
) -> Tuple[str, Dict[str, Any]]:
    project_id = uuid.uuid4().hex
    PROJECTS[project_id] = {
        "scenario_path": scenario_path,
        "spec": spec,
        "original_filename": original_name,
        "scenario_name": spec.get("name"),
        "created_label": time.strftime("%Y-%m-%d %H:%M:%S"),
        "preview_path": None,
        "preview_notice": None,
        "preview_version": None,
        "form_state": _build_form_state(),
    }
    return project_id, PROJECTS[project_id]


def _discard_project(project_id: str) -> None:
    project = PROJECTS.pop(project_id, None)
    if not project:
        return

    preview_path = project.get("preview_path")
    if isinstance(preview_path, Path):
        preview_path.unlink(missing_ok=True)

    scenario_path = project.get("scenario_path")
    if isinstance(scenario_path, Path):
        scenario_path.unlink(missing_ok=True)


def _current_project(
    explicit_id: Optional[str] = None,
) -> Tuple[Optional[str], Optional[Dict[str, Any]]]:
    project_id = explicit_id or session.get("project_id")
    if not project_id:
        return None, None

    project = PROJECTS.get(project_id)
    if project is None:
        if explicit_id is None:
            session.pop("project_id", None)
        return None, None

    return project_id, project


def _index_context(
    *,
    error: Optional[str] = None,
    preview_url: Optional[str] = None,
    preview_error: Optional[str] = None,
    preview_notice: Optional[str] = None,
    form_state: Optional[Dict[str, str]] = None,
    project_id: Optional[str] = None,
    status_message: Optional[str] = None,
) -> Dict[str, Any]:
    active_id, project = _current_project(project_id)

    if form_state is None:
        if project and project.get("form_state"):
            form_state = dict(project["form_state"])
        else:
            form_state = _build_form_state()
    else:
        form_state = _build_form_state(form_state)

    preview_source = preview_url
    preview_token = None
    notice = preview_notice
    if project:
        stored_preview = project.get("preview_path")
        if stored_preview and Path(stored_preview).exists():
            if not preview_source:
                preview_source = url_for("result_image", filename=Path(stored_preview).name)
            preview_token = project.get("preview_version")
        if not notice and project.get("preview_notice"):
            notice = project["preview_notice"]

    last_result = manager.get_last_result()
    last_image_url = None
    visualization_caption = None
    visualization_message = last_result.get("message") if last_result else None
    visualization_defaults = {
        "vector_mode": DEFAULT_VECTOR_MODE,
        "color_scale": DEFAULT_COLOR_SCALE,
        "quiver_skip": DEFAULT_QUIVER_SKIP,
        "draw_boundaries": True,
        "streamlines": DEFAULT_STREAMLINES,
    }
    if last_result:
        image_path = Path(last_result.get("image", ""))
        if image_path.exists():
            last_image_url = url_for("result_image", filename=image_path.name)
            visualization_caption = last_result.get("caption")
            visualization_message = None
        settings = last_result.get("settings", {})
        for key in visualization_defaults:
            if key in settings:
                visualization_defaults[key] = settings[key]

    project_summary = None
    if project and active_id:
        project_summary = {
            "id": active_id,
            "filename": project.get("original_filename"),
            "scenario_name": project.get("scenario_name"),
            "created_label": project.get("created_label"),
        }

    return {
        "running": manager.is_running,
        "error": error,
        "form_state": form_state,
        "preview_url": preview_source,
        "preview_error": preview_error,
        "preview_notice": notice,
        "preview_token": preview_token,
        "visualization_defaults": visualization_defaults,
        "last_visualization_url": last_image_url,
        "last_visualization_caption": visualization_caption,
        "visualization_message": visualization_message,
        "visualization_visible": bool(last_result),
        "DEFAULT_LOG_FLOOR": DEFAULT_LOG_FLOOR,
        "project": project_summary,
        "has_project": project_summary is not None,
        "project_download_url": url_for("project_scenario") if project_summary else None,
        "form_requires_file": project_summary is None,
        "status_message": status_message,
    }


def _store_scenario_file(geometry) -> Tuple[Path, dict, str, str]:
    suffix = Path(geometry.filename).suffix.lower()
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    original_name = geometry.filename or f"scenario_{timestamp}{suffix}"
    safe_name = secure_filename(original_name) or f"scenario_{timestamp}{suffix}"
    scenario_path = Path(app.config["UPLOAD_FOLDER"]) / f"{timestamp}_{safe_name}"

    scenario_path.parent.mkdir(parents=True, exist_ok=True)
    geometry.save(scenario_path)
    try:
        spec = json.loads(scenario_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        scenario_path.unlink(missing_ok=True)
        raise ValueError("Uploaded JSON is invalid.") from exc
    return scenario_path, spec, timestamp, original_name


def _collect_field_outputs(spec: dict, scenario_dir: Path) -> List[Dict[str, Any]]:
    outputs: List[Dict[str, Any]] = []
    for entry in spec.get("outputs", []):
        if entry.get("type") != "field_map":
            continue
        identifier = entry.get("id")
        path_value = entry.get("path")
        if not path_value and identifier:
            path_value = f"outputs/{identifier}.csv"
        if not path_value:
            continue
        path = Path(path_value)
        if not path.is_absolute():
            path = scenario_dir / path
        outputs.append({"id": identifier, "path": path})
    return outputs


@app.get("/")
def index() -> str:
    error = request.args.get("error")
    notice = request.args.get("notice")
    return render_template(
        "index.html",
        **_index_context(error=error, status_message=notice),
    )


@app.post("/upload")
def upload() -> Response:
    action = (request.form.get("action") or "run").lower()
    if manager.is_running:
        return redirect(url_for("index", error="A simulation is already running."))

    geometry = request.files.get("geometry_file")
    project_hint = request.form.get("project_id")
    current_project_id = session.get("project_id")

    project_id: Optional[str] = None
    project: Optional[Dict[str, Any]] = None
    scenario_path: Optional[Path] = None
    spec: Optional[Dict[str, Any]] = None

    if geometry is not None and geometry.filename:
        if not allowed_file(geometry.filename):
            return redirect(url_for("index", error="Unsupported file type. Use JSON scenarios."))

        suffix = Path(geometry.filename).suffix.lower()
        if suffix != ".json":
            return redirect(url_for("index", error="DXF uploads are not supported yet."))

        if current_project_id and current_project_id in PROJECTS:
            _discard_project(current_project_id)

        try:
            scenario_path, spec, _, original_name = _store_scenario_file(geometry)
        except ValueError as exc:
            return redirect(url_for("index", error=str(exc)))

        project_id, project = _register_project(scenario_path, spec, original_name)
        session["project_id"] = project_id
    else:
        project_id, project = _current_project(project_hint)
        if project is None:
            return redirect(url_for("index", error="Select or upload a scenario before continuing."))

        scenario_path = project.get("scenario_path")
        spec = project.get("spec")
        if not isinstance(scenario_path, Path) or not scenario_path.exists():
            return redirect(url_for("index", error="Scenario file is missing for this project."))
        if not isinstance(spec, dict):
            return redirect(url_for("index", error="Scenario metadata is unavailable for this project."))

    assert project_id is not None and project is not None
    session["project_id"] = project_id
    if scenario_path is None:
        scenario_path = project["scenario_path"]
    if spec is None:
        spec = project["spec"]

    form_state = _build_form_state(
        {
            "solver": request.form.get("solver"),
            "tol": request.form.get("tol"),
            "max_iters": request.form.get("max_iters"),
            "outputs": request.form.get("outputs"),
        }
    )
    project["form_state"] = dict(form_state)
    project["spec"] = spec
    project["scenario_path"] = scenario_path
    project["scenario_name"] = spec.get("name")

    if action == "preview":
        preview_path_existing = project.get("preview_path")
        if isinstance(preview_path_existing, Path):
            preview_path_existing.unlink(missing_ok=True)

        preview_filename = f"geometry_preview_{project_id}_{time.strftime('%Y%m%d-%H%M%S')}.png"
        preview_path = Path(app.config["RESULTS_FOLDER"]) / preview_filename
        try:
            render_geometry_preview(scenario_path, preview_path)
        except RuntimeError as exc:
            project["preview_path"] = None
            project["preview_version"] = None
            project["preview_notice"] = None
            return render_template(
                "index.html",
                **_index_context(
                    preview_error=str(exc),
                    form_state=form_state,
                    project_id=project_id,
                ),
            )

        project["preview_path"] = preview_path
        project["preview_notice"] = (
            f"Preview generated from {project.get('original_filename') or scenario_path.name}"
        )
        project["preview_version"] = time.time()

        return render_template(
            "index.html",
            **_index_context(
                preview_notice=project["preview_notice"],
                form_state=form_state,
                project_id=project_id,
            ),
        )

    solver = (form_state["solver"] or DEFAULT_SOLVER).lower()
    if solver not in {"cg", "sor"}:
        return redirect(url_for("index", error="Solver must be either CG or SOR."))

    try:
        tolerance = float(form_state["tol"])
    except ValueError:
        return redirect(url_for("index", error="Tolerance must be numeric."))

    try:
        max_iters = int(form_state["max_iters"])
    except ValueError:
        return redirect(url_for("index", error="Max iterations must be an integer."))

    outputs_value = form_state["outputs"].strip()

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

    run_stamp = time.strftime("%Y%m%d-%H%M%S")
    log_path = Path(app.config["RESULTS_FOLDER"]) / f"simulation_{run_stamp}.log"

    field_outputs = _collect_field_outputs(spec, scenario_path.parent)
    extra_downloads = [
        {
            "category": "result",
            "path": str(entry["path"]),
            "label": f"Field map CSV ({entry['id']})" if entry.get("id") else entry["path"].name,
        }
        for entry in field_outputs
    ]

    try:
        manager.start(
            command,
            scenario_path=scenario_path,
            log_path=log_path,
            extra_downloads=extra_downloads,
            field_outputs=[{"id": entry.get("id"), "path": str(entry["path"])} for entry in field_outputs],
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


@app.post("/project/reset")
def project_reset() -> Response:
    if manager.is_running:
        return redirect(url_for("index", error="Stop the running simulation before clearing the project."))

    project_id = session.pop("project_id", None)
    if project_id:
        _discard_project(project_id)

    manager.clear_last_result()
    return redirect(url_for("index", notice="Project cleared."))


@app.get("/project/scenario")
def project_scenario() -> Response:
    project_id, project = _current_project()
    if project is None:
        abort(404)

    scenario_path = project.get("scenario_path")
    if not isinstance(scenario_path, Path) or not scenario_path.exists():
        abort(404)

    download_name = project.get("original_filename") or scenario_path.name
    return send_file(scenario_path, as_attachment=True, download_name=download_name)


@app.get("/result-image/<path:filename>")
def result_image(filename: str) -> Response:
    safe_name = secure_filename(filename)
    if not safe_name:
        abort(400)

    file_path = Path(app.config["RESULTS_FOLDER"]) / safe_name
    if not file_path.exists():
        abort(404)

    return send_file(file_path, mimetype="image/png")


@app.get("/visualization.png")
def visualization_image() -> Response:
    last_result = manager.get_last_result()
    scenario_path = Path(last_result.get("scenario_path", ""))
    field_map_path = Path(last_result.get("field_map", ""))
    if not scenario_path.exists() or not field_map_path.exists():
        abort(404)

    vector_mode = request.args.get("vector", DEFAULT_VECTOR_MODE)
    if vector_mode not in {"linear", "log", "off"}:
        abort(400)

    try:
        quiver_skip = int(request.args.get("skip", DEFAULT_QUIVER_SKIP))
        if quiver_skip < 1:
            raise ValueError
    except ValueError:
        abort(400)

    color_scale = request.args.get("scale", DEFAULT_COLOR_SCALE)
    if color_scale not in {"linear", "log"}:
        abort(400)

    draw_boundaries = request.args.get("boundaries", "1") != "0"
    streamlines = request.args.get("stream", "0") == "1"

    try:
        log_floor_value = float(request.args.get("log_floor", DEFAULT_LOG_FLOOR))
        vector_floor_value = float(request.args.get("vector_floor", DEFAULT_LOG_FLOOR))
    except ValueError:
        abort(400)

    data = render_field_map_image(
        scenario_path,
        field_map_path,
        output_path=None,
        vector_mode=vector_mode,
        quiver_skip=quiver_skip,
        color_scale=color_scale,
        draw_boundaries=draw_boundaries,
        streamlines=streamlines,
        log_floor=log_floor_value,
        vector_log_floor=vector_floor_value,
    )

    if isinstance(data, Path):  # pragma: no cover - defensive guard
        payload = data.read_bytes()
    else:
        payload = data

    return Response(payload, mimetype="image/png")


if __name__ == "__main__":  # pragma: no cover - manual execution helper
    app.run(host="0.0.0.0", port=5000, debug=True)
