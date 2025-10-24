import io
import json
import queue
import sys
import time
from pathlib import Path
from typing import List

import pytest

pytest.importorskip("flask")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from python.gui import app_flask


@pytest.fixture(autouse=True)
def _reset_manager(tmp_path, monkeypatch):
    uploads = tmp_path / "uploads"
    results = tmp_path / "results"
    uploads.mkdir()
    results.mkdir()
    app_flask.app.config["UPLOAD_FOLDER"] = str(uploads)
    app_flask.app.config["RESULTS_FOLDER"] = str(results)
    app_flask.manager.reset()
    app_flask.PROJECTS.clear()
    yield
    app_flask.manager.reset()
    app_flask.PROJECTS.clear()


@pytest.fixture
def client():
    with app_flask.app.test_client() as test_client:
        yield test_client


@pytest.fixture
def visualization_ready(tmp_path):
    scenario = tmp_path / "viz_scenario.json"
    scenario.write_text(
        json.dumps(
            {
                "version": "0.2",
                "name": "visualisation fixture",
                "domain": {"Lx": 0.1, "Ly": 0.1},
                "sources": [],
            }
        ),
        encoding="utf-8",
    )

    field_map = tmp_path / "outputs" / "fixture_field.csv"
    field_map.parent.mkdir(parents=True, exist_ok=True)
    field_map.write_text(
        "x,y,Bx,By,Bmag\n"
        "0.0,0.0,0.1,0.0,0.1\n"
        "0.0,0.1,0.1,0.0,0.1\n"
        "0.1,0.0,0.1,0.0,0.1\n"
        "0.1,0.1,0.1,0.0,0.1\n",
        encoding="utf-8",
    )

    app_flask.manager._last_result = {  # type: ignore[attr-defined]
        "scenario_path": str(scenario),
        "field_map": str(field_map),
        "settings": {
            "vector_mode": app_flask.DEFAULT_VECTOR_MODE,
            "color_scale": app_flask.DEFAULT_COLOR_SCALE,
            "quiver_skip": app_flask.DEFAULT_QUIVER_SKIP,
            "draw_boundaries": True,
            "streamlines": False,
        },
    }

    return {"scenario": scenario, "field_map": field_map}


def test_upload_starts_simulation(monkeypatch, client):
    popen_calls: List[List[str]] = []

    class DummyProcess:
        def __init__(self) -> None:
            self.stdout = io.StringIO("10%\nFinished\n")
            self.returncode = 0
            self.terminated = False

        def wait(self) -> int:
            return self.returncode

        def terminate(self) -> None:
            self.terminated = True
            self.returncode = 143

    def fake_popen(cmd, **kwargs):  # noqa: ANN001 - signature mirrors subprocess
        popen_calls.append(cmd)
        return DummyProcess()

    monkeypatch.setattr(app_flask.subprocess, "Popen", fake_popen)

    payload = json.dumps({"version": "0.1"}).encode("utf-8")
    response = client.post(
        "/upload",
        data={
            "geometry_file": (io.BytesIO(payload), "test.json"),
            "solver": "cg",
            "tol": "1e-6",
            "max_iters": "5",
        },
        content_type="multipart/form-data",
    )

    assert response.status_code == 302

    deadline = time.time() + 1
    while time.time() < deadline and not popen_calls:
        time.sleep(0.01)

    assert popen_calls, "subprocess.Popen should be invoked"

    events = []
    deadline = time.time() + 2
    while time.time() < deadline:
        try:
            item = app_flask.manager.queue.get(timeout=0.1)
            events.append(item)
            if item.get("complete"):
                break
        except queue.Empty:
            pass

    assert any(evt.get("started") for evt in events)
    assert any(evt.get("complete") for evt in events)
    assert not app_flask.manager.is_running


def test_preview_generates_geometry_image(client):
    payload = json.dumps(
        {
            "version": "0.2",
            "domain": {"Lx": 0.2, "Ly": 0.1},
            "sources": [],
        }
    ).encode("utf-8")

    response = client.post(
        "/upload",
        data={
            "geometry_file": (io.BytesIO(payload), "preview.json"),
            "solver": "cg",
            "tol": "1e-6",
            "max_iters": "100",
            "action": "preview",
        },
        content_type="multipart/form-data",
    )

    assert response.status_code == 200
    body = response.data.decode("utf-8")
    assert "geometry-preview-image" in body
    results_dir = Path(app_flask.app.config["RESULTS_FOLDER"])
    generated = list(results_dir.glob("geometry_preview_*.png"))
    assert generated, "Expected a preview image to be created"
    with client.session_transaction() as session:
        project_id = session.get("project_id")
    assert project_id in app_flask.PROJECTS


def test_preview_invalid_json_reports_error(client):
    response = client.post(
        "/upload",
        data={
            "geometry_file": (io.BytesIO(b"not-json"), "invalid.json"),
            "action": "preview",
        },
        content_type="multipart/form-data",
        follow_redirects=True,
    )

    assert response.status_code == 200
    assert "Uploaded JSON is invalid" in response.get_data(as_text=True)


def test_run_after_preview_without_new_upload(monkeypatch, client):
    popen_calls: List[List[str]] = []

    class DummyProcess:
        def __init__(self) -> None:
            self.stdout = io.StringIO("10%\nFinished\n")
            self.returncode = 0
            self.terminated = False

        def wait(self) -> int:
            return self.returncode

        def terminate(self) -> None:
            self.terminated = True

    def fake_popen(cmd, **kwargs):  # noqa: ANN001 - signature mirrors subprocess
        popen_calls.append(cmd)
        return DummyProcess()

    monkeypatch.setattr(app_flask.subprocess, "Popen", fake_popen)

    scenario_payload = json.dumps({"version": "0.2", "domain": {"Lx": 0.1, "Ly": 0.1}, "sources": []}).encode(
        "utf-8"
    )

    preview_response = client.post(
        "/upload",
        data={
            "geometry_file": (io.BytesIO(scenario_payload), "reuse.json"),
            "solver": "cg",
            "tol": "1e-6",
            "max_iters": "50",
            "action": "preview",
        },
        content_type="multipart/form-data",
    )

    assert preview_response.status_code == 200
    with client.session_transaction() as session:
        project_id = session.get("project_id")

    assert project_id in app_flask.PROJECTS
    stored_path = app_flask.PROJECTS[project_id]["scenario_path"]

    response = client.post(
        "/upload",
        data={
            "project_id": project_id,
            "solver": "cg",
            "tol": "1e-6",
            "max_iters": "5",
            "action": "run",
        },
        content_type="multipart/form-data",
    )

    assert response.status_code == 302

    deadline = time.time() + 1
    while time.time() < deadline and not popen_calls:
        time.sleep(0.01)

    assert popen_calls, "subprocess.Popen should run without re-uploading"
    assert str(stored_path) in popen_calls[0]


def test_visualization_route_after_run(monkeypatch, client):
    stdout_queue: "queue.Queue[str | None]" = queue.Queue()

    class DummyStdout:
        def __iter__(self):
            return self

        def __next__(self):
            item = stdout_queue.get()
            if item is None:
                raise StopIteration
            return item

    class DummyProcess:
        def __init__(self) -> None:
            self.stdout = DummyStdout()
            self.returncode = 0

        def wait(self) -> int:
            return self.returncode

        def terminate(self) -> None:
            self.returncode = 143

    def fake_popen(cmd, **kwargs):  # noqa: ANN001 - signature mirrors subprocess
        return DummyProcess()

    monkeypatch.setattr(app_flask.subprocess, "Popen", fake_popen)

    scenario_payload = json.dumps(
        {
            "version": "0.2",
            "domain": {"Lx": 0.1, "Ly": 0.1},
            "sources": [
                {"type": "wire", "x": 0.0, "y": 0.0, "radius": 0.002, "I": 5.0}
            ],
            "outputs": [
                {
                    "type": "field_map",
                    "id": "domain_field",
                    "path": "outputs/test_field.csv",
                }
            ],
        }
    ).encode("utf-8")

    response = client.post(
        "/upload",
        data={
            "geometry_file": (io.BytesIO(scenario_payload), "scenario.json"),
            "solver": "cg",
            "tol": "1e-6",
            "max_iters": "100",
            "outputs": "domain_field",
        },
        content_type="multipart/form-data",
        follow_redirects=False,
    )

    assert response.status_code == 302

    deadline = time.time() + 2
    scenario_path: Path | None = None
    while time.time() < deadline:
        scenario = app_flask.manager._metadata.get("scenario_path")  # type: ignore[attr-defined]
        if isinstance(scenario, Path):
            scenario_path = scenario
            break
        time.sleep(0.01)

    assert scenario_path is not None

    field_path = scenario_path.parent / "outputs" / "test_field.csv"
    field_path.parent.mkdir(parents=True, exist_ok=True)
    field_path.write_text(
        "x,y,Bx,By,Bmag\n"
        "0.0,0.0,0.1,0.0,0.1\n"
        "0.05,0.0,0.1,0.0,0.1\n"
        "0.0,0.05,0.1,0.0,0.1\n"
        "0.05,0.05,0.1,0.0,0.1\n",
        encoding="utf-8",
    )

    stdout_queue.put("10%\n")
    stdout_queue.put("Frame 0: wrote field_map 'domain_field' to \"outputs/test_field.csv\"\n")
    stdout_queue.put("Finished\n")
    stdout_queue.put(None)

    events = []
    deadline = time.time() + 3
    while time.time() < deadline:
        try:
            item = app_flask.manager.queue.get(timeout=0.1)
            events.append(item)
            if item.get("complete"):
                break
        except queue.Empty:
            pass

    assert any("visualization" in evt for evt in events if isinstance(evt, dict))
    last_result = app_flask.manager.get_last_result()
    assert Path(last_result.get("field_map", "")).exists()

    viz_response = client.get(
        "/visualization.png",
        query_string={"vector": "off", "boundaries": "0", "skip": "2"},
    )
    assert viz_response.status_code == 200
    assert viz_response.mimetype == "image/png"
    payload = b"".join(viz_response.response)
    assert payload, "Expected PNG bytes from visualisation route"


def test_project_scenario_download(client):
    payload = json.dumps({"version": "0.2", "domain": {"Lx": 0.1, "Ly": 0.1}, "sources": []}).encode("utf-8")

    client.post(
        "/upload",
        data={
            "geometry_file": (io.BytesIO(payload), "download.json"),
            "action": "preview",
        },
        content_type="multipart/form-data",
    )

    with client.session_transaction() as session:
        project_id = session.get("project_id")

    assert project_id in app_flask.PROJECTS
    scenario_path = app_flask.PROJECTS[project_id]["scenario_path"]
    response = client.get("/project/scenario")
    assert response.status_code == 200
    assert response.data == scenario_path.read_bytes()
    assert "download" in response.headers.get("Content-Disposition", "")


def test_project_reset_clears_state(client):
    payload = json.dumps({"version": "0.2", "domain": {"Lx": 0.1, "Ly": 0.1}, "sources": []}).encode("utf-8")

    client.post(
        "/upload",
        data={
            "geometry_file": (io.BytesIO(payload), "reset.json"),
            "action": "preview",
        },
        content_type="multipart/form-data",
    )

    with client.session_transaction() as session:
        project_id = session.get("project_id")

    scenario_path = app_flask.PROJECTS[project_id]["scenario_path"]
    assert scenario_path.exists()

    response = client.post("/project/reset")
    assert response.status_code == 302
    assert "notice=Project+cleared." in response.headers.get("Location", "")

    with client.session_transaction() as session:
        assert "project_id" not in session

    assert project_id not in app_flask.PROJECTS
    assert not scenario_path.exists()


def test_project_reset_requires_idle(client):
    app_flask.manager._running = True  # type: ignore[attr-defined]
    try:
        response = client.post("/project/reset")
        assert response.status_code == 302
        assert "error=Stop+the+running+simulation" in response.headers.get("Location", "")
    finally:
        app_flask.manager._running = False  # type: ignore[attr-defined]


def test_visualization_route_rejects_invalid_params(visualization_ready, client):
    response = client.get("/visualization.png", query_string={"vector": "bad"})
    assert response.status_code == 400

    response = client.get("/visualization.png", query_string={"skip": "0"})
    assert response.status_code == 400

    response = client.get("/visualization.png", query_string={"scale": "bad"})
    assert response.status_code == 400


def test_visualization_route_accepts_streamlines(visualization_ready, client):
    response = client.get(
        "/visualization.png",
        query_string={
            "vector": "linear",
            "skip": "2",
            "scale": "log",
            "boundaries": "0",
            "stream": "1",
            "log_floor": "1e-9",
            "vector_floor": "1e-9",
        },
    )

    assert response.status_code == 200
    assert response.mimetype == "image/png"
    assert b"".join(response.response)


def test_stop_route_terminates_process(client):
    class DummyProcess:
        def __init__(self) -> None:
            self.terminated = False

        def terminate(self) -> None:
            self.terminated = True

    dummy = DummyProcess()
    app_flask.manager.reset()
    app_flask.manager._process = dummy  # type: ignore[attr-defined]
    app_flask.manager._running = True  # type: ignore[attr-defined]

    response = client.post("/stop")
    assert response.status_code == 202
    assert dummy.terminated

    message = app_flask.manager.queue.get_nowait()
    assert "Termination requested" in message["message"]


def test_progress_stream_emits_payloads(client):
    app_flask.manager.reset()
    app_flask.manager.queue.put({"message": "hello"})
    app_flask.manager.queue.put({"complete": True, "message": "done"})

    response = client.get("/progress")
    body = b"".join(response.response)

    assert b"hello" in body
    assert b"done" in body
