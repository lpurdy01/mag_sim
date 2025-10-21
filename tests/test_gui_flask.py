import io
import json
import queue
import time
import sys
from pathlib import Path
from typing import List

import pytest

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
    yield
    app_flask.manager.reset()


@pytest.fixture
def client():
    with app_flask.app.test_client() as test_client:
        yield test_client


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
