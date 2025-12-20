import re
import sys
import threading
from pathlib import Path

import pytest
from werkzeug.serving import make_server

pytest.importorskip("flask")
playwright_sync = pytest.importorskip("playwright.sync_api")
from playwright.sync_api import Error, expect  # type: ignore

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from python.gui import app_flask


@pytest.fixture(scope="session")
def sample_assets(tmp_path_factory):
    assets_dir = tmp_path_factory.mktemp("e2e_assets")

    scenario_src = Path("inputs/induction_gui_demo.json")
    scenario_copy = assets_dir / "induction_gui_demo.json"
    scenario_copy.write_text(scenario_src.read_text(encoding="utf-8"), encoding="utf-8")

    def _write_minimal_dxf(path: Path, layer: str) -> None:
        path.write_text(
            "\n".join(
                [
                    "0",
                    "SECTION",
                    "2",
                    "ENTITIES",
                    "0",
                    "LINE",
                    "8",
                    layer,
                    "10",
                    "0",
                    "20",
                    "0",
                    "11",
                    "1",
                    "21",
                    "1",
                    "0",
                    "ENDSEC",
                    "0",
                    "EOF",
                ]
            ),
            encoding="utf-8",
        )

    domain_dxf = assets_dir / "domain.dxf"
    rotor_bars_dxf = assets_dir / "rotor_bars.dxf"
    _write_minimal_dxf(domain_dxf, "DOMAIN")
    _write_minimal_dxf(rotor_bars_dxf, "ROTOR")

    return {
        "scenario": scenario_copy,
        "domain_dxf": domain_dxf,
        "rotor_bars_dxf": rotor_bars_dxf,
    }


@pytest.fixture()
def live_server(tmp_path_factory, monkeypatch):
    uploads = tmp_path_factory.mktemp("uploads")
    results = tmp_path_factory.mktemp("results")

    app_flask.app.config.update(
        TESTING=True,
        SERVER_NAME=None,
        UPLOAD_FOLDER=str(uploads),
        RESULTS_FOLDER=str(results),
        SECRET_KEY="e2e-secret",
    )
    app_flask.manager.reset()
    app_flask.PROJECTS.clear()

    bootstrap_dir = tmp_path_factory.mktemp("bootstrap")
    bootstrap_dxf = bootstrap_dir / "bootstrap.dxf"
    bootstrap_dxf.write_text(
        "\n".join(
            [
                "0",
                "SECTION",
                "2",
                "ENTITIES",
                "0",
                "LINE",
                "8",
                "BOOT",
                "10",
                "0",
                "20",
                "0",
                "11",
                "1",
                "21",
                "1",
                "0",
                "ENDSEC",
                "0",
                "EOF",
            ]
        ),
        encoding="utf-8",
    )

    test_client = app_flask.app.test_client()
    upload_response = test_client.post(
        "/dxf/upload",
        data={"dxf_files": (open(bootstrap_dxf, "rb"), bootstrap_dxf.name)},
        content_type="multipart/form-data",
    )
    session_cookie = None
    for cookie_header in upload_response.headers.getlist("Set-Cookie"):
        if cookie_header.startswith("session="):
            session_cookie = cookie_header.split(";", 1)[0].split("=", 1)[1]
            break

    if not session_cookie:
        raise RuntimeError("Bootstrap DXF upload did not return a session cookie")

    def _fake_run(self, command):  # type: ignore[no-untyped-def]
        scenario_path = self._metadata.get("scenario_path")
        log_path = self._metadata.get("log_path")
        process_cwd = Path(self._metadata.get("process_cwd", Path.cwd()))

        field_map_path = process_cwd / "outputs" / "induction_motor_field.csv"
        field_map_path.parent.mkdir(parents=True, exist_ok=True)
        field_map_path.write_text(
            "x,y,Bx,By,Bmag\n" "0.0,0.0,0.1,0.0,0.1\n" "0.1,0.0,0.1,0.0,0.1\n",
            encoding="utf-8",
        )

        if log_path:
            Path(log_path).write_text(
                "Frame 0: wrote field_map 'induction_motor_field' to \"outputs/induction_motor_field.csv\"\n"
                "Simulation complete.\n",
                encoding="utf-8",
            )

        self.queue.put({
            "started": True,
            "message": "Simulation launched.",
            "progress": 5,
        })
        self.queue.put({
            "message": "Frame 0: wrote field_map 'induction_motor_field' to \"outputs/induction_motor_field.csv\"",
            "progress": 90,
        })
        self._handle_field_map_event(
            frame=0, field_id="induction_motor_field", path_str="outputs/induction_motor_field.csv"
        )

        with self._lock:
            self._running = False
            self._process = None

        if scenario_path:
            self._emit_completion(success=True, message="Simulation complete.", scenario_path=scenario_path)

    monkeypatch.setattr(app_flask.SimulationManager, "_run_process", _fake_run, raising=False)

    server = make_server("127.0.0.1", 0, app_flask.app)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()

    yield f"http://127.0.0.1:{server.server_port}", session_cookie

    server.shutdown()
    thread.join()
    app_flask.manager.reset()
    app_flask.PROJECTS.clear()


@pytest.fixture()
def page(live_server, playwright):
    server_url, session_cookie = live_server
    try:
        browser = playwright.chromium.launch()
    except Exception as exc:  # pragma: no cover - environment guard
        pytest.skip(f"Playwright Chromium launch failed: {exc}")

    context = browser.new_context()
    context.add_cookies(
        [
            {
                "name": "session",
                "value": session_cookie,
                "url": server_url,
            }
        ]
    )
    page = context.new_page()
    page.goto(server_url)
    yield page
    context.close()
    browser.close()


@pytest.fixture(scope="session")
def playwright():
    with playwright_sync.sync_playwright() as p:
        yield p


def test_induction_workflow_end_to_end(page, sample_assets):
    expect(page.get_by_test_id("dxf-upload-input")).to_be_attached()
    expect(page.get_by_test_id("run-simulation")).to_be_attached()
