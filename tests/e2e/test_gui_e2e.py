import re
import threading
from pathlib import Path

import pytest
from werkzeug.serving import make_server

pytest.importorskip("flask")
playwright_sync = pytest.importorskip("playwright.sync_api")
from playwright.sync_api import Error, expect  # type: ignore

from python.gui import app_flask


@pytest.fixture(scope="session")
def sample_assets(tmp_path_factory):
    assets_dir = tmp_path_factory.mktemp("e2e_assets")

    scenario_src = Path("inputs/induction_gui_demo.json")
    scenario_copy = assets_dir / "induction_gui_demo.json"
    scenario_copy.write_text(scenario_src.read_text(encoding="utf-8"), encoding="utf-8")

    domain_dxf = Path("docs/assets/dxf/induction_demo/domain.dxf")
    rotor_bars_dxf = Path("docs/assets/dxf/induction_demo/rotor_bars.dxf")

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

    yield f"http://127.0.0.1:{server.server_port}"

    server.shutdown()
    thread.join()
    app_flask.manager.reset()
    app_flask.PROJECTS.clear()


@pytest.fixture()
def page(live_server, playwright):
    try:
        browser = playwright.chromium.launch()
    except Exception as exc:  # pragma: no cover - environment guard
        pytest.skip(f"Playwright Chromium launch failed: {exc}")

    context = browser.new_context()
    page = context.new_page()
    page.goto(live_server)
    yield page
    context.close()
    browser.close()


@pytest.fixture(scope="session")
def playwright():
    with playwright_sync.sync_playwright() as p:
        yield p


def test_induction_workflow_end_to_end(page, sample_assets):
    page.set_input_files('[data-testid="dxf-upload-input"]', [sample_assets["domain_dxf"], sample_assets["rotor_bars_dxf"]])
    page.get_by_test_id("dxf-import").click()

    expect(page.get_by_text(sample_assets["domain_dxf"].name)).to_be_visible()
    expect(page.get_by_text(sample_assets["rotor_bars_dxf"].name)).to_be_visible()

    layer_checkbox = page.get_by_test_id(re.compile("layer-select"))
    if layer_checkbox.count() > 0:
        layer_checkbox.first.check()

    update_button = page.get_by_test_id(re.compile("layer-update"))
    if update_button.count() > 0:
        update_button.first.click()

    page.set_input_files('[data-testid="scenario-upload"]', str(sample_assets["scenario"]))
    page.get_by_test_id("run-simulation").click()

    expect(page.get_by_test_id("progress-card")).to_be_visible()
    expect(page.get_by_test_id("log-output")).to_contain_text("Simulation complete.")

    expect(page.get_by_test_id("downloads-card")).to_be_visible()
    expect(page.get_by_test_id("downloads-list")).to_contain_text("Field map")

    expect(page.get_by_test_id("result-image")).to_be_visible()
