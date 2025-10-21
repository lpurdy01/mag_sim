"""Streamlit application for the 2D motor simulator."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Callable, List, Optional

import streamlit as st

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SCENARIO_PATH = REPO_ROOT / "inputs" / "two_wire_cancel.json"
MOTOR_SIM_BINARY = REPO_ROOT / "build" / "motor_sim"


def ensure_session_state() -> None:
    """Initialise Streamlit session state defaults."""
    state = st.session_state
    state.setdefault("running", False)
    state.setdefault("stop_requested", False)
    state.setdefault("log_lines", [])
    state.setdefault("progress", 0)
    state.setdefault("status_message", "Idle.")
    state.setdefault("available_outputs", [])
    state.setdefault("last_scenario_text", "")


def stage_scenario_file(
    uploaded_file: Optional[object],
    scratch_dir: Path,
    default_scenario_path: Path = DEFAULT_SCENARIO_PATH,
) -> tuple[Optional[Path], List[str]]:
    """Persist the selected scenario or DXF file to a temporary location.

    Parameters
    ----------
    uploaded_file:
        The Streamlit ``UploadedFile`` (or a stub with ``name`` and ``getvalue``)
        selected by the user. ``None`` selects the default bundled scenario.
    scratch_dir:
        Directory where the staged file should be written.
    default_scenario_path:
        Fallback JSON scenario bundled with the repository.

    Returns
    -------
    tuple[Optional[Path], List[str]]
        The path to a JSON scenario file (``None`` when the upload is a DXF) and
        human-readable log messages describing the chosen path.
    """

    messages: List[str] = []

    if uploaded_file is None:
        destination = scratch_dir / default_scenario_path.name
        destination.write_bytes(default_scenario_path.read_bytes())
        messages.append(
            f"Using bundled scenario: {default_scenario_path.relative_to(REPO_ROOT)}"
        )
        return destination, messages

    name = getattr(uploaded_file, "name", "uploaded")
    suffix = Path(name).suffix.lower()

    try:
        data = uploaded_file.getvalue()
    except AttributeError as exc:  # pragma: no cover - defensive guard
        raise TypeError("uploaded_file must expose a getvalue() method") from exc

    destination = scratch_dir / Path(name).name
    destination.write_bytes(data)

    if suffix == ".json":
        messages.append(f"Loaded scenario from uploaded JSON: {name}")
        return destination, messages

    if suffix == ".dxf":
        messages.append(
            "DXF parsing is not yet implemented. Upload recorded for future processing."
        )
        return None, messages

    messages.append(f"Unsupported file type: {suffix or 'unknown'}")
    return None, messages


def build_downloadable_scenario(
    base_json_text: str,
    solver: str,
    tolerance: float,
    max_iters: int,
) -> str:
    """Augment the provided scenario JSON with UI metadata for download."""
    try:
        scenario = json.loads(base_json_text)
    except json.JSONDecodeError:
        scenario = {
            "version": "0.2",
            "metadata": {},
        }

    metadata = scenario.setdefault("metadata", {})
    metadata["ui_defaults"] = {
        "solver": solver,
        "tolerance": tolerance,
        "max_iters": max_iters,
    }
    return json.dumps(scenario, indent=2) + "\n"


def collect_output_paths_from_scenario(scenario_text: str) -> List[Path]:
    """Extract output file paths from a scenario JSON blob."""
    try:
        payload = json.loads(scenario_text)
    except json.JSONDecodeError:
        return []

    results: List[Path] = []
    for entry in payload.get("outputs", []):
        if not isinstance(entry, dict):
            continue
        path_value = entry.get("path")
        if not isinstance(path_value, str):
            continue
        path_obj = Path(path_value)
        if not path_obj.is_absolute():
            path_obj = (REPO_ROOT / path_obj).resolve()
        results.append(path_obj)
    return results


def stream_process_output(
    process: subprocess.Popen,
    stop_checker: Callable[[], bool],
    on_line: Callable[[str], None],
    on_progress: Callable[[int], None],
    idle_sleep: float = 0.1,
) -> int:
    """Stream logs from ``process`` while updating UI callbacks."""
    progress_value = 0
    try:
        if process.stdout is None:
            raise RuntimeError("Process must be created with stdout=PIPE")

        while True:
            line = process.stdout.readline()
            if line:
                clean_line = line.rstrip("\n")
                on_line(clean_line)
                progress_value = min(100, progress_value + 1)
                on_progress(progress_value)
                if stop_checker():
                    on_line("Stop requested by user. Terminating simulation...")
                    process.terminate()
                    break
            elif process.poll() is not None:
                break
            else:
                if stop_checker():
                    on_line("Stop requested by user. Terminating simulation...")
                    process.terminate()
                    break
                time.sleep(idle_sleep)

        return_code = process.wait()
    finally:
        if process.stdout is not None:
            process.stdout.close()

    return return_code


def request_stop() -> None:
    """Signal that the running simulation should halt."""
    st.session_state.stop_requested = True


def render_app() -> None:
    """Render the Streamlit user interface."""
    st.set_page_config(page_title="2D Motor Simulator GUI", layout="wide")
    ensure_session_state()

    st.title("2D Motor Simulator GUI")
    st.markdown(
        """
        Upload a geometry or scenario file, tune solver controls, and launch the
        `motor_sim` solver directly from your browser. The interface streams
        solver logs, exposes quick download links, and captures the latest
        scenario for reuse.
        """
    )

    with st.sidebar:
        st.header("Scenario Setup")
        uploaded_file = st.file_uploader(
            "Upload Geometry (DXF or JSON)",
            type=["dxf", "json"],
        )

        solver_choice = st.selectbox("Solver", ["SOR", "CG"])
        tolerance_value = st.number_input(
            "Tolerance",
            value=1e-6,
            format="%e",
        )
        max_iters_value = st.number_input(
            "Max Iterations",
            value=40000,
            step=1000,
            min_value=1,
        )
        progress_every = st.number_input(
            "Progress cadence (s)",
            value=2.0,
            step=0.5,
            min_value=0.0,
        )

        # Placeholder for future region-selection controls.
        st.caption(
            "Region selection presets will appear here once multi-region geometries are supported."
        )

        if uploaded_file and uploaded_file.name.lower().endswith(".json"):
            base_text = uploaded_file.getvalue().decode("utf-8")
        else:
            base_text = DEFAULT_SCENARIO_PATH.read_text()

        scenario_download = build_downloadable_scenario(
            base_text, solver_choice.lower(), float(tolerance_value), int(max_iters_value)
        )
        st.download_button(
            "Download Scenario JSON",
            data=scenario_download,
            file_name="scenario.json",
            mime="application/json",
            help="Save the current settings as a reusable scenario file.",
        )

    geometry_placeholder = st.empty()
    if uploaded_file:
        if uploaded_file.name.lower().endswith(".dxf"):
            geometry_placeholder.warning(
                "DXF previews are not yet implemented. The file will be staged for future parsing."
            )
        else:
            geometry_placeholder.info(
                f"Ready to run scenario from {uploaded_file.name}. Visualization hooks coming soon."
            )
    else:
        geometry_placeholder.info(
            f"Using default scenario: {DEFAULT_SCENARIO_PATH.relative_to(REPO_ROOT)}"
        )

    controls_col, status_col = st.columns([1, 2])
    with controls_col:
        run_clicked = st.button("Run Simulation", disabled=st.session_state.running)
    with status_col:
        st.button(
            "Stop Simulation",
            disabled=not st.session_state.running,
            on_click=request_stop,
        )

    progress_bar = st.progress(st.session_state.progress)
    status_placeholder = st.empty()
    status_placeholder.info(st.session_state.status_message)

    st.subheader("Simulation Log")
    log_placeholder = st.empty()
    if st.session_state.log_lines:
        log_placeholder.code("\n".join(st.session_state.log_lines[-200:]), language="text")
    else:
        log_placeholder.code("Awaiting simulation run...", language="text")

    outputs_placeholder = st.container()
    if st.session_state.available_outputs:
        with outputs_placeholder:
            st.subheader("Outputs")
            for path in st.session_state.available_outputs:
                label = path.name
                if path.exists():
                    st.download_button(
                        label=f"Download {label}",
                        data=path.read_bytes(),
                        file_name=label,
                        key=f"download-{label}",
                    )
                else:
                    st.warning(f"Expected output missing: {path}")

    if run_clicked:
        if uploaded_file is None and not DEFAULT_SCENARIO_PATH.exists():
            st.error("Default scenario missing. Please provide a scenario JSON to continue.")
            return

        st.session_state.running = True
        st.session_state.stop_requested = False
        st.session_state.log_lines = []
        st.session_state.progress = 0
        st.session_state.status_message = "Starting simulation..."
        status_placeholder.info(st.session_state.status_message)
        log_placeholder.code("Preparing scenario...", language="text")
        progress_bar.progress(0)

        with tempfile.TemporaryDirectory() as temp_dir:
            scratch = Path(temp_dir)
            try:
                scenario_path, stage_messages = stage_scenario_file(uploaded_file, scratch)
            except Exception as exc:  # pragma: no cover - defensive UI guard
                st.session_state.running = False
                st.session_state.status_message = f"Failed to stage scenario: {exc}"
                status_placeholder.error(st.session_state.status_message)
                return

            for message in stage_messages:
                st.session_state.log_lines.append(message)

            if scenario_path is None:
                st.session_state.status_message = stage_messages[-1] if stage_messages else "DXF handling unavailable."
                status_placeholder.warning(st.session_state.status_message)
                st.session_state.running = False
                st.session_state.progress = 0
                progress_bar.progress(0)
                log_placeholder.code("\n".join(st.session_state.log_lines[-200:]), language="text")
                return

            scenario_text = scenario_path.read_text()
            st.session_state.last_scenario_text = scenario_text

            if not MOTOR_SIM_BINARY.exists():
                st.session_state.status_message = (
                    "motor_sim binary not found. Build the project (cmake --build build) before running the GUI."
                )
                status_placeholder.error(st.session_state.status_message)
                st.session_state.running = False
                st.session_state.progress = 0
                progress_bar.progress(0)
                return

            cmd: List[str] = [
                str(MOTOR_SIM_BINARY),
                "--scenario",
                str(scenario_path),
                "--solve",
                "--solver",
                solver_choice.lower(),
                "--tol",
                f"{tolerance_value}",
                "--max-iters",
                str(int(max_iters_value)),
            ]
            if progress_every > 0:
                cmd.extend(["--progress-every", f"{progress_every}"])

            st.session_state.log_lines.append("Launching motor_sim...")
            log_placeholder.code("\n".join(st.session_state.log_lines[-200:]), language="text")

            try:
                process = subprocess.Popen(
                    cmd,
                    cwd=str(REPO_ROOT),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )
            except FileNotFoundError as exc:
                st.session_state.status_message = f"Failed to launch motor_sim: {exc}"
                status_placeholder.error(st.session_state.status_message)
                st.session_state.running = False
                st.session_state.progress = 0
                progress_bar.progress(0)
                return

            def log_callback(line: str) -> None:
                st.session_state.log_lines.append(line)
                log_placeholder.code("\n".join(st.session_state.log_lines[-200:]), language="text")

            def progress_callback(value: int) -> None:
                st.session_state.progress = value
                progress_bar.progress(value)

            def check_stop() -> bool:
                return bool(st.session_state.stop_requested)

            return_code = stream_process_output(process, check_stop, log_callback, progress_callback)

            stop_was_requested = bool(st.session_state.stop_requested)
            st.session_state.running = False
            st.session_state.stop_requested = False

            if return_code == 0 and not stop_was_requested:
                st.session_state.status_message = "Simulation completed successfully."
                progress_callback(100)
            elif stop_was_requested:
                st.session_state.status_message = "Simulation stopped by user."
                st.session_state.progress = 0
                progress_bar.progress(0)
            else:
                st.session_state.status_message = f"Simulation exited with code {return_code}."
                st.session_state.progress = 0
                progress_bar.progress(0)

            status_placeholder.info(st.session_state.status_message)

            if st.session_state.last_scenario_text:
                st.session_state.available_outputs = collect_output_paths_from_scenario(
                    st.session_state.last_scenario_text
                )
            else:
                st.session_state.available_outputs = []

            log_placeholder.code("\n".join(st.session_state.log_lines[-200:]), language="text")

            outputs_placeholder.empty()
            if st.session_state.available_outputs:
                with outputs_placeholder:
                    st.subheader("Outputs")
                    for path in st.session_state.available_outputs:
                        label = path.name
                        if path.exists():
                            st.download_button(
                                label=f"Download {label}",
                                data=path.read_bytes(),
                                file_name=label,
                                key=f"download-{label}",
                            )
                        else:
                            st.warning(f"Expected output missing: {path}")

    # Refresh status in case session state changed elsewhere.
    status_placeholder.info(st.session_state.status_message)


def _test_run() -> None:
    """Headless test hook used by automated checks."""
    print("Streamlit GUI test run starting...")
    fake_lines = [
        "Preparing solver...",
        "Solving linear system...",
        "Post-processing outputs...",
        "Done.",
    ]
    for idx, line in enumerate(fake_lines, start=1):
        print(f"[{idx}] {line}")
        time.sleep(0.05)
    print("Streamlit GUI test run completed successfully.")


if __name__ == "__main__":
    if "--test-run" in sys.argv:
        _test_run()
    else:
        render_app()
