"""Unit tests for the Streamlit GUI helpers."""

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from typing import List

from python.gui import app_streamlit


class DummyUpload:
    def __init__(self, name: str, data: bytes) -> None:
        self.name = name
        self._data = data

    def getvalue(self) -> bytes:
        return self._data


class DummyStdout:
    def __init__(self, lines: List[str]) -> None:
        self._lines = lines
        self._index = 0
        self.closed = False

    def readline(self) -> str:
        if self._index < len(self._lines):
            line = self._lines[self._index]
            self._index += 1
            return line
        return ""

    def close(self) -> None:
        self.closed = True


class DummyProcess:
    def __init__(self, lines: List[str], return_code: int = 0) -> None:
        self.stdout = DummyStdout(lines)
        self._return_code = return_code
        self._terminated = False

    def poll(self) -> int | None:
        if self.stdout._index >= len(self.stdout._lines):
            return self._return_code
        return None

    def wait(self) -> int:
        return self._return_code

    def terminate(self) -> None:
        self._terminated = True

    @property
    def terminated(self) -> bool:
        return self._terminated


class StageScenarioTests(unittest.TestCase):
    def test_stage_default_scenario(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path, messages = app_streamlit.stage_scenario_file(
                uploaded_file=None,
                scratch_dir=Path(tmpdir),
            )
            self.assertIsNotNone(path)
            self.assertTrue(messages)
            assert path is not None
            self.assertTrue(path.exists())
            self.assertEqual(
                path.read_bytes(),
                app_streamlit.DEFAULT_SCENARIO_PATH.read_bytes(),
            )

    def test_stage_json_upload(self) -> None:
        payload = json.dumps({"version": "0.2", "metadata": {}}).encode("utf-8")
        upload = DummyUpload("custom.json", payload)
        with tempfile.TemporaryDirectory() as tmpdir:
            path, messages = app_streamlit.stage_scenario_file(upload, Path(tmpdir))
            self.assertIsNotNone(path)
            self.assertIn("uploaded JSON", " ".join(messages))
            self.assertEqual(path.read_bytes(), payload)

    def test_stage_dxf_upload(self) -> None:
        payload = b"0\nSECTION\nENDSEC\nEOF\n"
        upload = DummyUpload("shape.dxf", payload)
        with tempfile.TemporaryDirectory() as tmpdir:
            path, messages = app_streamlit.stage_scenario_file(upload, Path(tmpdir))
            self.assertIsNone(path)
            self.assertTrue(any("DXF" in message for message in messages))


class DownloadableScenarioTests(unittest.TestCase):
    def test_build_downloadable_scenario_injects_metadata(self) -> None:
        base = json.dumps({"version": "0.2", "metadata": {"foo": "bar"}})
        rendered = app_streamlit.build_downloadable_scenario(base, "cg", 1e-6, 12345)
        payload = json.loads(rendered)
        self.assertIn("ui_defaults", payload["metadata"])
        self.assertEqual(payload["metadata"]["ui_defaults"]["solver"], "cg")
        self.assertEqual(payload["metadata"]["ui_defaults"]["max_iters"], 12345)

    def test_build_downloadable_scenario_handles_invalid_json(self) -> None:
        rendered = app_streamlit.build_downloadable_scenario("not json", "sor", 1e-5, 42)
        payload = json.loads(rendered)
        self.assertEqual(payload["metadata"]["ui_defaults"]["solver"], "sor")


class StreamProcessTests(unittest.TestCase):
    def test_stream_process_output_collects_logs(self) -> None:
        process = DummyProcess(["line 1\n", "line 2\n"])
        logs: List[str] = []
        progress_updates: List[int] = []

        def stop_checker() -> bool:
            return False

        exit_code = app_streamlit.stream_process_output(
            process,
            stop_checker,
            logs.append,
            progress_updates.append,
            idle_sleep=0.0,
        )

        self.assertEqual(exit_code, 0)
        self.assertEqual(logs, ["line 1", "line 2"])
        self.assertTrue(progress_updates)
        self.assertFalse(process.terminated)

    def test_stream_process_output_honours_stop(self) -> None:
        process = DummyProcess(["line 1\n", "line 2\n"])
        logs: List[str] = []
        progress_updates: List[int] = []

        def stop_checker() -> bool:
            return len(progress_updates) >= 1

        app_streamlit.stream_process_output(
            process,
            stop_checker,
            logs.append,
            progress_updates.append,
            idle_sleep=0.0,
        )

        self.assertTrue(process.terminated)
        self.assertTrue(any("Stop requested" in line for line in logs))


class OutputCollectionTests(unittest.TestCase):
    def test_collect_output_paths_from_scenario(self) -> None:
        payload = {
            "outputs": [
                {"path": "outputs/demo.csv"},
                {"path": str(app_streamlit.REPO_ROOT / "outputs" / "absolute.csv")},
                {"path": 123},
            ]
        }
        paths = app_streamlit.collect_output_paths_from_scenario(json.dumps(payload))
        self.assertEqual(len(paths), 2)
        self.assertTrue(all(isinstance(p, Path) for p in paths))
        self.assertTrue(all(p.is_absolute() for p in paths))


class TestRunHook(unittest.TestCase):
    def test_test_run_executes_without_error(self) -> None:
        # The hook prints to stdout; we only need to ensure it completes.
        app_streamlit._test_run()
        # No assertion required—lack of exceptions is success.


if __name__ == "__main__":
    unittest.main()
