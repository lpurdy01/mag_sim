# GUI Module Agent Notes

This package contains the Streamlit-based user interface for the simulator.

- `app_streamlit.py` is the main entry point. Keep interactive logic in small,
  testable helpers – functions at module scope should accept pure-Python data
  so we can exercise them from `tests/test_gui_streamlit.py` without spinning up
  a Streamlit runtime.
- Stick to first-party Streamlit widgets. Avoid shipping compiled/bundled
  JavaScript or large static assets; render plots dynamically instead.
- Update `docs/user-guide/gui_streamlit.md` and developer notes when changing
  user-visible behaviour. Extend automated tests as part of any behavioural
  change.
- Remember that long-running subprocess calls must stream logs back to the UI.
  Use the helpers in this module to keep the interface responsive and honour
  the Stop button semantics.
- Do not commit binary artefacts (screenshots, renders, etc.). Generate them at
  runtime if needed for documentation.
