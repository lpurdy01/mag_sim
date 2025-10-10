# Next-Stage Plan Implementation Status

This document tracks progress against the "Next-Stage Development Plan — Time Series, Forces, and Advanced Probes" roadmap.

## Completed Milestones
- **VTK export and verification**: Implemented via `motorsim::write_vti_field_map` with regression coverage in `tests/output_quantity_test.cpp` and exercised in CI through the rotor ripple timeline artefacts. Exports now bundle combined `B`/`H` vector arrays plus geometry outline polydata companions to streamline ParaView workflows. The `python/verify_vtk.py` utility (documented in `docs/vtk_output.md`) runs during CI to sanity-check the generated `.vti` files.
- **Time-series infrastructure**: Scenario timelines, per-frame solving, and parallel execution are available (see `docs/time_series.md`). Rotors now group geometry, magnets, and sources under named pivots so timeline frames can rotate complete assemblies via the new `rotor_angles` overrides. The CI workflow drives `inputs/tests/rotor_ripple_test.json` with `--parallel-frames` to emit multi-frame datasets and CSV/VTI artefacts.
- **Maxwell stress probes**: Force/torque evaluation and documentation reside in `docs/torque_forces.md`, with validation fixtures in `tests/torque_validation_test.cpp` and rotor ripple torque ripple checks.
- **Back-EMF probes**: Flux integration across timeline frames emits per-interval voltages (`docs/back_emf.md`) and is now represented in the artefact bundle via `outputs/rotor_ripple_emf.csv`.
- **Visualization upgrades**: `python/visualize_scenario_field.py` exposes log-scale colormaps, boundary overlays, streamline controls, analytic overlays, and vector scaling modes (documented in `docs/visualization.md`). CI renders both static validation scenes and the new rotor ripple frames.

## CI artefact coverage
- Field plots: `ci_artifacts/rotor_ripple_frame0.png` and `ci_artifacts/rotor_ripple_frame2.png` illustrate timeline evolution alongside the existing two-wire/interface/iron-ring renders.
- Accuracy reports: `ci_artifacts/test_accuracy_report.txt` aggregates analytic/regression solver comparisons, now including force, torque, and back-EMF validation summaries captured from the dedicated regression binaries; additional torque/back-EMF CSVs are copied directly from `outputs/`.
- Sample VTK: `outputs/rotor_ripple_field_frame_000.vti` (and subsequent frames) are exported on every CI run and verified with `python/verify_vtk.py`. Companion outline polydata and `_labels.csv` mapping files are uploaded for ParaView inspection.

## Remaining work
- The roadmap milestones listed above are complete. Future efforts can focus on extending scenario libraries, refining solver performance, or integrating additional analytical benchmarks (e.g., full stator/rotor sweeps or 3D extrusions) as new research questions arise.
