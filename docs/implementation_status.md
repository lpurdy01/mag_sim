# Next-Stage Plan Implementation Status

This document tracks progress against the "Next-Stage Development Plan — Time Series, Forces, and Advanced Probes" roadmap.

## Stage progress — From Magnetostatics to Basic Motor Simulation

- **Stage 0.1 (SIMD-ready linops)**: Extracted the core linear algebra kernels
  (`applyA`, dot products, axpy/scal, residual norms) into `motorsim::linops`
  with optional vectorisation hints guarded by the new
  `MOTORSIM_ENABLE_SIMD_HINTS` CMake switch (`-DMOTORSIM_SIMD_HINTS`). Solver
  code now routes through these utilities, paving the way for SIMD or threaded
  back-ends without changing behaviour.
- **Stage 0.2 (Live progress polish)**: Expanded the runtime controls with
  `--progress-every`, `--snapshot-every`, and the new `--progress-history`
  output. The CLI now records residual samples to CSV (timeline runs append
  `_frame_###`), enabling CI to archive a long CG history directly from
  `motor_sim`.
- **Stage 1 (Three-phase synchronous baseline)**: The updated
  `python/gen_three_phase_stator.py` now emits a surface-mounted PM rotor,
  torque/back-EMF outputs, and timeline rotor angles in addition to the CI and
  hi-res stator profiles. `docs/three_phase_stator.md` covers the workflow,
  including the rotating-field regression (`python/check_three_phase_field.py`),
  ParaView `.pvd` series, torque CSVs with co-energy samples, and phase EMF
  reports. Maxwell-stress probes store co-energy so virtual-work torque checks
  remain within ≤10 % on the CI case, while the back-EMF infrastructure
  (`docs/back_emf.md`, `tests/back_emf_probe_test.cpp`) still verifies the
  120° phase separation and frequency scaling.
- **Stage 2.1 (RL circuit co-simulation)**: Scenario JSON now accepts
  lumped-element circuits with per-phase resistors, inductors, voltage sources,
  and coil links. Timeline frames drive the network via
  `"voltage_sources"`, the solver integrates branch currents with RK4 (while
  feeding back `-dλ/dt` from coil flux), and the Python generator emits the
  default three-phase star connection. `tests/circuit_rk_test.cpp` exercises the
  RK integrator against the analytic RL step response, and
  `docs/three_phase_stator.md` documents the voltage-driven workflow.
- **Stage 2.2 (Mechanical coupling)**: Introduced a light-weight rotor
  dynamics module that integrates speed and position via RK4 using the torque
  reported by Maxwell-stress probes. Scenarios can now declare inertial,
  damping, and load torque terms under `"mechanical"`, the generator emits a
  default PM rotor configuration with a constant load, and the runtime keeps
  timeline frames sequential to ensure warm-started field solves feed the
  coupled ODE. The integrator automatically stands down when timeline frames
  provide explicit rotor angles so deterministic synchronous demos (e.g. the CI
  rotating-field check) keep their scripted pose, while scenarios without those
  overrides step the coupled mechanics. `tests/mechanical_spinup_test.cpp`
  verifies the integrator against a constant-torque spin-up, and
  `docs/three_phase_stator.md` covers the new voltage + mechanical workflow.

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
- Rotor timeline bundles now rotate a multi-segment rotor and high-µ stator, emit outline polydata that ParaView can load without crashing, and log symmetry metrics (`opposition60_120`, `repeat0_180`) plus relative probe/back-EMF errors in the regression summary to track solver accuracy over time.

## Remaining work
- The roadmap milestones listed above are complete. Future efforts can focus on extending scenario libraries, refining solver performance, or integrating additional analytical benchmarks (e.g., full stator/rotor sweeps or 3D extrusions) as new research questions arise.
