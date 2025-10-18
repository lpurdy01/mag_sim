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
- **Stage 1 (Three-phase synchronous baseline)**: Restored
  `python/gen_three_phase_stator.py` to emit the original stator-only rotating
  field demo. The generator continues to supply CI and hi-res profiles, and
  `docs/three_phase_stator.md` focuses on the bore-angle regression
  (`python/check_three_phase_field.py`), ParaView `.pvd` series, and animation
  workflow without introducing rotor geometry. Maxwell-stress and back-EMF
  probes remain available for stator studies, with regression coverage in
  `tests/back_emf_probe_test.cpp`.
- **Stage 2.1 (RL circuit co-simulation)**: Scenario JSON now accepts
  lumped-element circuits with per-phase resistors, inductors, voltage sources,
  and coil links. Timeline frames drive the network via
  `"voltage_sources"`, the solver integrates branch currents with RK4 (while
  feeding back `-dλ/dt` from coil flux), and the new
  `python/gen_three_phase_pm_motor.py` emits the default three-phase star
  connection for the permanent-magnet motor demo. `tests/circuit_rk_test.cpp`
  exercises the RK integrator against the analytic RL step response, and
  `docs/three_phase_pm_motor.md` documents the voltage-driven workflow. The
  `current_region` source now records slot turns and copper fill fraction so the
  rasteriser deposits ampere-turns instead of raw current; coil links validate
  the turn count so the circuit and field models stay in sync, and the new
  `tests/current_region_turns_test.cpp` regression integrates the deposited
  density to keep the ampere-turn budget within 5% of the analytic target.
  Follow-up work introduced commutator metadata on `coil_link` entries so coil
  orientation can depend on rotor angle. `CircuitSimulator` now reads per-link
  segment tables, queries the current `dc_rotor` angle from the mechanical or
  timeline state, and applies the requested sign flip before depositing currents
  or integrating flux. The DC motor generator (`python/gen_dc_motor.py`) uses
  this plumbing to emit a brushed armature demo alongside a four-frame
  `commutator_test` fixture. `tests/dc_commutator_test.cpp` confirms the
  orientation switching logic, while `tests/dc_motor_spinup_test.cpp` exercises
  the full EM/circuit/mechanical loop with balanced stator and rotor ampere-turn
  budgets. Documentation lives in `docs/dc_commutated_motor.md`, and CI now runs
  the generator, solver, and mechanical validator to archive the new artefacts.
- **Stage 2.2 (Mechanical coupling)**: Introduced a light-weight rotor
  dynamics module that integrates speed and position via RK4 using the torque
  reported by Maxwell-stress probes. Scenarios can now declare inertial,
  damping, and load torque terms under `"mechanical"`; the PM motor generator
  provides a default configuration with a constant load, and the runtime keeps
  timeline frames sequential to ensure warm-started field solves feed the
  coupled ODE. The integrator automatically stands down when timeline frames
  provide explicit rotor angles so deterministic synchronous demos (e.g. the PM
  motor walkthrough) keep their scripted pose, while scenarios without those
  overrides step the coupled mechanics. `tests/mechanical_spinup_test.cpp`
  verifies the integrator against a constant-torque spin-up. Follow-up fixes
  ensure timeline exporters respect basenames that already end in `_frame`,
  restoring the expected `three_phase_frame_000.vti` artefact that CI archives
  for the stator rotating-field demo, and the new
  `tests/pm_motor_spinup_test.cpp` regression drives the full EM/circuit/mechanical
  loop to assert the rotor angle and speed rise over a short spin-up timeline
  while the CLI exposes a `mechanical_trace` output so CI and users can verify
  the rotor state history via `python/check_pm_spinup.py`. The PM motor guide now
  summarises the CI fixture parameters and documents lighter generator overrides
  for local smoke tests so developers can iterate without modifying the stored
  regression JSON. The generator now carves the magnet out of the rotor iron,
  assigns it μᵣ≈1.05, and trims the magnet strength (1×10⁵ A/m) alongside the
  phase drive (35 A peak warm-start currents, 20 V peak). The rotor bore field
  therefore stays on the same tens-of-millitesla scale as the stator-only demo
  instead of spiking into the 50–200 T range that occurred when the magnet
  inherited the 800× permeability of the surrounding steel. Follow-on cleanup
  reduced the CI tessellation counts (24 stator vertices, 18 for the bore, 12 for
  the rotor loop) so the committed spin-up fixture sits near 1.4k lines, and the
  `python/check_pm_spinup.py` helper now evaluates absolute angle/speed gains to
  accommodate rotors that accelerate in either direction.
- **Stage 3 (Frequency-domain induction path)**: Extended the material schema
  with per-material conductivities (`sigma`) and taught the rasteriser and grid
  container to track `sigma`, complex impressed currents, and an imaginary
  vector potential. A new harmonic solver assembles the coupled real/imaginary
  system and applies CG to the normal equations so frequency-domain eddy
  currents can be simulated without leaving the matrix-free framework. Utility
  routines compute complex \(\mathbf{B}\) and \(\mathbf{H}\) fields, and the
  regression `tests/skin_depth_test.cpp` validates skin-depth decay against the
  analytic \(e^{-x/\delta}\) profile for a half-space conductor, enforcing a
  ≤15% slope error on the bundled scenario.
- **Stage 4 (Transient magnetodynamics foundations)**: Introduced a
  Crank–Nicolson-style transient wrapper that augments the existing matrix-free
  operator with \(\sigma/\Delta t\) mass terms and reuses the CG machinery (and
  preconditioners) to march conductive regions forward in time. Timeline frames
  now opt into transient solves via the `"transient"` block, and the main loop
  detects when scripted rotor angles are absent so the coupled circuit/mechanical
  subsystems advance in lock-step with the EM step. The regression
  `tests/diffusion_test.cpp` drives the new
  `inputs/tests/diffusion_test.json` scenario to check magnetic diffusion into a
  conducting slab against the analytic erfc profile, enforcing a ≤20% envelope
  error while archiving the recovered surface field amplitude. Stage 4.2 now ships
  an induction spin-up path: `python/gen_three_phase_induction_motor.py` emits a
  conductive-bar rotor demo, `tests/induction_spinup_test.cpp` verifies the coupled
  transient/mechanical loop accelerates the rotor while staying below synchronous
  speed, and CI captures the accompanying VTK/mechanical artefacts.
- **Stage 5 (Solver UX/perf prep)**: Added a `--pc {none|jacobi|ssor}` CLI flag
  and matching scenario schema to control CG preconditioners. Jacobi and a
  matrix-free SSOR sweep are exposed as first-class options, documented in the
  solver guide alongside usage tips. The default remains unpreconditioned, but
  CI, the local workflow helper, and the docs now spell out when Jacobi or SSOR
  provide faster convergence on conductive or highly anisotropic cases.

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
- Additional UX polish (e.g., non-uniform grid bands and parallel CG kernels)
  can build on the new preconditioner hooks once the transient pipeline settles.
  Longer transient runs that co-simulate the RL network, slip-tuned voltage
  sources, and coarse-to-fine prolongation for induction demos remain on the
  wishlist alongside the optional mechanical enhancements from Stage 5.
