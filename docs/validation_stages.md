# Analytical Validation Ladder

The analytical validation ladder mirrors the "Next-Stage Development Plan — Time Series, Forces, and Advanced Probes" roadmap. Each rung pairs a solvable scenario with an automated regression test so we can quantify agreement with closed-form physics before escalating to full motor geometries.

## Stage summary

| Stage | Scenario / Test | Analytic reference | Acceptance check |
| ----- | ---------------- | ------------------ | ---------------- |
| 1 | `tests/analytic_wire_test` (`inputs/tests/analytic_wire_test.json`) | Biot–Savart field of an infinite line current | RMS relative error of sampled \|B\| at a 5&nbsp;cm radius ring below 25%. |
| 2 | `tests/magnet_strip_test` (`inputs/tests/magnet_strip_test.json`) | Surface-current model of a uniformly magnetised slab | RMS relative error of \|B\| along the midline < 15% and Hy inside the magnet near zero. |
| 2a | `tests/current_region_turns_test` (`inputs/tests/current_region_turns_test.json`) | Ampere-turn conservation for polygonal slot sources | Integrated current density matches the requested ampere-turns within 5%. |
| 3 | `tests/torque_validation_test` (`inputs/tests/torque_validation_test.json`) | Torque on a dipole in an external field | Maxwell stress tensor torque agrees with dipole and virtual-work estimates within 20% / 25%, respectively. |
| 4 | `tests/back_emf_probe_test` (`inputs/tests/back_emf_probe_test.json`) | Faraday law using discrete flux differences | Polygon and rectangular probes integrate flux exactly for synthetic fields and produce the expected EMF between frames. |
| 5 | `tests/rotor_ripple_test` (`inputs/tests/rotor_ripple_test.json`) | Quasi-static PM rotor in a stator field | Torque sign flips as the rotor sweeps 0°→180° and the simulated peak-to-peak ripple exceeds 0.3&nbsp;N·m·m⁻¹. |
| 6 | `tests/skin_depth_test` (`inputs/tests/skin_depth_test.json`) | Classical skin-depth decay in a conducting half-space | Linear-fit slope of \(|\mathbf{B}|\) vs depth matches \(-1/\delta\) within 15%. |
| 7 | `tests/diffusion_test` (`inputs/tests/diffusion_test.json`) | Semi-infinite magnetic diffusion following the error-function solution | Max relative error of sampled \(B_x\) vs \(\mathrm{erfc}\) profile stays below 20%. |
| 8 | `tests/induction_spinup_test` (`inputs/tests/induction_spinup_test.json`) | Slip-limited acceleration of a conductive rotor | Rotor angle/speed rise > 5 deg / 5 rad·s⁻¹, final speed < 80% synchronous, slip within (0.05, 0.95). |

All tests execute via `ctest` (see `CMakeLists.txt`) and run as part of the GitHub Actions workflow. They share compact domains so the Gauss–Seidel solver converges within a few seconds per frame.

## Scenario notes

- **Stage 1:** The analytic wire scenario writes CSV artefacts for regression plots. The `python/visualize_wire.py` helper overlays the Biot–Savart solution and is used in CI artefact generation.
- **Stage 2:** The magnet strip benchmark samples the field along `y=0` and compares against a 4,000-segment surface-current quadrature. It also confirms that the magnet interior is nearly demagnetised (Hy ≈ 0).
- **Stage 2a:** The current-region turns test integrates the rasterised `J_z` over a coil slot polygon to ensure the deposited ampere-turns stay within 5% of the analytic `orientation × I × turns × fill_fraction` target.
- **Stage 3:** The torque validation setup places a rectangular magnet between counter-wound conductors. The magnet experiences a uniform transverse field so the torque can be predicted from `τ = (M · A) × B`. The regression checks both the Maxwell stress integration and a virtual-work finite difference (`ΔW/Δθ`).
- **Stage 4:** The back-EMF test exercises polygonal and rectangular integration regions, verifies frame selection semantics, and validates the EMF series using synthetic flux ramps.
- **Stage 5:** The rotor ripple scenario reuses the Stage&nbsp;3 geometry but sweeps the magnetisation vector through four angles (0°, 60°, 120°, 180°) via timeline frames. The regression ensures torque polarity and amplitude evolve consistently with the expected sinusoid.
- **Stage 6:** The skin-depth fixture drives a uniform conductor slab with a harmonic current sheet. Sampling \(|\mathbf{B}|\) along the conductor normal and fitting a log-linear slope recovers the expected \(-1/\delta\) decay constant.
- **Stage 7:** The diffusion regression excites a semi-infinite conductor with a step current and marches the transient solve forward using the \(\sigma/\Delta t\) mass term. Averaging \(B_x\) across the slab height and comparing against the analytic \(\mathrm{erfc}\) profile enforces the ≤20% envelope error target.
- **Mechanical spin-up:** `tests/pm_motor_spinup_test.cpp` loads the coupled PM motor spin-up scenario, solving the EM field,
  RL circuits, and RK4 mechanical loop frame-by-frame. The regression asserts that rotor angle and speed rise over the timeline,
  mirroring the CI spin-up demo that writes `pm_motor_spinup_mechanical.csv` for further inspection.
- **Induction spin-up:** `tests/induction_spinup_test.cpp` exercises the transient Crank–Nicolson solver and the mechanical
  integrator on the conductive-bar rotor. It reuses the torque probe to advance the rotor, checks that acceleration is positive,
  and ensures the final speed remains below synchronous so slip stays realistic. CI archives the accompanying mechanical trace,
  VTK series, and outline polydata.

## Running the ladder locally

```bash
cmake --build build
ctest --output-on-failure -R "analytic_wire|magnet_strip|torque_validation|back_emf_probe|rotor_ripple"
```

The command above rebuilds the project and runs only the validation ladder tests. Use `ctest --output-on-failure --parallel $(nproc)` for the full suite to keep turnaround tight.

## Extending the ladder

Future validation efforts (e.g., laminated stators, eddy-current suppression) should add a dedicated scenario JSON, an automated regression, and a short documentation update. The ladder table gives a template for summarising accuracy targets and analytic references.
