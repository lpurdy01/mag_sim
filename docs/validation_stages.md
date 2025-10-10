# Analytical Validation Ladder

The analytical validation ladder mirrors the "Next-Stage Development Plan — Time Series, Forces, and Advanced Probes" roadmap. Each rung pairs a solvable scenario with an automated regression test so we can quantify agreement with closed-form physics before escalating to full motor geometries.

## Stage summary

| Stage | Scenario / Test | Analytic reference | Acceptance check |
| ----- | ---------------- | ------------------ | ---------------- |
| 1 | `tests/analytic_wire_test` (`inputs/tests/analytic_wire_test.json`) | Biot–Savart field of an infinite line current | RMS relative error of sampled \|B\| at a 5&nbsp;cm radius ring below 25%. |
| 2 | `tests/magnet_strip_test` (`inputs/tests/magnet_strip_test.json`) | Surface-current model of a uniformly magnetised slab | RMS relative error of \|B\| along the midline < 15% and Hy inside the magnet near zero. |
| 3 | `tests/torque_validation_test` (`inputs/tests/torque_validation_test.json`) | Torque on a dipole in an external field | Maxwell stress tensor torque agrees with dipole and virtual-work estimates within 20% / 25%, respectively. |
| 4 | `tests/back_emf_probe_test` (`inputs/tests/back_emf_probe_test.json`) | Faraday law using discrete flux differences | Polygon and rectangular probes integrate flux exactly for synthetic fields and produce the expected EMF between frames. |
| 5 | `tests/rotor_ripple_test` (`inputs/tests/rotor_ripple_test.json`) | Quasi-static PM rotor in a stator field | Torque sign flips as the rotor sweeps 0°→180° and the simulated peak-to-peak ripple exceeds 0.3&nbsp;N·m·m⁻¹. |

All tests execute via `ctest` (see `CMakeLists.txt`) and run as part of the GitHub Actions workflow. They share compact domains so the Gauss–Seidel solver converges within a few seconds per frame.

## Scenario notes

- **Stage 1:** The analytic wire scenario writes CSV artefacts for regression plots. The `python/visualize_wire.py` helper overlays the Biot–Savart solution and is used in CI artefact generation.
- **Stage 2:** The magnet strip benchmark samples the field along `y=0` and compares against a 4,000-segment surface-current quadrature. It also confirms that the magnet interior is nearly demagnetised (Hy ≈ 0).
- **Stage 3:** The torque validation setup places a rectangular magnet between counter-wound conductors. The magnet experiences a uniform transverse field so the torque can be predicted from `τ = (M · A) × B`. The regression checks both the Maxwell stress integration and a virtual-work finite difference (`ΔW/Δθ`).
- **Stage 4:** The back-EMF test exercises polygonal and rectangular integration regions, verifies frame selection semantics, and validates the EMF series using synthetic flux ramps.
- **Stage 5:** The rotor ripple scenario reuses the Stage&nbsp;3 geometry but sweeps the magnetisation vector through four angles (0°, 60°, 120°, 180°) via timeline frames. The regression ensures torque polarity and amplitude evolve consistently with the expected sinusoid.

## Running the ladder locally

```bash
cmake --build build
ctest --output-on-failure -R "analytic_wire|magnet_strip|torque_validation|back_emf_probe|rotor_ripple"
```

The command above rebuilds the project and runs only the validation ladder tests. Use `ctest --output-on-failure --parallel $(nproc)` for the full suite to keep turnaround tight.

## Extending the ladder

Future validation efforts (e.g., laminated stators, eddy-current suppression) should add a dedicated scenario JSON, an automated regression, and a short documentation update. The ladder table gives a template for summarising accuracy targets and analytic references.
