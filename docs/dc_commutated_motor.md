# DC motor with commutated armature

This scenario family complements the three-phase demos with a brushed DC
motor that exercises the new commutator primitives. The generator lives in
`python/gen_dc_motor.py` and emits both regression-sized fixtures and larger
profiles suitable for ParaView walkthroughs.

## Generator overview

```
python3 python/gen_dc_motor.py --profile ci --mode spinup --out inputs/dc_motor_spinup_ci.json
python3 python/gen_dc_motor.py --profile hires --mode spinup --out inputs/dc_motor_spinup_hires.json
python3 python/gen_dc_motor.py --profile ci --mode commutator_test --out inputs/tests/dc_commutator_test.json
```

Profiles follow the established CI/hires split. The `ci` configuration keeps the
mesh at 65×65, limits circular tessellation to 24–48 vertices, and rounds all
floating-point coordinates to three decimals so JSON fixtures stay readable. The
`hires` preset scales those counts for offline studies.

`mode` selects the timeline wiring:

* `spinup` – enables the RK4 mechanical simulator, runs a few dozen frames of
  torque-driven acceleration, and emits VTK, torque, and mechanical traces.
  The CI profile now stops after two electrical cycles (24 frames) so the
  regression rotor history remains strictly increasing; the hires preset keeps
  the longer six-cycle spin-up for offline inspection.
* `locked` – leaves the rotor at its initial angle. Handy for debugging static
  field maps or exporting outlines for documentation.
* `commutator_test` – generates a four-frame timeline with scripted rotor angles
  that the regression harness uses to validate the orientation switching logic.
  The armature inductor starts at 8 A so the deposited currents make the sign
  flips obvious.

## Geometry and materials

* Stator yoke: steel ring (`μᵣ ≈ 1800`) with bore radius 20 mm.
* Rotor core: laminated steel disk (`μᵣ ≈ 1200`) centred on the origin and
  attached to the `dc_rotor` mechanical entry.
* Field coils: two slot polygons at ±90° with 220 turns each, a 0.58 fill
  fraction, and default orientation of +1 (north pole) and −1 (south pole).
* Armature coil: two slots at 0° and 180° with 110 turns each, 0.52 fill, and
  base orientations +1 and −1 so the return path flows correctly.

The field and armature ampere-turns were tuned so the stator and rotor fields
live in the same 0.01–0.08 T envelope. The field circuit runs from 18 V through
1.6 Ω/80 mH; the armature is driven by 12 V through 0.75 Ω/32 mH. Both circuits
are defined in the same MNA network (`dc_drive`).

## Commutator schema

Each `coil_link` may now carry a `commutator` block:

```json
{
  "type": "coil_link",
  "inductor": "arm_L",
  "region": "armature_a",
  "turns": 110.0,
  "commutator": {
    "rotor": "dc_rotor",
    "default_orientation": 1.0,
    "segments": [
      {"start_deg": -90.0, "end_deg": 90.0, "orientation": 1.0},
      {"start_deg": 90.0, "end_deg": 270.0, "orientation": -1.0}
    ]
  }
}
```

The solver normalises the rotor angle, picks the segment whose range contains
that angle (wrapping across ±180° as needed), and multiplies the region’s base
orientation by the selected `orientation` factor. The `default_orientation`
value applies if no segment matches (for example when the rotor is exactly on a
boundary). The orientation update happens before each frame solve so both the
current deposition and flux linkage integrals respect the commutator state.

## Regression coverage

* `tests/dc_commutator_test.cpp` loads the `commutator_test` fixture, advances
  through the scripted angles, and asserts that the armature region orientations
  switch signs at ±90° exactly as the JSON specifies. The harness now checks
  that both polarities occur and that at least one sign transition happens, so
  the commutator crossing is exercised explicitly.
* `tests/dc_motor_spinup_test.cpp` reuses the spin-up scenario, performs full CG
  solves with torque probes, and requires ≥8° of electrical rotation and ≥5 rad/s
  of speed gain.

Both binaries are wired into CTest, `scripts/run_ci_checks.sh`, and the GitHub
Actions workflow.

## CI demo

`run_ci_checks.sh` and `.github/workflows/ci.yml` regenerate a fresh
`inputs/dc_motor_spinup_ci.json`, solve it, validate the mechanical trace with
`python/check_pm_spinup.py --rotor dc_rotor`, and archive the VTK series,
mechanical CSV, torque CSV, and outlines. That makes the brushed DC demo visible
alongside the three-phase stator, PM motor, and induction motor artefacts.

## Tuning notes

* Adjust the field/armature balance by changing the turn counts and resistances
  in `DEFAULT_PROFILES`. Because current regions record turn counts, the circuit
  parser enforces consistency across coil links.
* `commutator_test` keeps only four frames so the regression JSON stays compact;
  you can expand the `angles` list in the generator if you want additional
  validation positions.
* The mechanical defaults (0.65×10⁻³ kg·m² inertia, 0.08 N·m load) produce a
  gentle spin-up in CI. Increase the load torque or damping to explore steady
  states, or drop them to showcase faster acceleration in the hires profile.
