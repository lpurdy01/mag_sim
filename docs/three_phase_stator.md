# Three-Phase Stator Scenario

This scenario bundles a six-slot, two-pole three-phase stator that drives a
balanced rotating field against a permanent-magnet rotor. The Python generator
produces both the tiny CI-sized case and a high-resolution configuration by
tweaking a single profile flag.

## Quickstart

```bash
python3 python/gen_three_phase_stator.py --profile ci --out inputs/three_phase_stator_ci.json
./build/motor_sim --scenario inputs/three_phase_stator_ci.json --solve --vtk-series outputs/three_phase_ci.pvd --tol 5e-6 --max-iters 40000
python3 python/check_three_phase_field.py --pvd outputs/three_phase_ci.pvd --scenario inputs/three_phase_stator_ci.json
python3 python/animate_three_phase.py --pvd outputs/three_phase_ci.pvd --scenario inputs/three_phase_stator_ci.json --save three_phase_demo.mp4 --frame-png three_phase_demo.png
```

The generated scenario exports cell-centred VTK frames, a bore-average CSV,
stress-tensor torque CSVs (one per frame with co-energy samples), and
back-EMF measurements for each phase. Stage 2 adds a `circuits` block that
models the three stator phases as simple RL branches (`R≈0.4 Ω`, `L≈12 mH`)
tied together at a neutral node and a `mechanical` section that tracks the
surface-mounted rotor as a rigid body with inertia, damping, and a constant
load torque. Timeline entries now drive the simulation with sinusoidal phase
voltages via `"voltage_sources"`, after which the solver integrates coil
currents and advances the rotor speed/angle before each magnetostatic frame.

## Scaling up

The generator exposes two profiles:

- `ci`: 65×65 grid, 12 frames (one electrical cycle)
- `hires`: 401×401 grid, 120 frames per cycle, three electrical cycles

Switch profiles via `--profile hires` to emit the high-resolution configuration.
All other pipeline steps remain unchanged—simply re-run `motor_sim` and the
animation command on the new JSON.

## Outputs

- `outputs/three_phase_frame_###.vti`: per-frame cell-centred B/H fields.
- `outputs/three_phase_ci.pvd`: ParaView time-series index (generated via
  `--vtk-series`).
- `outputs/three_phase_outlines.vtp`: geometry polylines for overlaying slot
  and stator boundaries.
- `outputs/bore_angle.csv`: bore-average B components, magnitudes, and angles.
- `outputs/rotor_torque_frame_###.csv`: Maxwell-stress torque with accompanying
  co-energy samples for virtual-work checks.
- `outputs/phase_a_emf.csv`, `outputs/phase_b_emf.csv`, `outputs/phase_c_emf.csv`:
  back-EMF series derived from the positive slot of each phase.

## ParaView tips

1. Open the `.pvd` series to load the time-resolved field data.
2. Add `three_phase_outlines.vtp` as a separate source and enable it in the
   pipeline to overlay slot and stator geometry.
3. Use the “Glyph” filter on `three_phase_outlines.vtp` for quick directional
   cues, or switch the VTI representation to “Surface LIC” for streamline-like
   visuals.

## Animation

`python/animate_three_phase.py` renders a full-field animation that overlays the
cell-centred |B| map, quiver arrows, bore compass, labelled slot outlines, and
the driving phase currents. The CLI accepts:

- `--pvd`: VTK time-series index produced by `motor_sim --vtk-series`.
- `--scenario`: scenario JSON (required to extract timeline currents and bore
  polygon).
- `--save`: output path (MP4/GIF; binaries are uploaded as CI artefacts rather
  than committed).
- `--fps` and `--width`: tune playback speed and output resolution.
- `--html`: emit an interactive HTML player (uses the same data as the MP4).
- `--frame-png`: write a static render of the first frame (handy for docs or
  quick inspection).
- `--log-scale`: switch the |B| colour map to logarithmic scaling to emphasise
  the field in low-magnitude regions.

The animation is designed for the CI demo case and remains lightweight enough
for larger offline runs.

## Notes

- The stator slots are modelled as uniform current regions that draw their
  per-frame currents from `phase_currents` entries in the timeline.
- The bore-average sanity check (`python/check_three_phase_field.py`) unwraps
  the bore field angle, verifies monotonic rotation, enforces an R² > 0.95 fit
  against a straight line, and guards against magnitude collapse. The CI
  workflow runs it automatically.
- The surface-mounted permanent magnet rotor lags the electrical angle by the
  configured load angle so the stress-tensor torque stays positive. The
  `CoEnergy` column in the torque CSV enables a virtual-work cross-check by
  differencing neighbouring frames.
- Back-EMF probes integrate the magnetic flux magnitude inside the positive slot
  of each phase; the resulting CSV already carries the induced `emf` column, so
  no additional differentiation is required.
- The scenario's `circuits` section exposes per-phase voltage sources, series
  resistances, inductances, and coil links for each slot polygon. Timeline
  entries under `"voltage_sources"` set the instantaneous phase voltages,
  allowing the simulator to integrate coil currents with a fourth-order Runge–
  Kutta step before each magnetostatic solve.
- The `mechanical` section specifies rotor inertia, damping, and load torque and
  points at the torque probe feeding the co-simulation loop. The solver solves
  frames sequentially so it can advance the rotor pose with an RK4 step after
  each magnetostatic solve; parallel frame execution is automatically disabled
  when this block is present.
- Keep binary artefacts (MP4/VTI samples) out of git history. The CI workflow
  uploads a small bundle with the demo VTK frame, bore CSV, and animation.
