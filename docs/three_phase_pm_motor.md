# Three-Phase PM Motor Scenario

This walkthrough builds on the rotating-field stator demo by adding a
surface-mounted permanent-magnet rotor, per-phase RL circuits, and a lightweight
mechanical model. Use it to explore torque production, back-EMF, and circuit
response in a synchronous motor setting. The generator now supports both a
locked-rotor mode (matching earlier releases) and a spin-up mode that lets the
mechanical integrator advance the rotor angle from the electromagnetic torque.

## Quickstart (locked rotor)

```bash
python3 python/gen_three_phase_pm_motor.py --profile ci --mode locked --out inputs/three_phase_pm_motor_ci.json
./build/motor_sim --scenario inputs/three_phase_pm_motor_ci.json --solve --vtk-series outputs/pm_motor_ci.pvd --tol 5e-6 --max-iters 40000
python3 python/animate_three_phase.py --pvd outputs/pm_motor_ci.pvd --scenario inputs/three_phase_pm_motor_ci.json --save pm_motor_demo.mp4 --frame-png pm_motor_demo.png
```

The locked mode emits the same deterministic timeline used in earlier demos: it
drives sinusoidal phase voltages, integrates coil currents via the circuit
solver, and keeps the rotor phased to the stator field via explicit
`rotor_angles` entries.

## Spin-up mode

```bash
python3 python/gen_three_phase_pm_motor.py --profile ci --mode spinup --out inputs/three_phase_pm_motor_spinup_ci.json
./build/motor_sim --scenario inputs/three_phase_pm_motor_spinup_ci.json --solve --vtk-series outputs/pm_motor_spinup_ci.pvd --tol 5e-6 --max-iters 40000
python3 python/check_pm_spinup.py --mechanical outputs/pm_motor_spinup_mechanical.csv --scenario inputs/three_phase_pm_motor_spinup_ci.json
python3 python/animate_three_phase.py --pvd outputs/pm_motor_spinup_ci.pvd --scenario inputs/three_phase_pm_motor_spinup_ci.json --save pm_motor_spinup.mp4 --frame-png pm_motor_spinup.png
```

Spin-up removes the timeline rotor overrides so the RK4 mechanical integrator
advances the rotor based on the torque probe feedback. The generator also emits
`pm_motor_spinup_mechanical.csv`, a rotor angle/speed log that the
`python/check_pm_spinup.py` helper validates for monotonic acceleration.

## Profiles

`python/gen_three_phase_pm_motor.py` exposes the same `ci` and `hires` presets as
the stator generator. Use `--mode locked` (default) for the deterministic rotor
timeline or `--mode spinup` to enable mechanical integration. All profiles accept
`--cycles` and `--frames-per-cycle` overrides when you need shorter or longer
timelines:

- `ci`: 65×65 grid, 12 frames (one electrical cycle)
- `hires`: 401×401 grid, 120 frames per cycle, three electrical cycles

Pass `--profile hires` to emit the larger dataset. All downstream commands stay
the same.

## Scenario highlights

- **Rotor assembly** – A rigid body named `pm_rotor` groups the rotor iron
  polygon and the rectangular magnet block. Locked mode includes
  `"rotor_angles"` timeline entries to phase-lock the magnet to the rotating
  stator field by a configurable load angle, while spin-up mode omits them so
  the mechanical integrator owns the pose.
- **Circuits** – The `stator_three_phase` network models each phase as a series
  `R-L` branch with a driven voltage source. Coil links bind both slot polygons
  for a phase to the shared inductor so the circuit solver injects equal and
  opposite current densities.
- **Mechanical coupling** – The `mechanical` section specifies rotor inertia,
  viscous damping, constant load torque, and the torque probe that feeds the
  integrator. Spin-up mode keeps the timeline free of overrides so the solver
  performs an RK4 update after each frame to evolve angle and speed.
- **Probes** – Back-EMF loops measure the induced voltage inside every positive
  slot, while the torque probe wraps a circular loop just outside the rotor to
  collect Maxwell stress and co-energy samples.

## Outputs

- Locked mode uses the `pm_motor_*.csv/.vti/.pvd` naming scheme. Spin-up mode
  mirrors the same exports with a `pm_motor_spinup_*` prefix and adds
  `pm_motor_spinup_mechanical.csv` capturing time, angle, speed, and torque for
  each rotor sample.
- Torque CSVs include a co-energy column so virtual-work finite differences can
  be computed alongside the Maxwell stress data.

### CI spin-up fixture at a glance

The regression stored at `inputs/tests/pm_motor_spinup_test.json` is generated
verbatim from:

```bash
python3 python/gen_three_phase_pm_motor.py --profile ci --mode spinup \
  --out inputs/tests/pm_motor_spinup_test.json
```

The JSON looks bulky because the generator rasterises each circular boundary
with a few hundred segments to keep the Maxwell-stress torque probe smooth, and
because the timeline spans three electrical cycles (36 frames) so the mechanical
solver sees a full rotation. Key baked-in values for the CI fixture are:

| Quantity | Value |
| --- | --- |
| Grid | 65×65 Cartesian cells over a 0.14 m square |
| Electrical frequency | 60 Hz (12 frames per cycle, 3 cycles) |
| Rotor inertia / damping | 2.5×10⁻³ kg·m², 1.5×10⁻⁴ N·m·s |
| Load torque | 4 N·m opposing rotation |
| Magnet strength | 4×10⁵ A/m surface-mounted block |
| Phase drive | 200 A peak sinusoidal current with 325 V peak phase voltage |

Rather than editing the JSON manually, re-run the generator when you need to
tweak those knobs so the derived polygons stay consistent.

### Faster local experiments

The bundled test purposely remains heavy (~140 s on two vCPUs) so the torque and
mechanical traces match the long-run CI artefacts. For day-to-day iteration you
can emit smaller timelines without touching the committed fixture. A few useful
variants:

```bash
# Half a cycle on the CI grid (6 frames) for quick smoke checks
python3 python/gen_three_phase_pm_motor.py --profile ci --mode spinup \
  --frames-per-cycle 6 --cycles 1 --out scratch/spinup_quick.json

# Coarser mechanical sampling on a reduced grid
python3 python/gen_three_phase_pm_motor.py --profile ci --mode spinup \
  --cycles 1 --frames-per-cycle 8 --out scratch/spinup_57.json

# Re-run only the spin-up regression once you have a new fixture
ctest --test-dir build --output-on-failure -R pm_motor_spinup
```

Keep the official regression JSON in place so GitHub Actions exercises the full
solve, but feel free to adjust the generator overrides locally when validating
mechanical changes or iterating on torque probes.

## Tips

- Run `python/check_three_phase_field.py` against the generated `.pvd` to verify
  the bore angle rotates monotonically and maintains a strong |B| magnitude.
- Use `python/check_pm_spinup.py` on spin-up runs to assert that rotor angle and
  speed increase throughout the timeline.
- The mechanical integrator respects timeline rotor overrides; locked mode
  retains them while spin-up drops them so the solver can demonstrate the coupled
  evolution.
- Pair the torque CSV with the co-energy column to compute finite-difference
  virtual-work torque checks over neighbouring frames.
- Keep generated VTK/MP4 artefacts out of the repository. CI uploads small
  samples for review when the workflow executes the PM motor pipeline.
