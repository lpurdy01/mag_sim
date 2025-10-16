# Three-Phase PM Motor Scenario

This walkthrough builds on the rotating-field stator demo by adding a
surface-mounted permanent-magnet rotor, per-phase RL circuits, and a lightweight
mechanical model. Use it to explore torque production, back-EMF, and circuit
response in a synchronous motor setting.

## Quickstart

```bash
python3 python/gen_three_phase_pm_motor.py --profile ci --out inputs/three_phase_pm_motor_ci.json
./build/motor_sim --scenario inputs/three_phase_pm_motor_ci.json --solve --vtk-series outputs/pm_motor_ci.pvd --tol 5e-6 --max-iters 40000
python3 python/animate_three_phase.py --pvd outputs/pm_motor_ci.pvd --scenario inputs/three_phase_pm_motor_ci.json --save pm_motor_demo.mp4 --frame-png pm_motor_demo.png
```

The generator writes a fully coupled scenario that drives sinusoidal phase
voltages, integrates coil currents via the circuit solver, and advances the
rotor pose through a timeline of prescribed electrical angles. The mechanical
integrator automatically stands down when explicit `rotor_angles` are present,
ensuring deterministic comparisons between runs.

## Profiles

`python/gen_three_phase_pm_motor.py` exposes the same `ci` and `hires` presets as
the stator generator:

- `ci`: 65×65 grid, 12 frames (one electrical cycle)
- `hires`: 401×401 grid, 120 frames per cycle, three electrical cycles

Pass `--profile hires` to emit the larger dataset. All downstream commands stay
the same.

## Scenario highlights

- **Rotor assembly** – A rigid body named `pm_rotor` groups the rotor iron
  polygon and the rectangular magnet block. Timeline entries include
  `"rotor_angles"` to phase-lock the magnet to the rotating stator field by a
  configurable load angle.
- **Circuits** – The `stator_three_phase` network models each phase as a series
  `R-L` branch with a driven voltage source. Coil links bind both slot polygons
  for a phase to the shared inductor so the circuit solver injects equal and
  opposite current densities.
- **Mechanical coupling** – The `mechanical` section specifies rotor inertia,
  viscous damping, constant load torque, and the torque probe that feeds the
  integrator. With timeline overrides removed, the solver performs an RK4 update
  after each frame to evolve angle and speed.
- **Probes** – Back-EMF loops measure the induced voltage inside every positive
  slot, while the torque probe wraps a circular loop just outside the rotor to
  collect Maxwell stress and co-energy samples.

## Outputs

- `outputs/pm_motor_frame_###.vti`: per-frame cell-centred B/H snapshots.
- `outputs/pm_motor_ci.pvd`: ParaView series referencing the frame stack.
- `outputs/pm_motor_outlines.vtp`: slot, stator, and rotor outlines for overlay.
- `outputs/pm_motor_bore.csv`: bore-average field magnitude and angle.
- `outputs/pm_motor_torque.csv`: Maxwell stress torque with matching co-energy
  entries for virtual-work validation.
- `outputs/pm_motor_phase_[abc]_emf.csv`: back-EMF waveforms for each phase.

## Tips

- Run `python/check_three_phase_field.py` against the generated `.pvd` to verify
  the bore angle rotates monotonically and maintains a strong |B| magnitude.
- The mechanical integrator respects timeline rotor overrides; remove the
  `rotor_angles` entries to observe a free spin-up against the configured load
  torque and damping.
- Pair the torque CSV with the co-energy column to compute finite-difference
  virtual-work torque checks over neighbouring frames.
- Keep generated VTK/MP4 artefacts out of the repository. CI uploads small
  samples for review when the workflow executes the PM motor pipeline.
