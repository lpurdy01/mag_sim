# Force and Torque Probes

Magnetic machines often require estimates of the mechanical force or torque that
the air-gap fields apply to a rotor. `mag_sim` now exposes a **stress tensor
probe** that integrates the Maxwell stress tensor along a closed contour, giving
per-unit-length force and torque values directly from the solved field map.

## Maxwell Stress Tensor refresher

For magnetostatics the Maxwell stress tensor in Cartesian coordinates is

\[
T_{ij} = \frac{1}{\mu_0} \left(B_i B_j - \tfrac{1}{2} \delta_{ij} B^2\right),
\]

with the traction on a surface of unit normal **n** given by

\[
\mathbf{f} = T \cdot \hat{\mathbf{n}}.
\]

Integrating the traction around a closed 2D loop yields the net force per unit
length, while the torque about the out-of-plane axis follows from

\[
\tau_z = \oint (x f_y - y f_x)\,\mathrm{d}s.
\]

The solver evaluates these expressions numerically using midpoint sampling on
polygon edges. Bilinear interpolation recovers **B** at arbitrary contour points
so probes are not restricted to grid-aligned loops.

## JSON schema

Declare probes inside the scenario `outputs` array:

```json
{
  "type": "probe",
  "id": "rotor_mst",
  "probe_type": "force_and_torque",
  "method": "stress_tensor",
  "loop": {
    "type": "polygon",
    "vertices": [
      [-0.03, -0.02],
      [0.03, -0.02],
      [0.03, 0.02],
      [-0.03, 0.02]
    ]
  },
  "path": "outputs/rotor_mst.csv"
}
```

* `probe_type` chooses which quantities are of interest: `"force"`,
  `"torque"`, or `"force_and_torque"` (all three values are written either way).
* `method` currently supports only `"stress_tensor"`.
* `loop` supplies at least three vertices describing the closed contour. Provide
  either an array of `[x, y]` pairs or an object with `type="polygon"` and a
  `vertices` array as shown above.

The solver writes a CSV with a single row:

```text
Fx,Fy,Tz,CoEnergy
-1.234567890123e+02,5.678901234567e+01,-2.468013579240e-03,1.234500000000e-01
```

Values represent force (newtons per metre) and torque (newton-metres per metre)
about the global origin in SI units. The optional `CoEnergy` column captures the
magnetic co-energy integral for the entire slice, enabling finite-difference
virtual-work checks without recomputing the field.

When timelines are active the solver still emits per-frame CSVs with the
`_frame_###` suffix for detailed inspection, but it also aggregates the samples
into the requested base path (for example `outputs/dc_motor_torque.csv`). The
timeline CSV adds `time_s` and `frame_index` columns ahead of the stress tensor
values so downstream scripts can correlate torque with the simulation clock.

## Usage tips

* Keep probe contours well inside the simulation domain so bilinear
  interpolation never samples outside the grid.
* Dense grids yield smoother estimates. For coarse meshes consider slightly
  inflating the loop to avoid sampling immediately adjacent to discretisation
  artefacts.
* Combine with timeline frames to capture torque ripple across electrical angles
  or to cross-check against virtual-work calculations using the magnetic
  co-energy helper described below.

## Virtual-work cross-check

The solver now exposes `motorsim::compute_magnetic_coenergy`, which evaluates
the magnetic co-energy (including permanent-magnet contributions),

\[
W_m = \tfrac{1}{2}\int_{A} \mathbf{B} \cdot (\mathbf{H} + \mathbf{M}) \, \mathrm{d}A,
\]

over the 2D slice. When magnetisation is zero this reduces to the familiar
\(\tfrac{1}{2}\int \mathbf{B}\cdot\mathbf{H}\,\mathrm{d}A\). Combined with timeline frames at neighbouring rotor
angles the virtual-work estimate follows directly from a finite difference,

\[
\tau_z \approx \frac{W_m(\theta + \Delta\theta) - W_m(\theta - \Delta\theta)}{2\,\Delta\theta}.
\]

`tests/torque_validation_test.cpp` exercises the full pipeline by solving the
rotor dipole scenario at ±5° offsets, evaluating both the Maxwell stress torque
and the co-energy difference, and enforcing a ≤10 % agreement on the CI grid.
Any probe requesting torque automatically triggers `computeH()` and records the
co-energy alongside the stress-tensor integral, making the diagnostic available
for future report/CSV exports and the three-phase PM motor walkthrough.

The `tests/probe_output_test.cpp` fixture exercises the ingestion and evaluation
path with a synthetic field that generates a known downward force, providing a
regression guard for the new feature.

[Open in GUI](../developer-guide/dev-environment.md){ .md-button }
