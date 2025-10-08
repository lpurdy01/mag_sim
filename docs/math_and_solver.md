# Magnetostatic Solver Math Notes

This document summarizes the equations, discretizations, and validation strategy used by the minimal 2D magnetostatic solver implemented in this repository. It follows the scalar potential formulation with a single non-zero component \(A_z(x, y)\) and highlights practical implementation details for the structured grid solver.

## 1. Problem statement: 2D magnetostatics

We consider a 2D Cartesian cross-section where impressed currents are directed along the out-of-plane \(z\)-axis. The vector potential has only one non-zero component,
\[
\mathbf{A} = (0, 0, A_z(x, y)).
\]
The magnetic flux density is recovered by taking the curl of \(\mathbf{A}\):
\[
\mathbf{B} = \nabla \times (0, 0, A_z) = \bigl(\partial_y A_z,\; -\partial_x A_z,\; 0\bigr).
\]
Consequently,
\[
B_x = \partial_y A_z, \qquad B_y = -\partial_x A_z.
\]

## 2. Governing PDE

Spatially varying materials are represented through the magnetic permeability \(\mu(x, y) = \mu_0 \mu_r(x, y)\). The scalar potential satisfies the Poisson-like equation
\[
\nabla \cdot \left( \frac{1}{\mu} \nabla A_z \right) = -J_z,
\]
where \(J_z\) is the impressed current density in the \(z\) direction. In the special case of uniform permeability, \(\mu = \text{const}\), the equation reduces to
\[
\nabla^2 A_z = -\mu J_z.
\]

## 3. Discretization on a uniform grid

The solver operates on a structured, uniform grid of size \((n_x, n_y)\) with spacings \(\Delta x, \Delta y\). All quantities are stored at cell centres:

- \(A_z\): vector potential solution.
- \(J_z\): source term.
- \(\nu = 1/\mu\): inverse permeability.

To correctly model discontinuities in \(\mu\) at material interfaces, face-centred inverse permeabilities are built via harmonic averaging. For the east face shared by cells \(P\) and \(E\),
\[
\nu_E = \frac{2\, \nu_P \nu_E}{\nu_P + \nu_E},
\]
and similarly for west (W), north (N), and south (S) faces. This gives the standard five-point stencil for diffusion with variable coefficients:
\[
\begin{aligned}
-R_{i,j} &=
\frac{\nu_E (A_{i+1,j} - A_{i,j}) - \nu_W (A_{i,j} - A_{i-1,j})}{\Delta x^2} \\
&\phantom{=}+ \frac{\nu_N (A_{i,j+1} - A_{i,j}) - \nu_S (A_{i,j} - A_{i,j-1})}{\Delta y^2},
\end{aligned}
\]
which should equal \(J_{i,j}\) when the PDE is satisfied. The residual \(R_{i,j}\) is used to monitor convergence.

## 4. Boundary conditions

The current implementation applies Dirichlet conditions \(A_z = 0\) on the outer rectangular boundary, representing a truncated far field. Neumann conditions (zero normal flux) are deferred for future work but noted as a natural extension.

## 5. Iterative solvers

A Gauss–Seidel method with Successive Over-Relaxation (SOR) is implemented as the default solver. For every interior cell,
\[
A_{i,j}^{(\text{new})} = \frac{\tfrac{\nu_E}{\Delta x^2} A_{i+1,j} + \tfrac{\nu_W}{\Delta x^2} A_{i-1,j} + \tfrac{\nu_N}{\Delta y^2} A_{i,j+1} + \tfrac{\nu_S}{\Delta y^2} A_{i,j-1} + J_{i,j}}{\tfrac{\nu_E + \nu_W}{\Delta x^2} + \tfrac{\nu_N + \nu_S}{\Delta y^2}}.
\]
Relaxation is applied as
\[
A_{i,j} \leftarrow (1 - \omega) A_{i,j} + \omega A_{i,j}^{(\text{new})}, \qquad 1 < \omega < 2.
\]
An (unimplemented) conjugate gradient solver is left as a future improvement for the symmetric positive definite system.

## 6. Convergence criteria

Iterations continue until either the maximum iteration count is reached or the relative residual falls below the user-specified tolerance:
\[
\frac{\sqrt{\sum_{i,j} R_{i,j}^2}}{\sqrt{\sum_{i,j} J_{i,j}^2 + \varepsilon}} < \text{tol},
\]
with \(\varepsilon\) guarding against zero-current scenarios.

## 7. Recovering the magnetic flux density

Given \(A_z\), the magnetic flux density components are approximated via finite differences. In the interior we use central differences:
\[
B_x \approx \frac{A_{i,j+1} - A_{i,j-1}}{2\Delta y}, \qquad
B_y \approx -\frac{A_{i+1,j} - A_{i-1,j}}{2\Delta x}.
\]
One-sided differences are applied on the domain boundary where neighbours are unavailable.

## 8. Validation with an infinite straight wire

To validate the solver, we model an infinite straight wire carrying current \(I\). The analytic magnetic field magnitude is
\[
B(r) = \frac{\mu_0 I}{2 \pi r}.
\]
Numerically, the wire is approximated by a small circular region (radius \(r_c\)) with uniform current density,
\[
J_z = \frac{I}{\pi r_c^2}.
\]
Cells with centre radius \(r \le r_c\) receive this current density, while the rest are zero. Choosing \(r_c\) a few cells wide (e.g., three cell widths) reduces discretisation error from approximating the singular source. After solving for \(A_z\), we compute \(\mathbf{B}\) and sample points along a ring of radius \(r_{\text{sample}}\) to compare the simulated \(\lVert \mathbf{B} \rVert\) against the analytic expression. A relative error below 25% is deemed acceptable for the coarse grid and finite domain used in the automated test.

---

This reference captures the mathematical foundations and numerical choices embedded in the minimal solver. Future extensions (e.g., alternative boundary conditions, conjugate gradient solvers, nonlinear materials) should extend these notes accordingly.

For runtime and scaling heuristics, consult `docs/solver_performance.md`, which summarises
benchmark data from the bundled tooling.

## 9. Scenario Spec v0.1

The solver now ingests simulation descriptions authored as JSON documents. The
schema is intentionally small to keep the C++ dependency surface minimal while
leaving room for future extensions (materials with geometry, time timelines,
etc.). A valid document contains the following top-level members:

| Field       | Type   | Notes |
| ----------- | ------ | ----- |
| `version`   | string | Must be `"0.1"` for the current ingestor. |
| `units`     | string | `"SI"` only; currents in amperes, lengths in metres. |
| `domain`    | object | `{Lx, Ly, nx, ny}` define the rectangular grid centred on the origin; `dx = Lx / (nx-1)` and likewise for `dy`. |
| `materials` | array  | Each entry defines `{name, mu_r}`; v0.1 uses a single uniform material. |
| `regions`   | array  | Uniform background assignments: `{ "type": "uniform", "material": "air" }`. |
| `sources`   | array  | Currently limited to wires: `{ "type": "wire", "x", "y", "radius", "I" }` with cylindrical patches of uniform `J_z`. |
| `outputs`   | array  | Optional list of export requests. v0.1 supports `field_map` and `line_probe` records with stable `id`s, formats (CSV), and target paths. |

The C++ helper `loadScenarioFromJson` performs structural validation and
produces a `ScenarioSpec`. `rasterizeScenarioToGrid` then deposits current
density and background permeability directly onto a `Grid2D`. The default `main`
executable wires these steps together:

```text
JSON spec → ScenarioSpec → rasterise to Grid2D → solve A_z → compute B
```

Two helper layers keep authoring ergonomic:

1. `python/scenario_api.py` exposes dataclasses mirroring the schema and a
   `Scenario.save_json()` convenience. Agents can write scenarios in Python,
   validate them, and emit the JSON artefact in one go.
2. `motor_sim --scenario path/to.json --solve [--list-outputs] [--outputs ids]`
   loads the spec, solves the magnetostatic system, and emits any outputs
   declared in the JSON. Use `--list-outputs` to inspect available IDs and
   `--outputs id1,id2` or `--outputs none` to control which requests are
   fulfilled at runtime. The legacy `--write-midline` flag remains available for
   ad-hoc dumps.

The `sources` list supports multiple wire entries. Each is rasterised as a disk
with constant current density `J_z = I / (π r²)` applied to cells whose centres
fall inside the radius. `tests/two_wire_cancel_test.cpp` integrates these cells
to confirm that the deposited current matches the requested value to within
15%, providing a regression guard on the rasteriser.

Reserved future fields (e.g. a `timeline` array for time-varying studies) can be
introduced without breaking the base schema because the parser ignores unknown
members. Geometry primitives beyond uniform regions should extend the `regions`
array with new `type` variants when the solver grows material heterogeneity.

### 9.1 Output requests

Output definitions live alongside the physical description so scenarios are
self-documenting and reproducible. Two request flavours are implemented in v0.1:

* **Field maps** (`{"type":"field_map", "id":"domain_field", "quantity":"B", "path":"outputs/two_wire_field_map.csv"}`)
  dump the full `B` field over the grid. The CSV contains `x,y,Bx,By,Bmag`, and
  the path defaults to `outputs/<id>.csv` when omitted.
* **Line probes** (`{"type":"line_probe", "axis":"x", "value":0.0, "quantity":"Bmag"}`)
  sample a horizontal or vertical line aligned with the grid. Specify `axis`
  (`"x"` or `"y"`), the coordinate to lock, the field component (`Bx`, `By`, or
  `Bmag`), and an output path. The ingestor validates that the requested line
  lands on an existing grid column/row.

Python authors can build these records via
`scenario_api.FieldMapOutput`/`LineProbeOutput`. Downstream tooling such as
`python/visualize_scenario_field.py` consumes the emitted CSV to produce quick
look plots for scenario debugging.
