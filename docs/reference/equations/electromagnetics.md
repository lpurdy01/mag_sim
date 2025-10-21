# Electromagnetics

## Magnetostatic formulation

The magnetostatic solve assumes a single non-zero vector potential component \(A_z(x, y)\). The governing PDE is

$$
\nabla \cdot (\nu \nabla A_z) = -J_z,
$$

with \(\nu = 1/\mu\). For uniform materials the expression reduces to the familiar Poisson form

$$
\nabla^2 A_z = -\mu J_z.
$$

## Magnetic flux density

Recovered field components follow from central differences of \(A_z\):

$$
B_x \approx \frac{A_{i,j+1} - A_{i,j-1}}{2\Delta y}, \qquad
B_y \approx -\frac{A_{i+1,j} - A_{i-1,j}}{2\Delta x}.
$$

The magnetic field intensity combines the solved field with magnetisation contributions:

$$
\mathbf{H} = \nu \mathbf{B} - \frac{1}{\mu_r} \mathbf{M}.
$$

## Magnetisation sources

Permanent magnets appear as an equivalent bound current density

$$
J_{m,z} = \partial_x M_y - \partial_y M_x,
$$

which is added to \(J_z\) during rasterisation.
