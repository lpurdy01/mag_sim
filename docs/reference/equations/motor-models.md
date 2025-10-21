# Motor Models

These canonical relations summarise the quantities referenced across the machine scenario guides. All angles use electrical radians unless stated otherwise.

## Three-phase stator

Balanced phase currents follow

$$
i_a = I_{\text{peak}} \cos(\omega t), \quad
i_b = I_{\text{peak}} \cos\left(\omega t - \tfrac{2\pi}{3}\right), \quad
i_c = I_{\text{peak}} \cos\left(\omega t + \tfrac{2\pi}{3}\right).
$$

The resulting rotating magnetomotive force can be approximated near the bore by

$$
\mathcal{F}(\theta, t) \approx \tfrac{3}{2} N I_{\text{peak}} \cos(\theta - \omega t),
$$

where \(N\) is the number of turns per phase.

## Back-EMF sampling

When sampling induced voltage from a moving conductor the simulation reports

$$
e(t) = -\frac{d\Phi(t)}{dt},
$$

with \(\Phi(t)\) obtained by integrating \(B_n\) over the winding cross-section. See [Back EMF](../../developer-guide/math-and-solver/back-emf.md) for implementation notes.
