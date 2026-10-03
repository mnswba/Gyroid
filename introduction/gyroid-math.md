# Gyroid Surface: Mathematical Definitions

The gyroid is a triply periodic minimal surface (TPMS) discovered by Alan Schoen in 1970.
It is infinitely connected, contains no straight lines, and has no planes of mirror symmetry.

## 1. Implicit (nodal) approximation

The most common definition is the zero level set of a trigonometric function:

$$
F(x,y,z) = \sin x \cos y + \sin y \cos z + \sin z \cos x = 0
$$

Scaled to a unit cell of edge length $L$ (with $k = 2\pi/L$):

$$
F(x,y,z) = \sin(kx)\cos(ky) + \sin(ky)\cos(kz) + \sin(kz)\cos(kx) = 0
$$

- The surface is periodic with period $2\pi$ (or $L$) along $x$, $y$ and $z$.
- $F=0$ divides space into two congruent, interpenetrating labyrinths ($F>0$ and $F<0$).
- **Caveat:** this is a first-order Fourier approximation. The zero set is very close to,
  but not exactly, Schoen's minimal surface (its mean curvature is only approximately zero).

## 2. Level-set generalization (sheet / network solids)

$$
F(x,y,z) = t
$$

- $t = 0$: balanced surface, each labyrinth occupies 50% of the volume.
- $t \neq 0$: unbalanced, one labyrinth grows and the other shrinks (a "network" solid when
  one side $F \le t$ is filled).
- Sheet solid of thickness set by $|F| \le c$ (a thickened surface).
- Graded structures: let $t = t(x,y,z)$ or $k = k(x,y,z)$ vary in space.

## 3. Differential-geometric properties (exact gyroid)

For an exact minimal surface:

- Mean curvature: $H = \tfrac{1}{2}(\kappa_1 + \kappa_2) = 0$
- Gaussian curvature: $K = \kappa_1 \kappa_2 \le 0$ (negative except at isolated flat points)
- Unit normal of the implicit form: $\mathbf{n} = \nabla F / \lVert \nabla F \rVert$, with

$$
\nabla F = \begin{pmatrix}
\cos x \cos y - \sin z \sin x \\
\cos y \cos z - \sin x \sin y \\
\cos z \cos x - \sin y \sin z
\end{pmatrix}
$$

- For an implicit surface, mean curvature is
  $H = \dfrac{\nabla F^{T}\,\mathrm{Hess}(F)\,\nabla F - \lVert\nabla F\rVert^{2}\,\mathrm{tr}\,\mathrm{Hess}(F)}{2\lVert\nabla F\rVert^{3}}$
  (sign convention depends on normal orientation).

## 4. Weierstrass–Enneper representation (exact surface)

Minimal surfaces can be written as

$$
\mathbf{x}(\omega) = \mathrm{Re}\int^{\omega} e^{i\theta}
\begin{pmatrix}
\tfrac12(1-\tau^2) \\ \tfrac{i}{2}(1+\tau^2) \\ \tau
\end{pmatrix} R(\tau)\, d\tau
$$

where $\theta$ is the Bonnet (associate-family) angle and $R(\tau)$ is a function
on the underlying Riemann surface. For the gyroid, as reported in the literature
(Schoen 1970; Große-Brauckmann & Wohlgemuth 1996):

$$
R(\tau) = \frac{1}{\sqrt{\tau^8 - 14\,\tau^4 + 1}}, \qquad \theta \approx 38.0147^\circ
$$

Associate family: Schwarz D ($\theta = 90^\circ$), gyroid ($\theta \approx 38.01^\circ$),
Schwarz P ($\theta = 0^\circ$) share the same $R(\tau)$. The gyroid is the only member
of this family (besides P and D) that is embedded.

> Verify the exact $R(\tau)$ and $\theta$ against the references below before citing;
> they are recalled from the literature, not re-derived here.

## 5. Symmetry and topology

- Space group of the surface: $I4_132$ (No. 214), cubic, chiral.
  The two labyrinths are mirror images of each other (enantiomorphic).
- Genus 3 per primitive (body-centered) unit cell.
- Related nodal TPMS for comparison:
  - Schwarz P: $\cos x + \cos y + \cos z = 0$
  - Schwarz D: $\sin x \sin y \sin z + \sin x \cos y \cos z + \cos x \sin y \cos z + \cos x \cos y \sin z = 0$

## 6. Useful approximate relations

- Surface area per unit cell (edge $L$, $t=0$): $A \approx 3.09\,L^2$ (numerical value
  commonly quoted for the cubic cell; verify numerically for your use).
- Volume fraction vs. $t$ is roughly linear for $|t|$ up to about 1.

## References

1. A. H. Schoen, *Infinite periodic minimal surfaces without self-intersections*, NASA Technical Note D-5541, 1970.
2. K. Große-Brauckmann, M. Wohlgemuth, *The gyroid is embedded and has constant mean curvature companions*, Calc. Var. 4, 499–523, 1996.
3. S. Hyde et al., *The Language of Shape*, Elsevier, 1997.
