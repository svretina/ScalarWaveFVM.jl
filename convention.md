# Scalar-field convention used by the 1D solver

This file is the sign and normalization contract for the current one-particle
implementation. Symbols are

\[
\Pi=\partial_t\Phi,\qquad \Psi=\partial_x\Phi,\qquad
\gamma=(1-v^2)^{-1/2}.
\]

## Field equation and particle coupling

The code uses the signature

\[
\Box=-\partial_t^2+\partial_x^2
\]

and the sourced equation

\[
\Box\Phi={q\over\gamma}\delta(x-z(t)).
\]

and the particle action

\[
S_p=-\int (m_0+q\Phi)\,d\tau-\int V(z)\,dt.
\]

Thus the scalar contribution to the coordinate acceleration has the negative
sign

\[
\dot v=-{q\over m\gamma^2}(\Psi_r+v\Pi_r)
       -{V'(z)\over m\gamma^3}.
\]

The evolved dynamical mass and proper time obey

\[
\dot m=q(\Pi_r+v\Psi_r)-{q^2\over2\gamma},\qquad
\dot\tau={1\over\gamma}.
\]

All fields in these particle equations are regular remainder fields
reconstructed at the particle by `ODE.regular_fields_at_particle`.

## Retarded-matched puncture

With \(r=x-z(t)\), \(n=\operatorname{sign}(r)\), and accumulated proper time
\(\tau(t)\), the singular field is

\[
\Phi_s={q\gamma\over2}(|r|+vr)-{q\over2}\tau,
\]
\[
\Pi_s=-{q\gamma\over2}(1+vn),\qquad
\Psi_s={q\gamma\over2}(n+v).
\]

Velocity derivatives in the remainder source hold \(r\) and \(\tau\) fixed:

\[
\partial_v\Phi_s={q\gamma^3\over2}(r+v|r|),
\]
\[
\partial_v\Pi_s=-{q\gamma^3\over2}(v+n),\qquad
\partial_v\Psi_s={q\gamma^3\over2}(1+vn).
\]

The finite-volume equations insert the exact cell averages of
\(-a\,\partial_v(\Pi_s,\Psi_s,\Phi_s)\).

## Actual retarded solution

The puncture is not obtained by evaluating itself at retarded particle data.
For a source switched on at \(t=0\), the retarded time satisfies

\[
t-t_r=|x-z(t_r)|,
\]

and the scalar potential is

\[
\Phi_{\rm ret}(t,x)=-{q\over2}\int_0^{t_r}{dt'\over\gamma(t')}.
\]

`RetardedPotentials.Φ_ret`, `Π_ret`, and `Ψ_ret` implement these retarded-time
expressions and return zero outside the future light cone of the switch-on
event.

## Initialization and coupling controls

`ScalarField.initialize_remainder` provides three explicit choices:

- `:zero_remainder`: the full initial field is the selected global puncture
  extension. This is not the switched-on retarded solution.
- `:retarded`: a supplied full retarded field is sampled and the puncture is
  subtracted.
- `:custom`: user-supplied `(Πr, Ψr, Φr)` data are copied.

Coupling direction is controlled independently:

- `field_source=true` lets the trajectory drive the remainder field.
- `particle_feedback=true` lets the regular field drive the particle and
  includes retarded self-mass loss.

Thus a prescribed trajectory uses a field source without particle feedback;
turning off feedback does not silently remove the particle's field source.
