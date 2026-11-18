# Energy balance for the validated 1+1 remainder system

The convention is

\[
\Box=-\partial_t^2+\partial_x^2,\qquad
\Box\Phi={q\over\gamma}\delta(x-z),\qquad
\Pi=\partial_t\Phi,\quad\Psi=\partial_x\Phi.
\]

The evolved remainder satisfies

\[
\partial_t\Pi_r=\partial_x\Psi_r+S_\Pi,
\qquad
\partial_t\Psi_r=\partial_x\Pi_r+S_\Psi,
\]

with the exact finite-volume averages of

\[
S_\Pi={qa\gamma^3\over2}(v+n),\qquad
S_\Psi=-{qa\gamma^3\over2}(1+vn).
\]

## Continuum identity

For

\[
E_r={1\over2}\int(\Pi_r^2+\Psi_r^2)\,dx,
\]

the continuum balance is

\[
{dE_r\over dt}=[\Pi_r\Psi_r]_{x_L}^{x_R}
+\int(\Pi_rS_\Pi+\Psi_rS_\Psi)\,dx.
\]

## Exact periodic semidiscrete identity

Let \(Q_i=(\Pi_i,\Psi_i)^T\), \(A=\begin{pmatrix}0&-c^2\\-1&0\end{pmatrix}\),
and \(H=\operatorname{diag}(c^{-2},1)\). The finite-volume energy is

\[
E_h={h\over2}\sum_i Q_i^THQ_i.
\]

For Godunov or Rusanov flux, define the cell-average first-order flux

\[
F^{(1)}_{i+1/2}=A{Q_i+Q_{i+1}\over2}-{c\over2}(Q_{i+1}-Q_i)
\]

and let \(\widehat F_{i+1/2}\) be the actual flux evaluated on the MUSCL
face states. On the periodic grid used by the solver,

\[
\boxed{
{dE_h\over dt}=P_{S,h}-D_h+R_h
}
\]

where

\[
P_{S,h}=h\sum_iQ_i^THS_i,
\]

\[
D_h={c\over2}\sum_i(\Delta Q_i)^TH(\Delta Q_i)\ge0,
\]

and

\[
R_h=\sum_i(\Delta Q_i)^TH
\left(\widehat F_{i+1/2}-F^{(1)}_{i+1/2}\right).
\]

For piecewise-constant reconstruction, \(R_h=0\). With MUSCL,
\(R_h\) is the exact reconstruction correction; it is not assumed to have a
sign. `Diagnostics.semidiscrete_energy_terms` evaluates this decomposition
using the solver's actual reconstruction, numerical flux, cell averages, and
source discretization. The regression test verifies it directly against
\(h\sum_i Q_i^TH\dot Q_i\).

Nonperiodic boundaries require the corresponding numerical entropy-flux term.
The current production spatial operator is periodic, so its net boundary term
is zero.
