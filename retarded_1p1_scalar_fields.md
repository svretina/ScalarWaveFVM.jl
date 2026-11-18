# Retarded scalar field in 1+1D flat spacetime

This note collects the retarded-field calculations for a scalar point particle in flat 1+1-dimensional spacetime, using the sign and normalization convention used in the code discussion.

We use coordinates

$$
(t,x),
$$

with metric signature $(-+ )$, so that

$$
\Box = -\partial_t^2 + \partial_x^2 .
$$

The code convention used below is equivalent to writing the scalar wave equation as

$$
(\partial_t^2 - \partial_x^2)\Phi(t,x)
= - q \sqrt{1-v(t)^2}\,\delta(x-z(t)),
$$

where the particle worldline is

$$
x=z(t),
\qquad
v(t)=\dot z(t),
\qquad |v(t)|<1.
$$

This differs from the common curved-spacetime convention

$$
\Box \Phi = -4\pi \mu
$$

by an overall sign and normalization. With the $4\pi$ convention, the fields below are multiplied by $-4\pi$. The structural formulas are the same.

---

## 1. Retarded Green function in 1+1D

The retarded Green function for

$$
(\partial_t^2-\partial_x^2)G_{\rm ret}(t,x;t',x')
= \delta(t-t')\delta(x-x')
$$

is

$$
G_{\rm ret}(t,x;t',x')
=
\frac12
\Theta(t-t')
\Theta\bigl(t-t'-|x-x'|\bigr).
$$

The important point is that in 1+1D the Green function has support inside the past light cone, not only on the light cone. Consequently the retarded scalar field contains a history integral over the worldline.

---

## 2. Point-particle source

The proper-time source density is

$$
\mu(t,x)=q\int \delta(t-z^0(\tau))\delta(x-z^1(\tau))\,d\tau .
$$

Using coordinate time $t'$ as the worldline parameter,

$$
z^\alpha(t')=(t',z(t')),
\qquad
 d\tau = \sqrt{1-v(t')^2}\,dt',
$$

so the source density becomes

$$
\mu(t,x)=q\sqrt{1-v(t)^2}\,\delta(x-z(t)).
$$

In the code convention, the wave equation is therefore

$$
(\partial_t^2-\partial_x^2)\Phi
= -q\sqrt{1-v(t)^2}\,\delta(x-z(t)).
$$

---

## 3. Exact retarded potential with switch-on time $t_0$

Suppose the source is switched on at time $t_0$. The exact retarded potential is

$$
\Phi_{\rm ret}(t,x)
=
-\frac{q}{2}
\int_{t_0}^{t}
\sqrt{1-v(\tau)^2}\,
\Theta\bigl(t-\tau-|x-z(\tau)|\bigr)
\,d\tau .
$$

If the field point $(t,x)$ is causally connected to the worldline after $t_0$, define the retarded time $t_r$ by

$$
t-t_r-|x-z(t_r)|=0,
\qquad t_r\ge t_0.
$$

Then, because the worldline is timelike, the function

$$
f(\tau)=t-\tau-|x-z(\tau)|
$$

is monotone decreasing across the relevant branch, and the Heaviside factor selects

$$
t_0\le \tau\le t_r.
$$

Thus

$$
\boxed{
\Phi_{\rm ret}(t,x)
=
-\frac{q}{2}
\int_{t_0}^{t_r}
\sqrt{1-v(\tau)^2}\,d\tau .
}
$$

If there is no causal contact after the switch-on time, i.e.

$$
t-t_0-|x-z(t_0)|<0,
$$

then

$$
\Phi_{\rm ret}(t,x)=0.
$$

For the common choice $t_0=0$, this is exactly what the code implements:

```julia
function Φ_ret_exact(x, p, t)
    f(tr) = t - tr - abs(x - z(tr, p))
    prob = NonlinearProblem((tr, _) -> f(tr), zero(t))
    tr = solve(prob; reltol=1e-18, abstol=1e-18).u
    tr = max(zero(tr), tr)

    integrand(τ) = sqrt(1 - v(τ, p)^2)
    I, _ = quadgk(integrand, zero(tr), tr; rtol=1e-18, atol=1e-18)

    return -p.q * I / 2
end
```

A more explicit causal guard is recommended:

```julia
if t - abs(x - z(0, p)) <= 0
    return zero(t)
end
```

assuming the source is switched on at $t_0=0$.

---

## 4. Exact retarded derivatives

Let

$$
s_r:=\operatorname{sign}(x-z(t_r)),
\qquad
v_r:=v(t_r),
\qquad
\alpha_r:=\sqrt{1-v_r^2}.
$$

The retarded time satisfies

$$
F(t,x,t_r)=t-t_r-|x-z(t_r)|=0.
$$

Differentiate implicitly. First,

$$
\frac{\partial F}{\partial t}=1,
\qquad
\frac{\partial F}{\partial t_r}=-1+s_r v_r.
$$

Therefore

$$
0=\partial_t F
=1+(-1+s_rv_r)\frac{\partial t_r}{\partial t},
$$

so

$$
\boxed{
\frac{\partial t_r}{\partial t}
=
\frac{1}{1-s_rv_r}.
}
$$

Similarly,

$$
\frac{\partial F}{\partial x}=-s_r,
$$

so

$$
0=\partial_x F
=-s_r+(-1+s_rv_r)\frac{\partial t_r}{\partial x},
$$

and therefore

$$
\boxed{
\frac{\partial t_r}{\partial x}
=
-\frac{s_r}{1-s_rv_r}.
}
$$

Since

$$
\Phi_{\rm ret}(t,x)
=-\frac{q}{2}\int_{t_0}^{t_r}\alpha(\tau)\,d\tau,
\qquad
\alpha(\tau)=\sqrt{1-v(\tau)^2},
$$

we get

$$
\partial_t\Phi_{\rm ret}
=-\frac{q}{2}\alpha_r\frac{\partial t_r}{\partial t},
$$

and

$$
\partial_x\Phi_{\rm ret}
=-\frac{q}{2}\alpha_r\frac{\partial t_r}{\partial x}.
$$

Thus

$$
\boxed{
\Pi_{\rm ret}
:=\partial_t\Phi_{\rm ret}
=
-\frac{q}{2}\frac{\alpha_r}{1-s_rv_r}
}
$$

and

$$
\boxed{
\Psi_{\rm ret}
:=\partial_x\Phi_{\rm ret}
=
+\frac{q}{2}\frac{s_r\alpha_r}{1-s_rv_r}.
}
$$

These expressions are exact. For accelerated motion, it is essential that the quantities are evaluated at the retarded time $t_r$, not at the current time $t$.

---

## 5. Static particle check

For a static particle at the origin,

$$
z(t)=0,
\qquad v(t)=0,
\qquad \alpha=1.
$$

With switch-on time $t_0=0$, the retarded time is

$$
t_r=t-|x|.
$$

Therefore

$$
\boxed{
\Phi_{\rm ret}(t,x)
=
-\frac{q}{2}(t-|x|)\Theta(t-|x|).
}
$$

The derivatives inside the future light cone are

$$
\boxed{
\Pi_{\rm ret}=-\frac{q}{2},
\qquad
\Psi_{\rm ret}=\frac{q}{2}\operatorname{sign}(x).
}
$$

This is why a purely static cusp potential

$$
\Phi_0(x)=\frac{q}{2}|x|
$$

does not reproduce the unshifted retarded $\Pi$. The full static retarded potential contains the additional term

$$
-\frac{q}{2}t.
$$

---

## 6. Uniform motion and Lorentz-boosted potential

For uniform motion,

$$
z(t)=z_0+vt,
\qquad v=\text{constant},
\qquad \alpha=\sqrt{1-v^2}.
$$

The exact retarded potential can be written in current-time form as

$$
\boxed{
\Phi_{\rm ret}(t,x)
=
\frac{q}{2}\alpha
\left[
\frac{|x-z(t)|}{1-v\sigma}-t+t_0
\right],
}
$$

where

$$
\sigma:=\operatorname{sign}(x-z(t)).
$$

For $t_0=0$, this becomes

$$
\boxed{
\Phi_{\rm ret}(t,x)
=
\frac{q}{2}\alpha
\left[
\frac{|x-z(t)|}{1-v\sigma}-t
\right].
}
$$

This is the source of the singular/local ansatz used in the code:

```julia
@inline function Φs(x, z, v, q, t)
    return 0.5 * q * sqrt(one(v) - v*v) *
           (abs(x - z) / (one(v) - v * sign(x - z)) - t)
end
```

The corresponding derivatives are

$$
\boxed{
\Pi_s
=
-\frac{q}{2}\frac{\alpha}{1-v\sigma},
}
$$

and

$$
\boxed{
\Psi_s
=
+\frac{q}{2}\frac{\sigma\alpha}{1-v\sigma}.
}
$$

In code:

```julia
@inline function Πs(x, z, v, q)
    σ = sign(x - z)
    α = sqrt(one(v) - v*v)
    return -0.5 * q * α / (one(v) - v*σ)
end

@inline function Ψs(x, z, v, q)
    σ = sign(x - z)
    α = sqrt(one(v) - v*v)
    return 0.5 * q * σ * α / (one(v) - v*σ)
end
```

For uniform motion, these satisfy

$$
\partial_t\Phi_s=\Pi_s,
\qquad
\partial_x\Phi_s=\Psi_s.
$$

They also satisfy the first-order reduced wave system with the particle source:

$$
\partial_t\Phi_s=\Pi_s,
$$

$$
\partial_t\Psi_s-\partial_x\Pi_s=0,
$$

$$
\partial_t\Pi_s-\partial_x\Psi_s
= -q\alpha\,\delta(x-z(t)).
$$

Equivalently,

$$
(\partial_t^2-\partial_x^2)\Phi_s
= -q\alpha\,\delta(x-z(t)).
$$

---

## 7. Why boosting only the static cusp is insufficient in 1+1D

The static spatial cusp is

$$
\Phi_0(X)=\frac{q}{2}|X|.
$$

Boosting only this gives the cusp-only potential

$$
\Phi_{\rm cusp}
=\frac{q}{2}\frac{\alpha |x-z(t)|}{1-v\sigma}.
$$

This satisfies the second-order wave equation for uniform motion, because the missing term is homogeneous. However,

$$
\partial_t\Phi_{\rm cusp}\ne \Pi_{\rm ret}.
$$

The missing term is

$$
-\frac{q}{2}\alpha t.
$$

The correct local retarded potential is therefore

$$
\Phi_s
=\Phi_{\rm cusp}-\frac{q}{2}\alpha t.
$$

The reason is specific to 1+1D: the Green function has support inside the light cone, so the static retarded solution is not just $q|X|/2$; it is

$$
\Phi_{\rm static,ret}(T,X)=\frac{q}{2}(|X|-T).
$$

Boosting the full rest-frame retarded expression gives the unshifted retarded derivatives. Boosting only $|X|$ gives the shifted/cusp-only fields.

---

## 8. Accelerated motion and the singular/local ansatz

For accelerated motion, the exact retarded field is still

$$
\Phi_{\rm ret}(t,x)
= -\frac{q}{2}\int_{t_0}^{t_r}\sqrt{1-v(\tau)^2}\,d\tau,
$$

with derivatives evaluated at $t_r$:

$$
\Pi_{\rm ret}
=
-\frac{q}{2}\frac{\alpha_r}{1-s_rv_r},
\qquad
\Psi_{\rm ret}
=
+\frac{q}{2}\frac{s_r\alpha_r}{1-s_rv_r}.
$$

The local singular ansatz instead uses current-time data:

$$
\Phi_s(t,x)=
\frac{q}{2}\alpha(t)
\left[
\frac{|x-z(t)|}{1-v(t)\sigma}-t
\right],
$$

$$
\Pi_s(t,x)=
-\frac{q}{2}\frac{\alpha(t)}{1-v(t)\sigma},
$$

$$
\Psi_s(t,x)=
+\frac{q}{2}\frac{\sigma\alpha(t)}{1-v(t)\sigma},
$$

where

$$
\sigma=\operatorname{sign}(x-z(t)).
$$

For accelerated motion this is not the exact retarded field. It is a frozen-velocity local singular field. The evolved regular field must carry the difference between the exact retarded history and this local ansatz.

---

## 9. Regular/singular decomposition

Let

$$
U=(\Phi,\Pi,\Psi).
$$

Decompose

$$
U_{\rm total}=U_r+U_s.
$$

For uniform motion, $U_s$ satisfies the sourced wave equation exactly, so $U_r$ satisfies a homogeneous equation.

For accelerated motion, $U_s$ depends on $v(t)$. Its time derivative contains an additional term

$$
\dot v\,\partial_v U_s=a\,\partial_v U_s.
$$

Therefore the regular variables should receive source terms of the form

$$
\boxed{
S_r=-a\,\partial_v U_s
}
$$

provided the decomposition is

$$
U_{\rm total}=U_r+U_s.
$$

This is the sign convention used in the code discussion.

---

## 10. Derivatives with respect to velocity

Let

$$
r=x-z,
\qquad
\sigma=\operatorname{sign}(r),
\qquad
\alpha=\sqrt{1-v^2},
\qquad
D=1-v\sigma.
$$

The singular fields are

$$
\Phi_s=\frac{q}{2}\alpha\left(\frac{|r|}{D}-t\right),
$$

$$
\Pi_s=-\frac{q}{2}\frac{\alpha}{D},
$$

$$
\Psi_s=\frac{q}{2}\frac{\sigma\alpha}{D}.
$$

Holding $x,z,t,q$ fixed, the useful identity is

$$
\partial_v\left(\frac{\alpha}{1-v\sigma}\right)
=
\frac{\sigma-v}{\alpha(1-v\sigma)^2}.
$$

Therefore

$$
\boxed{
\partial_v\Phi_s
=
\frac{q}{2}
\left[
|r|\frac{\sigma-v}{\alpha(1-v\sigma)^2}
+
\frac{vt}{\alpha}
\right],
}
$$

$$
\boxed{
\partial_v\Pi_s
=
-\frac{q}{2}\frac{\sigma-v}{\alpha(1-v\sigma)^2}
=
\frac{q}{2}\frac{v-\sigma}{\alpha(1-v\sigma)^2},
}
$$

and

$$
\boxed{
\partial_v\Psi_s
=
\frac{q}{2}\frac{1}{\alpha(1-v\sigma)}.
}
$$

Thus the acceleration sources are

$$
S_\Phi=-a\partial_v\Phi_s,
\qquad
S_\Pi=-a\partial_v\Pi_s,
\qquad
S_\Psi=-a\partial_v\Psi_s.
$$

---

## 11. Finite-volume cell averages in the particle cell

Consider a cell containing the particle:

$$
x_L<z<x_R,
\qquad
h=x_R-x_L.
$$

Define

$$
h_L=z-x_L,
\qquad
h_R=x_R-z,
\qquad
h=h_L+h_R.
$$

On the left side,

$$
\sigma=-1,
\qquad D_-=1+v.
$$

On the right side,

$$
\sigma=+1,
\qquad D_+=1-v.
$$

### 11.1 Average of $\Phi_s$

Because $\Phi_s$ is linear on each side of the particle, the exact average is the sum of two trapezoids:

$$
\boxed{
\langle\Phi_s\rangle_i
=
\frac{1}{2h}
\left[
 h_L\Phi_s(x_L)
 + h_R\Phi_s(x_R)
 + h\Phi_s(z)
\right].
}
$$

The particle value is

$$
\Phi_s(z)=-\frac{q}{2}\alpha t.
$$

This term is essential. Omitting it gives the correct average only for the cusp-only potential, not for the full retarded-compatible $\Phi_s$.

### 11.2 Average of $\Pi_s$

Since $\Pi_s$ is piecewise constant,

$$
\boxed{
\langle \Pi_s\rangle_i
=
\frac{h_L\Pi_-+h_R\Pi_+}{h},
}
$$

with

$$
\Pi_-=-\frac{q}{2}\frac{\alpha}{1+v},
\qquad
\Pi_+=-\frac{q}{2}\frac{\alpha}{1-v}.
$$

### 11.3 Average of $\Psi_s$

Similarly,

$$
\boxed{
\langle \Psi_s\rangle_i
=
\frac{h_L\Psi_-+h_R\Psi_+}{h},
}
$$

with

$$
\Psi_-=-\frac{q}{2}\frac{\alpha}{1+v},
\qquad
\Psi_+=+\frac{q}{2}\frac{\alpha}{1-v}.
$$

---

## 12. Cell averages of the acceleration sources

The finite-volume acceleration source for a cell is

$$
\bar S_U=-a\langle \partial_v U_s\rangle.
$$

For the particle-containing cell, define

$$
I_U^L=\int_{x_L}^{z}\partial_v U_s\,dx,
\qquad
I_U^R=\int_z^{x_R}\partial_v U_s\,dx.
$$

Then

$$
\bar S_U=-\frac{a}{h}(I_U^L+I_U^R).
$$

For $\partial_v\Phi_s$,

$$
\boxed{
I_\Phi^L
=
\frac{q}{2\alpha}
\left[
-\frac{h_L^2}{2(1+v)}+vt\,h_L
\right],
}
$$

and

$$
\boxed{
I_\Phi^R
=
\frac{q}{2\alpha}
\left[
\frac{h_R^2}{2(1-v)}+vt\,h_R
\right].
}
$$

For $\partial_v\Pi_s$,

$$
\boxed{
I_\Pi^L
=
\frac{q h_L}{2\alpha(1+v)},
}
$$

and

$$
\boxed{
I_\Pi^R
=
-\frac{q h_R}{2\alpha(1-v)}.
}
$$

For $\partial_v\Psi_s$,

$$
\boxed{
I_\Psi^L
=
\frac{q h_L}{2\alpha(1+v)},
}
$$

and

$$
\boxed{
I_\Psi^R
=
\frac{q h_R}{2\alpha(1-v)}.
}
$$

Thus

$$
\boxed{
\bar S_\Phi
=
-\frac{a q}{2\alpha h}
\left[
-\frac{h_L^2}{2(1+v)}
+
\frac{h_R^2}{2(1-v)}
+vt(h_L+h_R)
\right],
}
$$

$$
\boxed{
\bar S_\Pi
=
-\frac{a q}{2\alpha h}
\left[
\frac{h_L}{1+v}
-
\frac{h_R}{1-v}
\right],
}
$$

and

$$
\boxed{
\bar S_\Psi
=
-\frac{a q}{2\alpha h}
\left[
\frac{h_L}{1+v}
+
\frac{h_R}{1-v}
\right].
}
$$

Away from the particle cell, the same formulas reduce to ordinary one-sided values because $\sigma$ is constant across the cell.

---

## 13. Initial data and switch-on time

This is the most important practical point for comparisons.

If the exact retarded solution is defined with lower integration limit $t_0=0$, then the source is switched on at $t=0$. Away from the initial particle location,

$$
U_{\rm ret}(0,x)=0.
$$

Since

$$
U_{\rm total}=U_r+U_s,
$$

the regular initial data consistent with this switched-on retarded solution is

$$
\boxed{
U_r(0,x)=U_{\rm ret}(0,x)-U_s(0,x)=-U_s(0,x).
}
$$

In finite-volume form,

$$
\boxed{
\bar U_{r,i}(0)=-\bar U_{s,i}(0).
}
$$

This gives

$$
\bar U_{{\rm total},i}(0)
=\bar U_{r,i}(0)+\bar U_{s,i}(0)=0,
$$

matching the switched-on retarded solution.

If instead one chooses

$$
U_r(0,x)=0,
$$

then

$$
U_{\rm total}(0,x)=U_s(0,x).
$$

This corresponds to starting from a preloaded instantaneous singular/local field, not to a source switched on at $t=0$. It is a different Cauchy problem. The corresponding total solution is the switched-on retarded solution plus a homogeneous solution generated by the initial mismatch.

Therefore:

- To compare with the exact retarded fields whose history integral starts at $t_0=0$, use

  $$
  U_r(0)=-U_s(0).
  $$

- If using

  $$
  U_r(0)=0,
  $$

  then do not compare against $U_{\rm ret}^{(t_0=0)}-U_s$ unless the additional homogeneous contribution is also included.

---

## 14. What to compare in a finite-volume code

Since the evolved variables are finite-volume cell averages, comparisons should also use cell averages.

For example, for $\Psi_r$, compare

$$
\boxed{
\bar\Psi_{r,i}^{\rm exact}
=
\frac{1}{h_i}\int_{x_{i-1/2}}^{x_{i+1/2}}
\left[\Psi_{\rm ret}(t,x)-\Psi_s(t,x)\right]dx
}
$$

against the numerical cell average $\bar\Psi_{r,i}^{\rm num}$.

Similarly, the total field comparison should use

$$
\boxed{
\bar\Psi_{{\rm total},i}^{\rm num}
=
\bar\Psi_{r,i}^{\rm num}+\bar\Psi_{s,i}
}
$$

and compare it against

$$
\boxed{
\bar\Psi_{{\rm ret},i}^{\rm exact}
=
\frac{1}{h_i}\int_{x_{i-1/2}}^{x_{i+1/2}}
\Psi_{\rm ret}(t,x)\,dx.
}
$$

Pointwise analytic values are adequate only in smooth regions away from the particle and away from causal fronts.

---

## 15. Causal fronts and convergence measurements

For a source switched on at $t_0=0$, the retarded solution has causal fronts near

$$
x=z(0)\pm t.
$$

The regular field

$$
U_r=U_{\rm ret}-U_s
$$

inherits sharp features at these fronts. Therefore global convergence norms over regions that include the fronts will generally not show clean second-order convergence, even if the scheme is second-order in smooth regions.

For convergence diagnostics, use one of the following:

1. Exclude fixed-width neighborhoods of the causal fronts.
2. Exclude a guard region around the particle position.
3. Compare to cell-averaged analytic data.
4. Use self-convergence with three resolutions.
5. Report separate smooth-region and global convergence rates.

A typical smooth-region mask is

$$
|x-z(t)|>\epsilon,
\qquad
|x-(z(0)-t)|>\epsilon,
\qquad
|x-(z(0)+t)|>\epsilon,
$$

with $\epsilon$ chosen large enough to exclude the numerical front transition region.

---

## 16. Summary

The exact retarded potential in the code convention is

$$
\boxed{
\Phi_{\rm ret}(t,x)
=
-\frac{q}{2}
\int_{t_0}^{t_r}\sqrt{1-v(\tau)^2}\,d\tau,
}
$$

where

$$
\boxed{
t-t_r-|x-z(t_r)|=0.
}
$$

The exact retarded derivatives are

$$
\boxed{
\Pi_{\rm ret}
=
-\frac{q}{2}\frac{\sqrt{1-v(t_r)^2}}{1-s_rv(t_r)},
}
$$

and

$$
\boxed{
\Psi_{\rm ret}
=
+\frac{q}{2}\frac{s_r\sqrt{1-v(t_r)^2}}{1-s_rv(t_r)},
}
$$

with

$$
s_r=\operatorname{sign}(x-z(t_r)).
$$

The local singular ansatz is exact for uniform motion and is given by

$$
\boxed{
\Phi_s
=
\frac{q}{2}\sqrt{1-v^2}
\left[
\frac{|x-z|}{1-v\operatorname{sign}(x-z)}-t
\right],
}
$$

$$
\boxed{
\Pi_s
=
-\frac{q}{2}\frac{\sqrt{1-v^2}}{1-v\operatorname{sign}(x-z)},
}
$$

$$
\boxed{
\Psi_s
=
+\frac{q}{2}\frac{\operatorname{sign}(x-z)\sqrt{1-v^2}}{1-v\operatorname{sign}(x-z)}.
}
$$

For accelerated motion, $U_s$ is a frozen-velocity singular approximation. The regular variables should be evolved with acceleration sources

$$
\boxed{
S_r=-a\partial_v U_s.
}
$$

For a retarded solution switched on at $t_0=0$, the consistent regular initial data is

$$
\boxed{
U_r(0)=-U_s(0).
}
$$

The alternative choice $U_r(0)=0$ corresponds to a preloaded local-field Cauchy problem, not to the switched-on retarded solution.
