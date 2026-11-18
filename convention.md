# Scalar Particle Coupled to a 1+1D Scalar Field: Choice A vs Choice B

## Setup

Worldline:
$$
z=z(t),\qquad v=\dot z(t),\qquad \gamma=\frac{1}{\sqrt{1-v^2}}.
$$

Point-particle source:
$$
\rho(t,x)=\frac{q}{\gamma}\,\delta(x-z(t)).
$$

External mechanical potential:
$$
-\,V(z)
$$
inside the particle Lagrangian.

If a split field is used:
$$
\Phi=\Phi_s+\Phi_r,
$$
where $\Phi_s$ is the prescribed singular/lab-frame piece and $\Phi_r$ is the evolved remainder field.

---

## Choice A: action-field convention

This is the convention with a **plus** interaction term.

### Lagrangian
$$
L
=
\int dx\,\frac12\big(\Phi_t^2-\Phi_x^2\big)
+\int dx\,\rho\,\Phi
-m_0\sqrt{1-v^2}
-V(z).
$$

Equivalent particle form:
$$
L_{\rm part}=-(m_0-q\Phi(t,z))\sqrt{1-v^2}-V(z).
$$

### Field equation
$$
\Phi_{tt}-\Phi_{xx}=\rho.
$$

### Static / Poisson equation
$$
-\Phi_{xx}=\rho,
\qquad\text{equivalently}\qquad
\Phi_{xx}=-\rho.
$$

### Particle energy
$$
E_p=(m_0-q\Phi)\gamma+V(z).
$$

With a split $\Phi=\Phi_s+\Phi_r$:
$$
E_p=(m_0-q\Phi_s-q\Phi_r)\gamma+V(z).
$$

### Effective-mass form
If the singular piece is absorbed into an effective mass:
$$
m_{\rm eff}=m_0-q\Phi_s^{\rm reg}(z),
$$
then
$$
E_p=m_{\rm eff}\gamma-q\Phi_r(z)\gamma+V(z).
$$

### Potential-energy interpretation
$$
U_\Phi=-q\Phi.
$$

For $q>0$, larger $\Phi$ means **lower** particle potential energy.

### Particle momentum equation
$$
\frac{d}{dt}\Big[(m_0-q\Phi)\gamma v\Big]
=
-\,V'(z)+\frac{q}{\gamma}\Phi_x(t,z).
$$

### Particle acceleration equation
$$
\ddot z
=
\frac{
-\,V'(z)(1-v^2)^{3/2}
+q(1-v^2)\big(\Phi_x+v\Phi_t\big)
}{m_0-q\Phi(t,z)}.
$$

If first-order variables are used,
$$
\Pi=\Phi_t,\qquad \Psi=\Phi_x,
$$
then
$$
\ddot z
=
\frac{
-\,V'(z)(1-v^2)^{3/2}
+q(1-v^2)\big(\Psi+v\Pi\big)
}{m_0-q\Phi}.
$$

### Static nonrelativistic limit
$$
m_0\ddot z\approx -V'(z)+q\Phi_x.
$$

### Particle energy-rate law
$$
\frac{dE_p}{dt}=-\frac{q}{\gamma}\,\Phi_t(t,z(t)).
$$

With remainder variables:
$$
\dot E_p\supset -\frac{q}{\gamma}\Pi_r(z)
$$
for the regular-field contribution.

### Field energy
$$
E_{\rm field}=\int dx\,\frac12\big(\Phi_t^2+\Phi_x^2\big).
$$

### Field energy flux
For
$$
\Pi=\Phi_t,\qquad \Psi=\Phi_x,
$$
the 1D energy flux is
$$
\mathcal F=-\Pi\Psi.
$$

### Field energy balance
$$
\frac{dE_{\rm field}}{dt}=-\mathcal F_{\rm out}+q\Pi(z).
$$

### Total energy
$$
E_{\rm tot}=E_p+E_{\rm field}
$$
for the closed coupled system, up to numerical error and boundary fluxes.

---

## Choice B: lab-potential convention

This is the convention with a **minus** interaction term, where the field behaves like an ordinary potential-energy-per-charge variable.

### Lagrangian
$$
L
=
\int dx\,\frac12\big(\Phi_t^2-\Phi_x^2\big)
-\int dx\,\rho\,\Phi
-m_0\sqrt{1-v^2}
-V(z).
$$

Equivalent particle form:
$$
L_{\rm part}=-(m_0+q\Phi(t,z))\sqrt{1-v^2}-V(z).
$$

### Field equation
$$
\Phi_{tt}-\Phi_{xx}=-\rho.
$$

### Static / Poisson equation
$$
-\Phi_{xx}=-\rho,
\qquad\text{equivalently}\qquad
\Phi_{xx}=\rho.
$$

### Particle energy
$$
E_p=(m_0+q\Phi)\gamma+V(z).
$$

With a split $\Phi=\Phi_s+\Phi_r$:
$$
E_p=(m_0+q\Phi_s+q\Phi_r)\gamma+V(z).
$$

### Effective-mass form
If the singular piece is absorbed into an effective mass:
$$
m_{\rm eff}=m_0+q\Phi_s^{\rm reg}(z),
$$
then
$$
E_p=m_{\rm eff}\gamma+q\Phi_r(z)\gamma+V(z).
$$

### Potential-energy interpretation
$$
U_\Phi=+q\Phi.
$$

For $q>0$, larger $\Phi$ means **higher** particle potential energy.

### Particle momentum equation
$$
\frac{d}{dt}\Big[(m_0+q\Phi)\gamma v\Big]
=
-\,V'(z)-\frac{q}{\gamma}\Phi_x(t,z).
$$

### Particle acceleration equation
$$
\ddot z
=
\frac{
-\,V'(z)(1-v^2)^{3/2}
-q(1-v^2)\big(\Phi_x+v\Phi_t\big)
}{m_0+q\Phi(t,z)}.
$$

If first-order variables are used,
$$
\Pi=\Phi_t,\qquad \Psi=\Phi_x,
$$
then
$$
\ddot z
=
\frac{
-\,V'(z)(1-v^2)^{3/2}
-q(1-v^2)\big(\Psi+v\Pi\big)
}{m_0+q\Phi}.
$$

### Static nonrelativistic limit
$$
m_0\ddot z\approx -V'(z)-q\Phi_x.
$$

### Particle energy-rate law
$$
\frac{dE_p}{dt}=\frac{q}{\gamma}\,\Phi_t(t,z(t)).
$$

With remainder variables:
$$
\dot E_p\supset \frac{q}{\gamma}\Pi_r(z)
$$
for the regular-field contribution.

### Field energy
$$
E_{\rm field}=\int dx\,\frac12\big(\Phi_t^2+\Phi_x^2\big).
$$

### Field energy flux
For
$$
\Pi=\Phi_t,\qquad \Psi=\Phi_x,
$$
the 1D energy flux is
$$
\mathcal F=-\Pi\Psi.
$$

### Field energy balance
$$
\frac{dE_{\rm field}}{dt}=-\mathcal F_{\rm out}-q\Pi(z).
$$

### Total energy
$$
E_{\rm tot}=E_p+E_{\rm field}
$$
for the closed coupled system, up to numerical error and boundary fluxes.

---

## Side-by-side dictionary

### Interaction sign
Choice A:
$$
L_{\rm int}=+\int dx\,\rho\Phi
$$

Choice B:
$$
L_{\rm int}=-\int dx\,\rho\Phi
$$

### Wave equation
Choice A:
$$
\Box\Phi=\rho
$$

Choice B:
$$
\Box\Phi=-\rho
$$

### Static equation
Choice A:
$$
\Phi_{xx}=-\rho
$$

Choice B:
$$
\Phi_{xx}=\rho
$$

### Potential-energy interpretation
Choice A:
$$
U_\Phi=-q\Phi
$$

Choice B:
$$
U_\Phi=+q\Phi
$$

### Effective mass from singular piece
Choice A:
$$
m_{\rm eff}=m_0-q\Phi_s^{\rm reg}
$$

Choice B:
$$
m_{\rm eff}=m_0+q\Phi_s^{\rm reg}
$$

### Static force law
Choice A:
$$
F_\Phi=+q\Phi_x
$$

Choice B:
$$
F_\Phi=-q\Phi_x
$$

### Remainder-field power term
Choice A:
$$
\dot E_p\supset -\frac{q}{\gamma}\Pi_r(z)
$$

Choice B:
$$
\dot E_p\supset +\frac{q}{\gamma}\Pi_r(z)
$$

---

## Recommended bookkeeping with a split field

If $\Phi_s$ is prescribed and $\Phi_r$ is evolved, the cleanest setups are:

### Choice A
Keep bare mass:
$$
m=m_0,
\qquad
E_p=m_0\gamma-q\Phi_s^{\rm reg}\gamma-q\Phi_r\gamma+V,
$$
or absorb only the singular piece:
$$
m=m_0-q\Phi_s^{\rm reg},
\qquad
E_p=m\gamma-q\Phi_r\gamma+V.
$$

### Choice B
Keep bare mass:
$$
m=m_0,
\qquad
E_p=m_0\gamma+q\Phi_s^{\rm reg}\gamma+q\Phi_r\gamma+V,
$$
or absorb only the singular piece:
$$
m=m_0+q\Phi_s^{\rm reg},
\qquad
E_p=m\gamma+q\Phi_r\gamma+V.
$$

### What to avoid
Do not absorb the same field piece into $m$ and also add it explicitly again in the energy or force diagnostics. That double-counts that piece.
