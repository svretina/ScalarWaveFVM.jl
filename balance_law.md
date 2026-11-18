# **Energy balance for your evolved remainder fields**

For your actual coupled system, the natural energy diagnostic is

\\[ \boxed{ E\_{\mathrm r}(t) = \frac{1}{8\pi}\int\_\Omega \left(\Pi\_{\mathrm r}^2+ |\boldsymbol\Psi\_{\mathrm r}|^2\right)d^d x, \qquad d=1,3. } \\]

**This energy is not conserved by itself.** Its evolution is driven by the acceleration-dependent effective sources in your field equations.

We can derive an exact balance law involving only the evolved remainder variables, the analytical source coefficients, and the particle trajectory. We can then combine it with the particle energy equation to construct a quantity that should remain constant during your numerical evolution.

There is an important distinction:

- The **remainder-energy balance** follows exactly from your established field evolution equations.
- The **physical total conserved energy** additionally contains singular–remainder interaction energy and renormalized particle contributions.

These provide two different numerical tests. The first can be implemented immediately; the second tests whether your complete coupled system represents the correct physical scalar self-force.

I will use the consistent PPV normalization throughout:

\\[ \boxed{ \Box\Phi=-4\pi\rho,\qquad \rho=\frac q\gamma\delta^{(d)}(\mathbf x-\mathbf z), \qquad \Box=-\partial_t^2+\nabla^2. } \\]

## **1. Energy equation for the evolved fields**

Your decomposition is

\\[ \Phi=\Phi\_{\mathrm s}+\Phi\_{\mathrm r}, \qquad \Pi=\Pi\_{\mathrm s}+\Pi\_{\mathrm r}, \qquad \boldsymbol\Psi=\boldsymbol\Psi\_{\mathrm s}+\boldsymbol\Psi\_{\mathrm r}. \\]

Your remainder equations are

\\[ \boxed{ \begin{aligned} \partial_t\Pi\_{\mathrm r} &=\nabla\cdot\boldsymbol\Psi\_{\mathrm r} -\mathbf a\cdot\nabla\_{\mathbf v}\Pi\_{\mathrm s},\\\ \partial_t\boldsymbol\Psi\_{\mathrm r} &=\nabla\Pi\_{\mathrm r} -(\nabla\_{\mathbf v}\boldsymbol\Psi\_{\mathrm s})\mathbf a,\\\ \partial_t\Phi\_{\mathrm r} &=\Pi\_{\mathrm r} -\mathbf a\cdot\nabla\_{\mathbf v}\Phi\_{\mathrm s}. \end{aligned} } \tag{1} \\]

Define

\\[ u\_{\mathrm r} =\frac{\Pi\_{\mathrm r}^2+ |\boldsymbol\Psi\_{\mathrm r}|^2}{8\pi}. \\]

Differentiate:

\\[ \partial_tu\_{\mathrm r} =\frac1{4\pi} \left[ \Pi\_{\mathrm r}\partial_t\Pi\_{\mathrm r} +\boldsymbol\Psi\_{\mathrm r}\cdot \partial_t\boldsymbol\Psi\_{\mathrm r} \right]. \\]

Substitute (1):

\\[ \begin{aligned} 4\pi\partial_tu\_{\mathrm r} ={}& \Pi\_{\mathrm r}\nabla\cdot\boldsymbol\Psi\_{\mathrm r} +\boldsymbol\Psi\_{\mathrm r}\cdot\nabla\Pi\_{\mathrm r}\\\ &-\Pi\_{\mathrm r} \mathbf a\cdot\nabla\_{\mathbf v}\Pi\_{\mathrm s}\\\ &-\boldsymbol\Psi\_{\mathrm r}\cdot (\nabla\_{\mathbf v}\boldsymbol\Psi\_{\mathrm s})\mathbf a. \end{aligned} \\]

The first two terms form a divergence:

\\[ \Pi\_{\mathrm r}\nabla\cdot\boldsymbol\Psi\_{\mathrm r} +\boldsymbol\Psi\_{\mathrm r}\cdot\nabla\Pi\_{\mathrm r} = \nabla\cdot(\Pi\_{\mathrm r}\boldsymbol\Psi\_{\mathrm r}). \\]

Define the remainder energy flux

\\[ \mathbf S\_{\mathrm r} =-\frac1{4\pi}\Pi\_{\mathrm r}\boldsymbol\Psi\_{\mathrm r}. \\]

We obtain the exact local law

\\[ \boxed{ \partial_tu\_{\mathrm r} +\nabla\cdot\mathbf S\_{\mathrm r} = -\mathcal W\_{\mathrm{acc}}, } \tag{2} \\]

where

\\[ \boxed{ \mathcal W\_{\mathrm{acc}} = \frac1{4\pi} \left[ \Pi\_{\mathrm r} \mathbf a\cdot\nabla\_{\mathbf v}\Pi\_{\mathrm s} + \boldsymbol\Psi\_{\mathrm r}\cdot (\nabla\_{\mathbf v}\boldsymbol\Psi\_{\mathrm s})\mathbf a \right]. } \tag{3} \\]

This is the central result for your code.

**The acceleration-dependent effective sources exchange energy with the remainder fields at a rate&#x20;**\\(-\mathcal W\_{\mathrm{acc}}\\)**.**

There is no delta distribution in this energy equation because the point-particle source has already been eliminated from the remainder evolution equations.

The \\(\Phi\_{\mathrm r}\\) equation does not explicitly enter the quadratic field-energy density, but it is still needed to maintain the first-order constraints and for any scalar-dependent particle mass prescription.

## **2. Integrated remainder-energy balance**

Let

\\[ E\_{\mathrm r}(t)=\int\_\Omega u\_{\mathrm r}\\,d^dx, \\]

and define the outward remainder-field flux

\\[ \mathcal F\_{\mathrm r}(t) = \oint\_{\partial\Omega} \mathbf S\_{\mathrm r}\cdot\mathbf n\\,dA. \\]

Also define

\\[ \mathcal P\_{\mathrm{acc}}(t) = \int\_\Omega\mathcal W\_{\mathrm{acc}}\\,d^dx. \\]

Then

\\[ \boxed{ \frac{dE\_{\mathrm r}}{dt} +\mathcal F\_{\mathrm r} +\mathcal P\_{\mathrm{acc}} =0. } \tag{4} \\]

Integrating in time gives

\\[ \boxed{ E\_{\mathrm r}(t) +\int\_{t_0}^{t}\mathcal F\_{\mathrm r}(s)\\,ds +\int\_{t_0}^{t}\mathcal P\_{\mathrm{acc}}(s)\\,ds = E\_{\mathrm r}(t_0). } \tag{5} \\]

This is an exact conserved balance quantity for the remainder-field subsystem, provided the volume integrals and inner-worldline limits are well-defined.

It is directly applicable in both 1D and 3D.

For point-particle fields with singular limiting behavior, apply the identity first on a punctured domain or to the same finite-volume regularization used by your solver. An internal flux contribution must be included if the limiting remainder field carries nonvanishing energy flux through a shrinking surface around the particle.

## **3. Explicit one-dimensional balance law**

In genuine 1+1-dimensional scalar theory, the constant-velocity singular solution is

\\[ \Phi\_{\mathrm s}=-2\pi q\gamma|x-z|, \\]

so that

\\[ \Psi\_{\mathrm s}=-2\pi q\gamma\\,\operatorname{sgn}(x-z), \qquad \Pi\_{\mathrm s}=-v\Psi\_{\mathrm s}. \\]

For \\(x\ne z\\), its velocity derivatives at fixed \\(x-z\\) are

\\[ \frac{\partial\Psi\_{\mathrm s}}{\partial v} =-2\pi q\gamma^3v\\,\operatorname{sgn}(x-z), \\]

\\[ \frac{\partial\Pi\_{\mathrm s}}{\partial v} =2\pi q\gamma^3\\,\operatorname{sgn}(x-z). \\]

Substituting these into (3) gives

\\[ \boxed{ \mathcal W\_{\mathrm{acc}}^{1D} = \frac{qa\gamma^3}{2} \operatorname{sgn}(x-z) \left(\Pi\_{\mathrm r}-v\Psi\_{\mathrm r}\right). } \tag{6} \\]

Consequently, on \\(x_L\le x\le x_R\\),

\\[ \boxed{ \begin{aligned} \frac{dE\_{\mathrm r}^{1D}}{dt} ={}& \frac1{4\pi} [\Pi\_{\mathrm r}\Psi\_{\mathrm r}]\_{x_L}^{x_R}\\\ &-\frac{qa\gamma^3}{2} \int\_{x_L}^{x_R} \operatorname{sgn}(x-z) (\Pi\_{\mathrm r}-v\Psi\_{\mathrm r})\\,dx. \end{aligned} } \tag{7} \\]

Here

\\[ E\_{\mathrm r}^{1D} =\frac1{8\pi}\int\_{x_L}^{x_R} (\Pi\_{\mathrm r}^2+\Psi\_{\mathrm r}^2)\\,dx. \\]

The resulting time-integrated diagnostic is

\\[ \boxed{ \begin{aligned} \mathcal C\_{\mathrm r}^{1D}(t) ={}&E\_{\mathrm r}^{1D}(t)\\\ &-\frac1{4\pi}\int\_{t_0}^{t} [\Pi\_{\mathrm r}\Psi\_{\mathrm r}]\_{x_L}^{x_R}\\,ds\\\ &+\frac q2\int\_{t_0}^{t}a\gamma^3 \int\_{x_L}^{x_R} \operatorname{sgn}(x-z) (\Pi\_{\mathrm r}-v\Psi\_{\mathrm r})\\,dx\\,ds. \end{aligned} } \tag{8} \\]

For a sufficiently resolved evolution satisfying the corresponding continuum balance,

\\[ \boxed{ \mathcal C\_{\mathrm r}^{1D}(t)=E\_{\mathrm r}^{1D}(t_0). } \\]

This is a particularly useful test for your 1D code because every term is expressed in evolved variables and known particle quantities.

The static 1D field does not decay at infinity, so the finite computational domain and its boundary fluxes are essential.

## **4. Explicit three-dimensional balance law**

In 3+1 dimensions, define

\\[ \mathbf r=\mathbf x-\mathbf z,\qquad H=\mathbf r^2+\gamma^2(\mathbf v\cdot\mathbf r)^2. \\]

Your singular potential is

\\[ \Phi\_{\mathrm s}=\frac q{\sqrt H}. \\]

Its first-order variables are

\\[ \boxed{ \begin{aligned} \boldsymbol\Psi\_{\mathrm s} &=-q\frac{\mathbf r+\gamma^2(\mathbf v\cdot\mathbf r)\mathbf v}{H^{3/2}},\\\ \Pi\_{\mathrm s}&=-\mathbf v\cdot\boldsymbol\Psi\_{\mathrm s}. \end{aligned} } \tag{9} \\]

For implementation, the required velocity derivatives can be calculated explicitly.

Let

\\[ b=\mathbf v\cdot\mathbf r,\quad \mathbf N=\mathbf r+\gamma^2b\mathbf v,\quad \mathbf C=\gamma^2\mathbf r+\gamma^4b\mathbf v. \\]

Then

\\[ \frac{\partial H}{\partial v_j}=2bC_j, \\]

and

\\[ \frac{\partial N_i}{\partial v_j} =\gamma^2b\delta\_{ij} +\gamma^2v_i r_j +2\gamma^4bv_iv_j. \\]

Therefore,

\\[ \boxed{ \frac{\partial\Psi\_{\mathrm s,i}}{\partial v_j} = -q\left[ \frac{ \gamma^2b\delta\_{ij} +\gamma^2v_i r_j +2\gamma^4bv_iv_j }{H^{3/2}} -\frac{3bN_iC_j}{H^{5/2}} \right]. } \tag{10} \\]

And since \\(\Pi\_{\mathrm s}=-v_i\Psi\_{\mathrm s,i}\\),

\\[ \boxed{ \frac{\partial\Pi\_{\mathrm s}}{\partial v_j} = -\Psi\_{\mathrm s,j} -v_i\frac{\partial\Psi\_{\mathrm s,i}}{\partial v_j}. } \tag{11} \\]

Substitute these into

\\[ \mathcal W\_{\mathrm{acc}}^{3D} = \frac{a_j}{4\pi} \left[ \Pi\_{\mathrm r} \frac{\partial\Pi\_{\mathrm s}}{\partial v_j} +\Psi\_{\mathrm r,i} \frac{\partial\Psi\_{\mathrm s,i}}{\partial v_j} \right]. \\]

The exact balance law is

\\[ \boxed{ \begin{aligned} \frac{dE\_{\mathrm r}^{3D}}{dt} ={}& \frac1{4\pi} \oint\_{\partial\Omega} \Pi\_{\mathrm r}\boldsymbol\Psi\_{\mathrm r} \cdot\mathbf n\\,dA\\\ &-\frac1{4\pi} \int\_\Omega a_j \left[ \Pi\_{\mathrm r} \frac{\partial\Pi\_{\mathrm s}}{\partial v_j} +\Psi\_{\mathrm r,i} \frac{\partial\Psi\_{\mathrm s,i}}{\partial v_j} \right]d^3x. \end{aligned} } \tag{12} \\]

The first term measures energy transport through the computational boundaries. The second measures energy transferred between the remainder fields and the moving analytical puncture.

Equations (10)–(12) can be evaluated using precisely the same weakly integrated source terms as your field solver. In fact, doing so is preferable to introducing a separate numerical approximation for the energy-source integral.

## **5. Include the particle and external potential&#x20;**\\(V(\mathbf z)\\)

Now consider the particle state evolved by your RK integrator:

\\[ \dot{\mathbf z}=\mathbf v, \qquad \dot{\mathbf v}=\mathbf a. \\]

Let the finite, renormalized scalar-force inputs be

\\[ T=\Pi\_{\mathrm{ren},p}, \qquad \mathbf G=\boldsymbol\Psi\_{\mathrm{ren},p}. \\]

If your proposed interpolation closure is adopted, these are the interpolated numerical values of \\(\Pi\_{\mathrm r}\\) and \\(\boldsymbol\Psi\_{\mathrm r}\\). The physical validity of that closure remains a separate issue.

With a mechanical external potential \\(V(\mathbf z)\\), the consistent Cartesian acceleration equation is

\\[ \boxed{ \mathbf a = \frac{q}{M\gamma^2}(\mathbf G+\mathbf vT) -\frac1{M\gamma} (\mathbf I-\mathbf v\mathbf v^T)\nabla V. } \tag{13} \\]

This follows from the relativistic momentum equation

\\[ \frac{d}{dt}(\gamma M\mathbf v) = \frac q\gamma\mathbf G-\nabla V \\]

together with the scalar mass equation

\\[ \boxed{ \dot M=-q(T+\mathbf v\cdot\mathbf G). } \tag{14} \\]

The particle's energy, including the external mechanical potential, is

\\[ \boxed{ E\_{\mathrm p}=\gamma M+V(\mathbf z). } \tag{15} \\]

Differentiating and substituting (13)–(14) yields

\\[ \boxed{ \frac{dE\_{\mathrm p}}{dt} =-\frac q\gamma T. } \tag{16} \\]

This is the particle-side energy balance for the stated closure. It holds algebraically for your numerical particle equations if both the acceleration and mass equations are implemented consistently.

If instead your model holds \\(M\\) constant, equation (16) is not valid. The particle-energy rate must then be calculated from the actual evolution equations.

### **Combine the particle and remainder energies**

Adding (4) and (16),

\\[ \boxed{ \frac{d}{dt} \left(E\_{\mathrm r}+\gamma M+V\right) + \mathcal F\_{\mathrm r} + \mathcal P\_{\mathrm{acc}} + \frac q\gamma T =0. } \tag{17} \\]

Therefore the following cumulative quantity is conserved by the *stated coupled equations*:

\\[ \boxed{ \begin{aligned} \mathcal C\_{\mathrm{coupled}}(t) ={}&E\_{\mathrm r}(t)+\gamma(t)M(t)+V(\mathbf z(t))\\\ &+\int\_{t_0}^{t} \left[ \mathcal F\_{\mathrm r}(s) +\mathcal P\_{\mathrm{acc}}(s) +\frac q{\gamma(s)}T(s) \right]ds. \end{aligned} } \tag{18} \\]

Thus

\\[ \boxed{ \mathcal C\_{\mathrm{coupled}}(t) = E\_{\mathrm r}(t_0)+\gamma_0M_0+V(\mathbf z_0). } \tag{19} \\]

This is the most direct **coupled energy-balance diagnostic** for a code evolving \\((\Phi\_{\mathrm r},\Pi\_{\mathrm r},\boldsymbol\Psi\_{\mathrm r},\mathbf z,\mathbf v,M)\\) with (13)–(14).

It needs no explicit singular-field value at the particle. All the required puncture information enters through the same velocity-derivative source coefficients already used in the field equations.

However, \\(\mathcal C\_{\mathrm{coupled}}\\) is a cumulative balance invariant constructed from the evolution equations. It is **not automatically the physical conserved total energy**: that stronger statement requires accounting for singular–remainder interaction energy and matching the particle force and mass renormalization to the underlying scalar theory.

### **If your code uses constant mass**

For a constant \\(M\\), still using acceleration equation (13),

\\[ \frac{d}{dt}(\gamma M+V) =q\gamma\left( \mathbf v\cdot\mathbf G+\mathbf v^2T \right). \\]

The corresponding coupled diagnostic is instead

\\[ \boxed{ \begin{aligned} \mathcal C\_{\mathrm{coupled}}^{M=\mathrm{const}}(t) ={}&E\_{\mathrm r}+\gamma M+V\\\ &+\int\_{t_0}^{t} \Big[ \mathcal F\_{\mathrm r} +\mathcal P\_{\mathrm{acc}}\\\ &\hspace{31mm} -q\gamma(\mathbf v\cdot\mathbf G+\mathbf v^2T) \Big]ds. \end{aligned} } \tag{20} \\]

This is constant if the corresponding evolution equations are satisfied, subject to the same qualifications about inner boundaries and regularization.

## **6. What is the actual physical conserved energy?**

This is the crucial distinction for your self-force research.

The quantity \\(E\_{\mathrm r}+\gamma M+V\\) is not generally the physical total energy, even after mass renormalization.

To see why, expand the total field energy:

\\[ E\_\Phi =\frac1{8\pi}\int\_\Omega \left[ (\Pi\_{\mathrm s}+\Pi\_{\mathrm r})^2+ |\boldsymbol\Psi\_{\mathrm s}+\boldsymbol\Psi\_{\mathrm r}|^2 \right]d^dx. \\]

This gives

\\[ \boxed{ E\_\Phi=E\_{\mathrm r}+E\_{\mathrm s}+E\_{\mathrm{cross}}, } \tag{21} \\]

with

\\[ E\_{\mathrm s} =\frac1{8\pi}\int\_\Omega (\Pi\_{\mathrm s}^2+|\boldsymbol\Psi\_{\mathrm s}|^2)d^dx \\]

and

\\[ \boxed{ E\_{\mathrm{cross}} =\frac1{4\pi}\int\_\Omega \left( \Pi\_{\mathrm s}\Pi\_{\mathrm r} +\boldsymbol\Psi\_{\mathrm s}\cdot\boldsymbol\Psi\_{\mathrm r} \right)d^dx. } \tag{22} \\]

The cross term is particularly important. It represents energy associated with the interaction of your instantaneous boosted Coulomb field and the evolved remainder.

We can verify the exact exchange structure.

Define

\\[ C\_{\mathrm r} = \Pi\_{\mathrm r}\mathbf a\cdot\nabla\_{\mathbf v}\Pi\_{\mathrm s} + \boldsymbol\Psi\_{\mathrm r}\cdot (\nabla\_{\mathbf v}\boldsymbol\Psi\_{\mathrm s})\mathbf a \\]

and

\\[ C\_{\mathrm s} = \Pi\_{\mathrm s}\mathbf a\cdot\nabla\_{\mathbf v}\Pi\_{\mathrm s} + \boldsymbol\Psi\_{\mathrm s}\cdot (\nabla\_{\mathbf v}\boldsymbol\Psi\_{\mathrm s})\mathbf a. \\]

For a consistently regularized point source, the three field-energy sectors satisfy

\\[ \begin{aligned} \partial_tu\_{\mathrm r}+\nabla\cdot\mathbf S\_{\mathrm r} &=-\frac{C\_{\mathrm r}}{4\pi},\\\ \partial_tu\_{\mathrm s}+\nabla\cdot\mathbf S\_{\mathrm s} &=\rho\Pi\_{\mathrm s}+\frac{C\_{\mathrm s}}{4\pi},\\\ \partial_tu\_{\mathrm{cross}} +\nabla\cdot\mathbf S\_{\mathrm{cross}} &=\rho\Pi\_{\mathrm r} +\frac{C\_{\mathrm r}-C\_{\mathrm s}}{4\pi}. \end{aligned} \tag{23} \\]

Here

\\[ \mathbf S\_{\mathrm{cross}} =-\frac1{4\pi} (\Pi\_{\mathrm s}\boldsymbol\Psi\_{\mathrm r} +\Pi\_{\mathrm r}\boldsymbol\Psi\_{\mathrm s}). \\]

Adding these equations cancels every acceleration-dependent source:

\\[ \boxed{ \partial_tu\_\Phi+\nabla\cdot\mathbf S\_\Phi =\rho(\Pi\_{\mathrm s}+\Pi\_{\mathrm r}) =\rho\Pi. } \tag{24} \\]

This is an important consistency result for your ansatz.

**Your acceleration-dependent effective sources do not create or destroy physical energy. They redistribute energy among the instantaneous Coulomb field, the remainder field, and their cross term.**

The actual conserved physical quantity has the schematic form

\\[ \boxed{ E\_{\mathrm{physical}} = E\_{\mathrm r} +\gamma M+V +E\_{\mathrm{completion}}^{\mathrm{ren}}, } \tag{25} \\]

where \\(E\_{\mathrm{completion}}^{\mathrm{ren}}\\) contains the appropriately renormalized singular-field energy, cross energy, and any corresponding finite particle counterterm.

In 1D on a finite domain, these terms can often be handled directly with a specified regularization of the particle force. In 3D, the individual singular contributions diverge and must be combined with the bare particle energy before the regulator is removed.

We cannot simply set \\(E\_{\mathrm{completion}}^{\mathrm{ren}}=0\\). Nor is it established that the acceleration-source work integral in (18) equals that completion energy.

The physical energy balance ultimately takes the form

\\[ \boxed{ \frac{dE\_{\mathrm{physical}}}{dt} = -\oint\_{\partial\Omega} \mathbf S\_\Phi\cdot\mathbf n\\,dA, } \tag{26} \\]

with the **total** field flux, including singular and cross contributions.

## **7. What I recommend monitoring in your code**

I would implement two separate diagnostics.

| Diagnostic                                        | What it tests                                                                                                          |
| ------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| Remainder balance, equation (5)                   | Whether the discretized field equations correctly account for their acceleration-dependent sources and boundary fluxes |
| Coupled cumulative balance, equation (18) or (20) | Whether the numerical particle and field RHS satisfy their stated energy-exchange identities                           |

These are available directly from your evolved variables, with no need to reconstruct the singular total field at the particle.

For either diagnostic, calculate a normalized residual such as

\\[ \boxed{ \mathcal R_E(t)= \frac{ |\mathcal C(t)-\mathcal C(t_0)| }{ E\_{\mathrm{scale}} }, } \tag{27} \\]

where \\(E\_{\mathrm{scale}}\\) is a fixed positive reference energy, for example the initial remainder energy plus the absolute initial particle energy.

For the physical conservation test, additionally reconstruct the **total flux at the outer boundary** and investigate the renormalized completion term in (25).

Two qualifications matter:

- If your code does not evolve the scalar mass using (14), use the energy rate implied by its actual particle equations, not equation (18).
- If the continuum remainder is singular at the worldline, the integrated balance must include the corresponding inner-boundary flux or its finite-volume analogue.

### **Final result**

The immediately applicable balance law for the remainder fields you evolve is

\\[ \boxed{ \frac{d}{dt}\frac1{8\pi} \int\_\Omega \left(\Pi\_{\mathrm r}^2+ |\boldsymbol\Psi\_{\mathrm r}|^2\right)d^dx = \frac1{4\pi} \oint\_{\partial\Omega} \Pi\_{\mathrm r}\boldsymbol\Psi\_{\mathrm r}\cdot\mathbf n\\,dA - \frac1{4\pi} \int\_\Omega \left[ \Pi\_{\mathrm r}\mathbf a\cdot\nabla\_{\mathbf v}\Pi\_{\mathrm s} + \boldsymbol\Psi\_{\mathrm r}\cdot (\nabla\_{\mathbf v}\boldsymbol\Psi\_{\mathrm s})\mathbf a \right]d^dx. } \\]

This holds in both one and three spatial dimensions, with their respective singular fields and spatial measures.

For your numerical project, I would use this as the **first mandatory energy-convergence test**. It is derived directly from the equations you actually evolve.

The next, physically stronger step is to construct \\(E\_{\mathrm{completion}}^{\mathrm{ren}}\\) explicitly for your boosted-Coulomb ansatz, so that you have a genuinely conserved instantaneous total energy rather than only a cumulative balance diagnostic. That is also where the consistency of your proposed PPV particle-force closure can be tested against global energy conservation.&#x20;

[image](https://www.google.com/s2/favicons?domain=https://pmc.ncbi.nlm.nih.gov\&sz=32)

PubMed Central (PMC)
