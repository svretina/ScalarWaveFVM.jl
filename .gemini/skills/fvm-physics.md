---
name: FVM Physics Core
description: Local domain skills for 1D FVM solver with particle coupling.
---
# FVM and Field-Particle Coupling

- The domain is 1D scalar wave field coupled to a relativistic particle.
- The balance law involves discrete sums representing integrals.
- Work term uses $W_{src,h} = h \sum (\bar\Pi_r \bar S_\Pi + \bar\Psi_r \bar S_\Psi)$.
- Interpolation must be used to find acceleration from velocity over time.
- Uses `PCHIPInterpolation` or `AkimaInterpolation` from `DataInterpolations` and `DataInterpolations.derivative`.
