# ScalarWaveFVM

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://svretina.github.io/ScalarWaveFVM.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://svretina.github.io/ScalarWaveFVM.jl/dev/)
[![Build Status](https://github.com/svretina/ScalarWaveFVM.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/svretina/ScalarWaveFVM.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/svretina/ScalarWaveFVM.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/svretina/ScalarWaveFVM.jl)


# Validated 1+1 scalar self-force reference solver

The current reference implementation uses a retarded-matched puncture,
proper-time and dynamical-mass evolution, MUSCL finite volumes, and independent
controls for field sourcing and particle feedback. The sign and normalization
contract is documented in [`convention.md`](convention.md); continuum and exact
semidiscrete energy identities are documented in
[`balance_law.md`](balance_law.md).

Run the complete validation suite with:

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

The legacy one- and two-particle plotting/energy modules were retired because
they encoded the previous mass sign, obsolete puncture signatures, and
incompatible normalization. The supported diagnostic API is `Diagnostics`;
retarded-field plots remain available through `enable_plotting!()`.
 
