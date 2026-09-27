# mhdturbFoam-hardened

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23003542.svg)](https://doi.org/10.5281/zenodo.23003542)

**A hardened redistribution of the FOSSEE `mhdturbFoam` turbulent MHD solver for
OpenFOAM 6, with three numerical defects fixed and a verification suite.**

> **This is a derived work, not an original solver.** The physical formulation
> and the numerical structure come from OpenFOAM's `mhdFoam`; the turbulent
> skeleton comes from the FOSSEE `mhdturbFoam` derivative. The contribution of
> this package is the correction and the packaging — see
> [`NOTICE`](NOTICE) and [`docs/corrections.md`](docs/corrections.md).

The solver solves the incompressible Navier–Stokes equations coupled to the
**full magnetic-induction** MHD formulation:

$$
\frac{\partial \mathbf{u}}{\partial t} + \nabla\cdot(\mathbf{u}\mathbf{u})
= -\nabla p + \nabla\cdot(\nu_\text{eff}\nabla\mathbf{u})
- \frac{1}{\rho}\nabla\cdot\left(\frac{\mathbf{B}\mathbf{B}}{\mu_0}\right)
+ \frac{1}{\rho}\nabla\left(\frac{|\mathbf{B}|^2}{2\mu_0}\right),
\qquad \nabla\cdot\mathbf{u} = 0
$$

$$
\frac{\partial \mathbf{B}}{\partial t} + \nabla\cdot(\mathbf{u}\mathbf{B})
- \nabla\cdot(\mathbf{B}\mathbf{u}) = \frac{1}{\mu_0\sigma}\nabla^2\mathbf{B},
\qquad \nabla\cdot\mathbf{B} = 0
$$

The Lorentz force is evaluated in **Maxwell-stress form from the total field**,
exactly as in `mhdFoam`. The divergence constraint on $\mathbf{B}$ is enforced
by a magnetic flux pressure `pB` (B-PISO loop).

## What is corrected here

The upstream derivative (`mhdturbFoam`, FOSSEE 2019) changes `mhdFoam` by
replacing the laminar `fvm::laplacian(nu, U)` with `turbulence->divDevReff(U)`
and adding a transport/turbulence model. In doing so it introduced three
numerical defects. All three are fixed in this package; the fixes are the only
differences from the upstream source:

1. **Restored outer PISO corrector loop.** The upstream version replaced
   `while (piso.correct())` with a plain block, so `rAU`, `HbyA` and `phiHbyA`
   were computed once per time step and `piso.finalInnerIter()` was never true —
   the `pFinal` solver dictionary was unreachable and every pressure solve
   stopped at `relTol 0.05`. With the loop restored, `pFinal` is reached in the
   last corrector of each step (see the verification table below: final pressure
   residual $\approx 10^{-6}$ instead of $\approx 10^{-3}$).
2. **Reference level for the singular `pB` problem.** In a closed or fully
   periodic domain every `pB` boundary is `zeroGradient` or `cyclic`, so the
   `pB` equation is a purely Neumann Poisson problem and is singular. The solver
   now reads an optional `pBRefCell`/`pBRefValue` from the `BPISO` dictionary
   and calls `pBEqn.setReference(...)`. Without it the run dies with a floating
   point exception. The entries are optional and default to cell 0 / value 0, so
   open domains with a `fixedValue` `pB` outlet are unaffected.
3. **`fvOptions` coupling in the momentum equation.** `fvOptions(U)` is now a
   source term of `UEqn` and `fvOptions.constrain(UEqn)` is called. Without the
   latter, `fv::meanVelocityForce` never initialises `rAPtr_` and cannot drive a
   periodic channel.

## Layout

```
mhdturbFoam-hardened/
├── src/                        solver source (build with wmake)
│   ├── mhdturbFoam.C           corrected solver
│   ├── createFields.H          fields, transport, turbulence, fvOptions
│   ├── createControl.H         piso / bpiso controls
│   ├── createPhiB.H            magnetic face flux
│   ├── magneticFieldErr.H      div(B) diagnostic
│   ├── readBPISOControls.H     BPISO helpers
│   └── Make/                   wmake files
├── cases/                      verification, with objective metrics
│   ├── channelHartmann/        upstream reference case: SA Cv1=30, B = (0 8 0)
│   ├── channelHartmannSpanwise/ same, spanwise B = (0 0 -8)
│   └── channelHartmannPeriodic/ streamwise-periodic, force-driven (fvOptions)
├── examples/                   integration/smoke tests, no exact solution
│   └── swirlDuctCoarse/        polygonal duct, k-omega SST, B = (0 0 -10), ~30 s
├── docs/
│   ├── corrections.md          the three defects in detail, with evidence
│   └── verification.md         generated metrics (verify.sh --report)
├── verify.sh                   runs every case and reports the metrics
├── NOTICE                      attribution and licensing
├── LICENSE                     GPL-3.0
├── CITATION.cff
└── .zenodo.json
```

## Requirements

* OpenFOAM 6 (openfoam.org), including `libturbulenceModels`,
  `libincompressibleTurbulenceModels`, `libincompressibleTransportModels` and
  `libfvOptions`.
* A C++11 compiler.

## Build

```bash
source /path/to/OpenFOAM-6/etc/bashrc
cd src
wmake
```

The executable is installed in `$FOAM_USER_APPBIN/mhdturbFoam` and can be used
as a drop-in replacement for the upstream solver (`application mhdturbFoam;` in
`system/controlDict`).

## Run

Each case is a standard OpenFOAM case with `Allrun`/`Allclean`:

```bash
cd cases/channelHartmann
./Allrun
```

All cases use the geometry and boundary-condition layout of the standard
OpenFOAM 6 `mhdFoam` tutorial `electromagnetics/mhdFoam/hartmann` (channel
$0 \le x \le 20$, $-1 \le y \le 1$, $0 \le z \le 0.1$, mesh $100\times40\times1$,
4000 cells, patches `inlet`/`outlet`/`lowerWall`/`upperWall`/`frontAndBack`),
the same geometry as the upstream `mhdturbFoam` reference case:

| case | turbulence | `B` | driving | steps |
|---|---|---|---|---|
| `channelHartmann` | **RAS Spalart–Allmaras, `Cv1 30`** | `(0 8 0)` transverse | inlet/outlet pressure | 1000 |
| `channelHartmannSpanwise` | **RAS Spalart–Allmaras, `Cv1 30`** | `(0 0 -8)` spanwise | inlet/outlet pressure | 1000 |
| `channelHartmannPeriodic` | laminar | `(0 8 0)` transverse | `meanVelocityForce` | 5000 |

`channelHartmann` is a **faithful reproduction of the upstream reference
case**: same mesh, same `SpalartAllmaras` closure with the MHD-modified
coefficient `Cv1 30` recommended by Dietiker & Hoffmann (2003), same fields
(`nut`, `nuTilda`), same `B = (0 8 0)`, same `endTime 1` with `dt 0.001`.
The `nutkWallFunction` requires the walls to be of patch type `wall`, so
`system/blockMeshDict` declares `lowerWall`/`upperWall` as `type wall;` —
matching the patch types of the upstream `polyMesh`.

The periodic case is deliberately kept laminar: its purpose is to exercise the
`pB` reference-level fix and the `fvOptions` coupling, and laminar keeps those
5000 steps under a minute.

Turbulence closures are standard OpenFOAM models used unchanged; the MHD
modification of Spalart–Allmaras is a **coefficient** set in the case
dictionary (`Cv1`), not a change to the model code. Spalart–Allmaras,
k-epsilon and k-omega SST all work.

## Verification

`./verify.sh --report` runs every case from scratch and rewrites
[`docs/verification.md`](docs/verification.md). All three cases reach `End`
without floating-point errors:

| case | exit | steps | reached `End` | max continuity error | max magnetic flux divergence error | max final p residual | s |
|---|---|---|---|---|---|---|---|
| channelHartmann | 0 | 1000 | 1 | 9.99335e-11 | 9.99608e-11 | 9.99636e-07 | 13 |
| channelHartmannSpanwise | 0 | 1000 | 1 | 9.99619e-12 | 0 | 9.99402e-07 | 10 |
| channelHartmannPeriodic | 0 | 5000 | 1 | 9.9951e-10 | 9.66615e-13 | 9.99939e-07 | 45 |

`cases/` is verification: each case has an objective metric. `examples/` is
different — those are integration/smoke tests with no exact solution, kept out
of `verify.sh` on purpose. `examples/swirlDuctCoarse` adds a **k-omega SST**
closure, the `swirlFlowRateInletVelocity` inlet and the polygonal duct wall.
It runs in ~30 s serial and is documented in its own
[`README`](examples/swirlDuctCoarse/README.md).

## Limitations

* The **full-induction formulation fails when the applied field is
  non-uniform**. The two Maxwell-stress terms are each of order
  $B_0^2/(\mu_0\rho)$ and their cancellation is exact only for a uniform
  $\mathbf{B}_0$; for a non-uniform applied field the residual spurious force
  dominates the physical one. For non-uniform fields use the quasi-static
  electric-potential solver
  [`mhdturbFoamQS`](https://github.com/wellfonseca/mhdturbFoamQS), which removes
  the cancellation by construction.
* Boundary conditions on the electric field or current density cannot be
  imposed; only the magnetic field is prescribed, as in `mhdFoam`.
* The applied field is prescribed through the initial condition of `B`, so a
  driven steady state with a changing field is not supported.

## Licence and attribution

Distributed under the **GNU General Public License v3.0** (see [`LICENSE`](LICENSE)),
which is mandatory given the upstream licences. Any redistribution must
preserve [`NOTICE`](NOTICE). See also [`CITATION.cff`](CITATION.cff).

If you use this package, please cite the upstream works listed in `NOTICE`
together with the archived release of this package.

## Citation

Archived on Zenodo:

* **Concept DOI** (all versions, always the latest):
  [10.5281/zenodo.23003542](https://doi.org/10.5281/zenodo.23003542)
* **Version 0.1.0**:
  [10.5281/zenodo.23003543](https://doi.org/10.5281/zenodo.23003543)

Cite the **concept DOI** to point at the latest version, or the **version DOI**
to refer to exactly the results archived as 0.1.0.

```
Fonseca, W. da S. (2026). mhdturbFoam-hardened: a corrected redistribution of
the FOSSEE mhdturbFoam turbulent MHD solver for OpenFOAM 6 (Version 0.1.0)
[Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23003543
```
