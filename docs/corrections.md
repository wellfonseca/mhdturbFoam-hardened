# Corrections applied to the upstream `mhdturbFoam`

`mhdturbFoam-hardened` differs from the FOSSEE `mhdturbFoam` (2019) only in the
three numerical corrections below and in the packaging. Each correction is
accompanied by the evidence that motivated it.

The reference for what the code *should* do is OpenFOAM 6's `mhdFoam`
(`$WM_PROJECT_DIR/applications/solvers/electromagnetics/mhdFoam`), from which
`mhdturbFoam` descends.

---

## 1. Restored outer PISO corrector loop

### Upstream code

```cpp
// --- PISO loop
 {
    volScalarField rAU(1.0/UEqn.A());
    ...
    pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));
    ...
}
```

The `while (piso.correct())` of `mhdFoam` was replaced by a plain block.

### Consequences

* `rAU`, `HbyA` and `phiHbyA` are built **once per time step** instead of once
  per corrector, so the pressure–velocity coupling is effectively first order.
* `piso.correct()` is never called, therefore `corrPISO_` stays 0 and
  `piso.finalPISOIter()` is always false. `piso.finalInnerIter()`
  (= `finalNonOrthogonalIter() && finalPISOIter()`) is then never true, so the
  `pFinal` dictionary is unreachable and every pressure solve stops at the
  `relTol` of the `p` dictionary.

### Evidence

On the `channelHartmann` case of this package, the upstream code left the final
pressure residual between `1e-3` and `1e-2` (the `p` dictionary `relTol 0.05`,
reduced by the number of iterations), while the hardened code reaches the
`pFinal` tolerance `1e-6`:

```
max final pressure residual, hardened  : 9.99771e-07
```

On the periodic channel of the companion investigation the same defect made the
continuity error grow from `5e-4` to `1471`, and the run died with a floating
point exception. This is recorded verbatim in the development copy of the
solver:

> "O laco externo (piso.correct(), nCorrectors em fvSolution) e' indispensavel:
> rAU, HbyA e phiHbyA tem de ser recalculados em TODOS os corretores. Com um
> unico corretor o acoplamento pressao-velocidade fica de primeira ordem e a
> solucao divergiu mesmo com Courant << 1 (verificado no caso channelHartmann:
> sem o laco o erro de continuidade cresceu de 5e-4 para 1471 e o run morreu
> com FPE)."

### Fix

```cpp
while (piso.correct())
{
    volScalarField rAU(1.0/UEqn.A());
    ...
}
```

With the loop restored, `piso.finalInnerIter()` is true exactly in the last
non-orthogonal correction of the last corrector, which is when `pFinal` must be
used — the same structure as `mhdFoam` and `pisoFoam`.

---

## 2. Reference level for the singular `pB` Poisson problem

### Upstream code

```cpp
while (bpiso.correctNonOrthogonal())
{
    fvScalarField pBEqn
    (
        fvm::laplacian(rABf, pB) == fvc::div(phiB)
    );

    pBEqn.solve(mesh.solver(pB.select(bpiso.finalInnerIter())));
    ...
}
```

There is no `setReference` for `pB`.

### Consequences

When every `pB` boundary is `zeroGradient` or `cyclic` — a closed domain, or a
domain that is fully periodic in the resolved directions, such as a
streamwise-periodic channel with `empty` front and back — the discrete operator
is a purely Neumann Laplacian and therefore singular: the solution is defined
only up to a constant. `pB` then drifts and the run dies with a floating point
exception. Open domains that fix `pB` with `fixedValue` at an outlet are not
affected, which is why the defect went unnoticed upstream.

### Fix

```cpp
const dictionary& BPISOdict = mesh.solutionDict().subDict("BPISO");
const label  pBRefCell  = BPISOdict.lookupOrDefault<label>("pBRefCell", 0);
const scalar pBRefValue = BPISOdict.lookupOrDefault<scalar>("pBRefValue", 0.0);
...
pBEqn.setReference(pBRefCell, pBRefValue);
```

The entries are optional and default to cell 0 / value 0, so existing open
cases keep working unchanged. A case that needs them sets, in
`system/fvSolution`:

```
BPISO
{
    nCorrectors     3;
    pBRefCell       0;
    pBRefValue      0;
}
```

The `channelHartmannPeriodic` case in this package requires them and runs 5000
time steps to completion with a maximum magnetic flux divergence error of
`9.8e-12`.

---

## 3. `fvOptions` coupling in the momentum equation

### Upstream code

```cpp
fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
  - fvc::div(phiB, 2.0*DBU*B)
  + fvc::grad(DBU*magSqr(B))
);

UEqn.relax();
```

`fvOptions` is included in `createFields.H` but never used.

### Consequences

Any `fvOptions` source configured by the user is silently ignored. In
particular `fv::meanVelocityForce`, the standard way of driving a periodic
channel, has no effect, and because `fvOptions.constrain(UEqn)` is never called
its `rAPtr_` is never initialised from `1/UEqn.A()` — so even if the source were
added the force would be wrong.

### Fix

```cpp
fvVectorMatrix UEqn
(
    ...
 ==
    fvOptions(U)
);

UEqn.relax();

// initialises rAPtr_ and accumulates gradP0_ in fv::meanVelocityForce
fvOptions.constrain(UEqn);
...
fvOptions.correct(U);
```

The `channelHartmannPeriodic` case drives a streamwise-periodic channel with
`meanVelocityForce` at `Ubar = (1 0 0)` and holds it for 5000 steps at a
constant pressure gradient.

---

## What is *not* changed

* The Maxwell-stress Lorentz term
  `- fvc::div(phiB, 2.0*DBU*B) + fvc::grad(DBU*magSqr(B))` — inherited verbatim
  from `mhdFoam`.
* The B-PISO loop structure, `magneticFieldErr.H`, `createPhiB.H`,
  `createControl.H`, `readBPISOControls.H`.
* The `createFields.H` of the upstream derivative (transport, turbulence,
  `fvOptions`).
* No magnetic modification is applied to the turbulence closures. The upstream
  `README.md` recommends a modified Spalart–Allmaras model for MHD; that
  recommendation is **not** implemented here, and the standard OpenFOAM
  closures are used unmodified.
