/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    mhdturbFoam

Group
    grpElectroMagneticsSolvers

Description
    Solver for magnetohydrodynamics (MHD): incompressible, turbulent flow of a
    conducting fluid under the influence of a magnetic field.

    An applied magnetic field H acts as a driving force,
    at present boundary conditions cannot be set via the
    electric field E or current density J. The fluid viscosity nu,
    conductivity sigma and permeability mu are read in as uniform
    constants.

    A fictitous magnetic flux pressure pH is introduced in order to
    compensate for discretisation errors and create a magnetic face flux
    field which is divergence free as required by Maxwell's equations.

    However, in this formulation discretisation error prevents the normal
    stresses in UB from cancelling with those from BU, but it is unknown
    whether this is a serious error.  A correction could be introduced
    whereby the normal stresses in the discretised BU term are replaced
    by those from the UB term, but this would violate the boundedness
    constraint presently observed in the present numerics which
    guarantees div(U) and div(H) are zero.

    This is the "hardened" variant of the FOSSEE mhdturbFoam derivative of
    OpenFOAM's mhdFoam.  See docs/corrections.md and the NOTICE file.  Three
    numerical defects of the upstream derivative are fixed here:

      1. the outer PISO corrector loop (while (piso.correct())) is restored,
         so that rAU, HbyA and phiHbyA are recomputed at every corrector and
         the pFinal solver dictionary is actually reached;
      2. the singular Poisson problem for the magnetic flux pressure pB is
         given a reference cell/value through the BPISO dictionary, which is
         required in closed or fully periodic domains;
      3. fvOptions are coupled into the momentum equation (source term and
         constrain()), which is required by meanVelocityForce and by any other
         momentum source used to drive the flow.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pisoControl.H"
#include "fvOptions.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for magnetohydrodynamics (MHD):"
        " incompressible, turbulent flow of a conducting fluid"
        " under the influence of a magnetic field."
    );

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"

    // Reference level for the pB Poisson problem (closed or fully periodic
    // domains).  The entry is optional and defaults to cell 0 / value 0, so
    // ordinary open domains with a fixed-value pB boundary are unaffected.
    const dictionary& BPISOdict = mesh.solutionDict().subDict("BPISO");
    const label pBRefCell = BPISOdict.lookupOrDefault<label>("pBRefCell", 0);
    const scalar pBRefValue = BPISOdict.lookupOrDefault<scalar>("pBRefValue", 0.0);

    #include "createFields.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        {
            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
              + fvm::div(phi, U)
              + turbulence->divDevReff(U)
              - fvc::div(phiB, 2.0*DBU*B)
              + fvc::grad(DBU*magSqr(B))
             ==
                fvOptions(U)
            );

            UEqn.relax();

            // Required by fv::meanVelocityForce, which initialises rAPtr_ from
            // 1/UEqn.A() and accumulates gradP0_ inside constrain(); without
            // this call rAPtr_ stays empty and the run aborts.
            fvOptions.constrain(UEqn);

            if (piso.momentumPredictor())
            {
                solve(UEqn == -fvc::grad(p));
            }

            // --- PISO loop
            //
            // The outer corrector loop (piso.correct(), nCorrectors in
            // fvSolution) is essential: rAU, HbyA and phiHbyA must be
            // recomputed at EVERY corrector.  With a single corrector the
            // pressure-velocity coupling is only first order and the solution
            // diverges even at Courant numbers well below one.  This was
            // verified on the channelHartmann case: without the loop the
            // continuity error grew from 5e-4 to 1471 and the run died with a
            // floating point exception.
            //
            // Restoring the loop also restores the correct meaning of
            // piso.finalInnerIter(): it is true only in the last non-orthogonal
            // correction of the last PISO corrector, which is exactly when the
            // pFinal dictionary must be used.
            while (piso.correct())
            {
                volScalarField rAU(1.0/UEqn.A());
                surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
                volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
                surfaceScalarField phiHbyA
                (
                    "phiHbyA",
                    fvc::flux(HbyA)
                  + rAUf*fvc::ddtCorr(U, phi)
                );

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(p, U, phiHbyA, rAUf);

                while (piso.correctNonOrthogonal())
                {
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
                    );

                    pEqn.setReference(pRefCell, pRefValue);

                    pEqn.solve(mesh.solver(p.select(piso.finalInnerIter())));

                    if (piso.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                #include "continuityErrs.H"

                U = HbyA - rAU*fvc::grad(p);
                U.correctBoundaryConditions();
            }

            fvOptions.correct(U);
        }


        // --- B-PISO loop
        while (bpiso.correct())
        {
            fvVectorMatrix BEqn
            (
                fvm::ddt(B)
              + fvm::div(phi, B)
              - fvc::div(phiB, U)
              - fvm::laplacian(DB, B)
            );

            BEqn.solve();

            volScalarField rAB(1.0/BEqn.A());
            surfaceScalarField rABf("rABf", fvc::interpolate(rAB));

            phiB = fvc::flux(B);

            while (bpiso.correctNonOrthogonal())
            {
                fvScalarMatrix pBEqn
                (
                    fvm::laplacian(rABf, pB) == fvc::div(phiB)
                );

                // In a closed or fully periodic domain every pB boundary is
                // zeroGradient or cyclic, so the operator is a purely Neumann
                // Poisson problem and therefore singular: the level must be
                // fixed (up to a constant), exactly as for p.  With an open
                // outlet carrying fixedValue pB this is a harmless no-op.
                pBEqn.setReference(pBRefCell, pBRefValue);

                pBEqn.solve(mesh.solver(pB.select(bpiso.finalInnerIter())));

                if (bpiso.finalNonOrthogonalIter())
                {
                    phiB -= pBEqn.flux();
                }
            }

            #include "magneticFieldErr.H"
        }
        laminarTransport.correct();
        turbulence->correct();

        runTime.write();

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
