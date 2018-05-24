/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Application
    pisoFoam_Mod

Description
    Transient solver for incompressible flow.
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulenceModel.H"
#include "orthogonalSnGrad.H"

#include "IFstream.H"
#include "OFstream.H"
#include "IOmanip.H" // for input/ouput format control
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readTransportProperties.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"
    //#include "readTimeControls.H"

    #include "TaylorGreenFiles/readAndDeclareVariables.H"
    #include "TaylorGreenFiles/createErrorFields.H"

    #include "TaylorGreenFiles/initialize.H"
    #include "TaylorGreenFiles/errorNorm.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    //Rhie Chow interpolation stuff
    const surfaceVectorField ed = mesh.delta()()/mag(mesh.delta()());
    Foam::fv::orthogonalSnGrad<scalar> faceGradient(mesh);

    surfaceVectorField gradp_avg_f = linearInterpolate(fvc::grad(p));
    surfaceVectorField gradpDiff_f
        = -(gradp_avg_f & ed)*ed + (faceGradient.snGrad(p))*ed;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "CourantNo.H"

        // --- PISO loop
        for (int corr=0; corr<nCorr; corr++)
        {
            // step 2: Momentum predictor
            fvVectorMatrix UEqn
            (
                fvm::ddt(U)
              + fvm::div(phi, U)
              + turbulence->divDevReff(U)
            );

            //- Momentum solution adding the pressure gradient as a source term
            solve
            (
                UEqn == -fvc::grad(p)
            );

            // step 3:
            volScalarField rAU(1.0/UEqn.A());
            surfaceScalarField rAUf("rAUf",linearInterpolate(rAU));

            // volScalarField velDiff("velDiff", mag(U - (rAU*(UEqn.H() -fvc::grad(p)))));
            // Info << max(velDiff) << endl;

            surfaceScalarField rAUTf("rAUTf",rAUf/runTime.deltaT());

            gradp_avg_f = linearInterpolate(fvc::grad(p));
            gradpDiff_f = -(gradp_avg_f & ed)*ed + (faceGradient.snGrad(p))*ed;

            surfaceVectorField U_old_f = linearInterpolate(U.oldTime());
            surfaceScalarField phi_old = phi.oldTime();

            phi = (fvc::interpolate(U) & mesh.Sf())
                - (rAUf*gradpDiff_f & mesh.Sf())
                // + rAUTf*(phi_old - (U_old_f& mesh.Sf()));
                + fvc::ddtPhiCorr(rAU, U, phi);

            // surfaceScalarField temporalPhiCorr(
            //                                    "temporalPhiCorr",
            //                                    (
            //                                         fvc::ddtPhiCorr(rAU, U, phi)
            //                                         - rAUTf*(phi_old - (U_old_f& mesh.Sf())))
            //                                    );

            // Info << max(temporalPhiCorr) << endl;

            // Step 4:
            //- resetting pressure correction
            pp.internalField() = scalar(0.0);
            pp.correctBoundaryConditions();

            // Non-orthogonal pressure corrector loop
            // for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
           // Pressure corrector
            fvScalarMatrix ppEqn
            (
                - fvm::laplacian(rAUf, pp,"laplacian(pDiff,pp)")
                + fvc::div(phi)
            );

            ppEqn.setReference(pRefCell, pRefValue);
            ppEqn.solve();

            phi += ppEqn.flux();

            p += pp;
            p.correctBoundaryConditions();

            U -= rAU*fvc::grad(pp);
            U.correctBoundaryConditions();

            // PRIME loop
            for (int nPrime=0; nPrime < nPrimeIterations; nPrime++)
            {
                // Info << "nPrime: " << nPrime+1 << endl;
                fvVectorMatrix UprimeEqn
                (
                    fvm::ddt(U)
                  + fvm::div(phi, U)
                  + turbulence->divDevReff(U)
                );

                rAU = 1.0/UprimeEqn.A();

                U = rAU*(UprimeEqn.H() -fvc::grad(p));

                gradp_avg_f = linearInterpolate(fvc::grad(p));
                gradpDiff_f = -(gradp_avg_f & ed)*ed + (faceGradient.snGrad(p))*ed;

                rAUf = linearInterpolate(rAU);
                rAUTf = rAUf/runTime.deltaT();

                gradp_avg_f = linearInterpolate(fvc::grad(p));
                gradpDiff_f = -(gradp_avg_f & ed)*ed + (faceGradient.snGrad(p))*ed;

                phi = (fvc::interpolate(U) & mesh.Sf())
                    - (rAUf*gradpDiff_f & mesh.Sf())
                    // + rAUTf*(phi_old - (U_old_f& mesh.Sf()));
                    + fvc::ddtPhiCorr(rAU, U, phi);

                pp.internalField() = scalar(0.0);
                pp.correctBoundaryConditions();

                fvScalarMatrix ppPrimeEqn
                (
                - fvm::laplacian(rAUf, pp,"laplacian(pDiff,pp)")
                + fvc::div(phi)
                );

                ppPrimeEqn.setReference(pRefCell, pRefValue);

                ppPrimeEqn.solve(mesh.solver("pp"));

                phi += ppPrimeEqn.flux();

                p += pp;
                p.correctBoundaryConditions();

                U -= rAU*fvc::grad(pp);
                U.correctBoundaryConditions();
            }
            // end of PRIME loop

        }// end of corrector loop

        #include "continuityErrs.H"

        #include "TaylorGreenFiles/errorNorm.H"
        #include "TaylorGreenFiles/globalProperties.H"

        turbulence->correct();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    } // end of runTime loop

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
