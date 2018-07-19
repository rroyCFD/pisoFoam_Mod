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

    Algorithom is implemented followed by the outline in book:
        The Finite Volume Method in Computational Fluid Dynamics
        Section 15.7 (SIMPLE family of algorithoms)

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

    // #include "TaylorGreenFiles/readAndDeclareVariables.H"
    // #include "TaylorGreenFiles/createErrorFields.H"
    // #include "TaylorGreenFiles/initialize.H"
    // #include "TaylorGreenFiles/errorNorm.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    //Rhie Chow interpolation stuff
    const surfaceVectorField ed = mesh.delta()()/mag(mesh.delta()());
    Foam::fv::orthogonalSnGrad<scalar> faceGradient(mesh);

    surfaceVectorField gradpDiff_f
        =
        -(linearInterpolate(fvc::grad(p)) & ed)*ed
        + (faceGradient.snGrad(p))*ed;

    // update Fields
    Info <<"\nUpdating boundary fields..." << endl;
    p.correctBoundaryConditions();
    U.correctBoundaryConditions();

    Info<< "\nStarting time loop\n" << endl;
    // run-time coverge and advance loop
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "readTimeControls.H" // for fixed courant number
        #include "CourantNo.H"
        #include "setDeltaT.H" // for fixed courant number, variable timeStep

        // store old values for temporal discretization,
        // & temporal correction to phi in Rhie-Chow flux calculation.
        U.storeOldTime();
        phi.storeOldTime();

        // --- PISO loop
        for (int corr=1; corr<=nCorr; corr++)
        {
            #include "UEqn.H"


                #include "ppEqn.H"

                // PRIME loop
                for (int nPrime=1; nPrime <= nPrimeIterations; nPrime++)
                {
                    #include "PRIME.H"
                }
                // end of PRIME loop

        }// end of corrector loop

        #include "continuityErrs.H"

        // #include "TaylorGreenFiles/errorNorm.H"
        // #include "TaylorGreenFiles/globalProperties.H"

        turbulence->correct();

        //update nuEff
        nuEff = turbulence->nuEff();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    } // end of runTime loop

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
