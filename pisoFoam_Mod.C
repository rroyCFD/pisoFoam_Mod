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
    //#include "readTimeControls.H" // for fixed courant number

    #include "TaylorGreenFiles/readAndDeclareVariables.H"
    #include "TaylorGreenFiles/createErrorFields.H"
    #include "TaylorGreenFiles/initialize.H"
    #include "TaylorGreenFiles/errorNorm.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    //Rhie Chow interpolation stuff
    const surfaceVectorField ed = mesh.delta()()/mag(mesh.delta()());
    Foam::fv::orthogonalSnGrad<scalar> faceGradient(mesh);

    //surfaceVectorField gradp_avg_f = linearInterpolate(fvc::grad(p));
    surfaceVectorField gradpDiff_f
        =
        -(linearInterpolate(fvc::grad(p)) & ed)*ed
        + (faceGradient.snGrad(p))*ed;


    // run-time coverge and advance loop
    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readPISOControls.H"
        #include "CourantNo.H"
        //#include "setDeltaT.H" // for fixed courant number, variable timeStep

        // --- PISO loop
        for (int corr=0; corr<nCorr; corr++)
        {
            #include "UEqn.H"
            #include "ppEqn.H"

            // PRIME loop
            for (int nPrime=0; nPrime < nPrimeIterations; nPrime++)
            {
                #include "PRIME.H"
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