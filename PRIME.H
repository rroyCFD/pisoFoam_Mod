// PRIME loop

//---------------------------------------------------------------------------------//
// Step 6: Using the latest available velocity and pressure fields,
// calculate the coefficients of the momentum equation and solve it explicitly.
fvVectorMatrix UprimeEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
  + fvm::SuSp(-fvc::div(phi),U)
);

rAU = 1.0/UprimeEqn.A();

// explicit calculation of velocity with mass conserved velocity, phi and updated pressure
U = rAU*(UprimeEqn.H() - fvc::grad(p));  // eq. 15.176


//------------------------------------------------------------------------------------------------//
// Step 7: Update the mass flow rate at the cell faces using the Rhie-Chow interpolation technique.
rAUf = linearInterpolate(rAU);

gradpDiff_f =   -(linearInterpolate(fvc::grad(p)) & ed)*ed
                + (faceGradient.snGrad(p))*ed;

phi = (fvc::interpolate(U) & mesh.Sf())
    - (rAUf*gradpDiff_f & mesh.Sf())
    + fvc::ddtPhiCorr(rAU, U, phi);


//------------------------------------------------------------------------------------------------//
// Step 8: Using the new mass flow rates, assemble the pressure correction equation (Eq. 15.183)
// and solve it to obtain a pressure correction field.
pp.internalField() = scalar(0.0);
pp.correctBoundaryConditions();

// Non-orthogonal pressure corrector loop
for (int nonOrth=0; nonOrth <= nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix ppPrimeEqn
    (
        fvm::laplacian(rAUf, pp,"laplacian(pDiff,pp)")
        ==
        fvc::div(phi)
    );

    // set pp at refCell such that after solving pp,subsequent update of pressure (p),
    // pressure value at pRefCell is within solver tolerance of prescribed pRefValue.

    // ppPrimeEqn.setReference(pRefCell, pRefValue);
    ppPrimeEqn.setReference(pRefCell, (pRefValue - getRefCellValue(p, pRefCell)));

    // final pressure coorection solved upto the tightest convergence tolerance
    if(corr == nCorr && nPrime == nPrimeIterations && nonOrth == nNonOrthCorr)
    {
        ppPrimeEqn.solve(mesh.solver("ppFinal"));
    }
    else
    {
        ppPrimeEqn.solve(mesh.solver("pp"));
    }
    // Info << "solved ppPrimeEqn" << endl;


    //------------------------------------------------------------------------------------------------//
    // Step 9: Update pressure, velocity, and mass flow rate fields using expressions similar
    // to the ones given in Eq. (15.101).
    // update fields at the end of non-orthogonal correction loop
    if(nonOrth == nNonOrthCorr)
    {
        p += pp;
        p.correctBoundaryConditions();

        U -= rAU*fvc::grad(pp);
        U.correctBoundaryConditions();

        phi -= ppPrimeEqn.flux();
    }

}// end of non-orthogonal correction loop

// >>> Go to step 6 and repeat based on the desired number of PRIME steps.
