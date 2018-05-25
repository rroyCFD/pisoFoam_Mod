// PRIME loop

//---------------------------------------------------------------------------------//
// Step 6: Using the latest available velocity and pressure fields,
// calculate the coefficients of the momentum equation and solve it explicitly.

// Info << "nPrime: " << nPrime+1 << endl;
fvVectorMatrix UprimeEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
  + fvm::SuSp(-fvc::div(phi),U)
);

rAU = 1.0/UprimeEqn.A();

// explicit calculation of velocity with mass conserved velocity, phi and updated pressure
U = rAU*(UprimeEqn.H() -fvc::grad(p));


//------------------------------------------------------------------------------------------------//
// Step 7: Update the mass flow rate at the cell faces using the Rhie-Chow interpolation technique.
rAUf = linearInterpolate(rAU);
rAUTf = rAUf/runTime.deltaT();

gradpDiff_f =   -(linearInterpolate(fvc::grad(p)) & ed)*ed
                + (faceGradient.snGrad(p))*ed;

phi = (fvc::interpolate(U) & mesh.Sf())
    - (rAUf*gradpDiff_f & mesh.Sf())
    // + rAUTf*(phi_old - (U_old_f& mesh.Sf())); // Euler
    + fvc::ddtPhiCorr(rAU, U, phi);


//------------------------------------------------------------------------------------------------//
// Step 8: Using the new mass flow rates, assemble the pressure correction equation (Eq. 15.183)
// and solve it to obtain a pressure correction field.
pp.internalField() = scalar(0.0);
pp.correctBoundaryConditions();

fvScalarMatrix ppPrimeEqn
(
- fvm::laplacian(rAUf, pp,"laplacian(pDiff,pp)")
+ fvc::div(phi)
);

ppPrimeEqn.setReference(pRefCell, pRefValue);

ppPrimeEqn.solve(mesh.solver("pp"));


//------------------------------------------------------------------------------------------------//
// Step 9: Update pressure, velocity, and mass flow rate fields using expressions similar
// to the ones given in Eq. (15.101).
phi += ppPrimeEqn.flux();

p += pp;
p.correctBoundaryConditions();

U -= rAU*fvc::grad(pp);
U.correctBoundaryConditions();

// >>>>>
// Go to step 6 and repeat based on the desired number of PRIME steps.
