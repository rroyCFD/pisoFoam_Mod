//------------------------------------------------------------------------------------------------//
// Step 3: Update the mass flow rate at the cell faces using the Rhie-Chow interpolation
// technique (Eq. 15.100) to obtain a momentum satisfying mass flow rate field m^∗f.
volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("rAUf",linearInterpolate(rAU));

gradpDiff_f =   -(linearInterpolate(fvc::grad(p)) & ed) * ed
                + (faceGradient.snGrad(p)) * ed;

phi = (fvc::interpolate(U) & mesh.Sf())
    - (rAUf*gradpDiff_f & mesh.Sf())
    + fvc::ddtPhiCorr(rAU, U, phi);

// Step 4: Using the new mass flow rates, assemble the pressure correction
// equation (Eq. 15.98) and solve it to obtain a pressure correction field pp.
//- resetting pressure correction
pp.internalField() = scalar(0.0);
pp.correctBoundaryConditions();

// Non-orthogonal pressure corrector loop
for (int nonOrth=0; nonOrth <= nNonOrthCorr; nonOrth++)
{
    // Pressure corrector
    fvScalarMatrix ppEqn
    (
        fvm::laplacian(rAUf, pp,"laplacian(pDiff,pp)")
        ==
        fvc::div(phi)
    );

    ppEqn.setReference(pRefCell, pRefValue);
    ppEqn.solve();
    // Info << "solved ppEqn" << endl;

    if(nonOrth == nNonOrthCorr)
    {
        // Step 5: Update the pressure and velocity fields at the cell centroids
        // and the mass flow rate at the cell faces to obtain continuity-satisfying fields using Eq. (15.101).
        // update fields at the end of non-orthogonal correction loop
        p += pp;
        p.correctBoundaryConditions();

        U -= rAU*fvc::grad(pp);
        U.correctBoundaryConditions();

        phi -= ppEqn.flux();
    }

} // end of non-orthogonal correction loop
