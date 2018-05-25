// step 3:
volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rAUf("rAUf",linearInterpolate(rAU));

// volScalarField velDiff("velDiff", mag(U - (rAU*(UEqn.H() -fvc::grad(p)))));
// Info << max(velDiff) << endl;

surfaceScalarField rAUTf("rAUTf",rAUf/runTime.deltaT());

gradpDiff_f =   -(linearInterpolate(fvc::grad(p)) & ed)*ed
                + (faceGradient.snGrad(p))*ed;

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