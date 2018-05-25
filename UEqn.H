// step 2: Momentum predictor
fvVectorMatrix UEqn
(
    fvm::ddt(U)
  + fvm::div(phi, U)
  + turbulence->divDevReff(U)
  + fvm::SuSp(-fvc::div(phi),U)
);

//- Momentum solution adding the pressure gradient as a source term
solve
(
    UEqn == -fvc::grad(p)
);