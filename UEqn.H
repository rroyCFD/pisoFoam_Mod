// step 2: Solve implicitly the momentum equation given by Eq. (15.70)
// to obtain a new velocity field v*.
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