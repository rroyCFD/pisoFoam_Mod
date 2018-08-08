// Rhie-Chow interpolation stuff
    const surfaceVectorField ed = mesh.delta()()/mag(mesh.delta()());
    Foam::fv::orthogonalSnGrad<scalar> faceGradient(mesh);

    surfaceVectorField gradpDiff_f
        =
        - (linearInterpolate(fvc::grad(p)) & ed) * ed
        + (faceGradient.snGrad(p)) * ed;
