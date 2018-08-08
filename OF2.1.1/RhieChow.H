// Define ed
surfaceVectorField ed
(
    IOobject
    (
        "ed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("ed", dimLength, vector::zero)
);

{
    const volVectorField& C = mesh.C();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    forAll(owner, facei)
    {
        ed[facei] = C[neighbour[facei]] - C[owner[facei]];
    }

    forAll(ed.boundaryField(), patchi)
    {
        ed.boundaryField()[patchi] = mesh.boundary()[patchi].delta();
    }

    ed /= mag(ed); // normalization
}


// Rhie-Chow interpolation stuff
Foam::fv::orthogonalSnGrad<scalar> faceGradient(mesh);

surfaceVectorField gradpDiff_f
    =
    - (linearInterpolate(fvc::grad(p)) & ed) * ed
    + (faceGradient.snGrad(p)) * ed;
