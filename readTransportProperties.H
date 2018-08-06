Info<< "Reading transportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED, //MUST_READ,
        IOobject::NO_WRITE,
        false
    )
);

dimensionedScalar nu
(
    transportProperties.lookup("nu")
);
Info << "Viscosity: " << nu.value() << endl;

// ************************************************************************* //
