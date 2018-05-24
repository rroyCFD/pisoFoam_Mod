    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    const volScalarField::GeometricBoundaryField& pbf = p.boundaryField();
    wordList pbt = pbf.types();

    volScalarField pp
    (
        IOobject
        (

            "pp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", p.dimensions(), 0.0),
        pbt
    );


    // Set pp boundary values
    forAll(pp.boundaryField(), patchi)
    {
        if (isType<fixedValueFvPatchScalarField>(pp.boundaryField()[patchi]))
        {
           fixedValueFvPatchScalarField& ppbound =
                           refCast<fixedValueFvPatchScalarField>(pp.boundaryField()[patchi]);

           ppbound == scalarField(ppbound.size(),0.0);
        }

    }

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "createPhi.H"

    const dictionary& pisoDict = mesh.solutionDict().subDict("PISO");
    const bool pRefOn = pisoDict.lookupOrDefault("pRefOn", false);
    label  pRefCell =  pisoDict.lookupOrDefault<scalar>("pRefCell", 0);
    scalar pRefValue = pisoDict.lookupOrDefault<scalar>("pRefValue", 0);

    if(pRefOn)
    {
        Info << "pRefOn: " << pRefOn
            <<"  pRefCell: "<< pRefCell
            << "  pRefValue: "<< pRefValue
            <<endl;
    }

    setRefCell(p, pisoDict, pRefCell, pRefValue);

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );
