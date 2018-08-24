    Info<< "\nReading field p" << endl;
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

    const dictionary& pisoDict = mesh.solutionDict().subDict("PISO");
    const bool   pRefOn    = pisoDict.lookupOrDefault<bool>("pRefOn", false);
    label  pRefCell  = pisoDict.lookupOrDefault<scalar>("pRefCell", 0);
    scalar pRefValue = pisoDict.lookupOrDefault<scalar>("pRefValue", 0);

    if(pRefOn)
    {
        Info <<"pRefOn: " << pRefOn << " pRefCell: "<< pRefCell << " pRefValue: "<< pRefValue  << endl;

        if(mag(pRefValue - p[pRefCell]) > VSMALL)
        {
            if(Pstream::master())
            {
                Pout << "Un-corrected pressure at pRefCell: " << p[pRefCell] << endl;
            }

            p += (pRefValue -p[pRefCell])*dimensionedScalar("unitPressure", p.dimensions(), 1.0);

            if(Pstream::master())
            {
                Pout << "Corrected    pressure at pRefCell: " << p[pRefCell] << endl;
            }
        } // end of pressurecorrection loop

        Info << "writing corrected pressure field" << endl;
        p.correctBoundaryConditions();
        p.write();
    }

     // If the field fieldRef needs referencing find the reference cell nearest
    setRefCell(p, pisoDict, pRefCell, pRefValue);


    //-----------------------------------------------------------------------------------------//
    Info << "\nCreating pressure-coorection field pp" << endl;
    const volScalarField::GeometricBoundaryField& pbf = p.boundaryField();
    wordList pbt = pbf.types();

    // prepare the boundary condition types for pressure-correction field
    wordList ppbt = pbt;
    forAll(ppbt, b)
    {
        if(ppbt[b] == "timeVaryingMappedFixedValue")
        {
            ppbt[b] = "fixedValue";
        }
    }

    volScalarField pp
    (
        IOobject
        (

            "pp",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT, // NO_READ,
            IOobject::AUTO_WRITE // NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", p.dimensions(), 0.0),
        ppbt
    );
    Info << "Pressure correction fields boundary types:\n" << pp.boundaryField().types() << endl;


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

    //-----------------------------------------------------------------------------------------//
    Info<< "\nReading field U" << endl;
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

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );


    //-----------------------------------------------------------------------------------------//
    Info<< "\nCreating field nuEff" << endl;
    volScalarField nuEff
    (
        IOobject
        (
            "nuEff",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        //turbulence->nuEff() // can't find in lookupObjectRegistry, in BC implementation
        mesh,
        dimensionedScalar("nuEff",dimensionSet(0, 2, -1, 0, 0, 0, 0), 0.0)
    );
    nuEff = turbulence->nuEff();
    Info << "max(nuEff): " << max(nuEff) << endl;