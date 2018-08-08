// Read whether the case is Taylor Green Vortex
const Switch TaylorGreenVortex
(
    transportProperties.lookupOrDefault("TaylorGreenVortex", false)
);

//
// const Switch init
// (
//     transportProperties.lookupOrDefault("init", false)
// );
//

// Read the dimension of the case
const word TaylorGreenDimension
(
    transportProperties.lookupOrDefault<word>("TaylorGreenDimension", "")
);

//  Read initial Velocity V_0
dimensionedScalar Uinit
(
    "Uinit",
    dimVelocity,
    transportProperties.lookupOrDefault<scalar>("Uinit", 1.0)
    // transportProperties.lookup("Uinit")
);

//  Read length scale
dimensionedScalar L
(
    "L",
    dimLength,
    transportProperties.lookupOrDefault<scalar>("L", 1.0)
    // transportProperties.lookup("L")
);

if(TaylorGreenVortex)
{
    Info << "TaylorGreenVortex: " << TaylorGreenVortex << endl;
    Info << "TaylorGreenDimension: " << TaylorGreenDimension << endl;
    Info << "Velocity Scale: " << Uinit.value() << endl;
    Info << "Length Scale: " << L.value() << endl;
}

// Declare variables -------------------------------------------- //
dimensionedScalar epsilon
(
    "epsilon",
    dimensionSet(0, 0, 0, 0, 0),
    0.0
);

dimensionedScalar Ek
(
    "Ek",
    dimensionSet(0, 0, 0, 0, 0),
    0.0
);

// characteristic time
dimensionedScalar tc
(
    "tc",
    dimensionSet(0, 0, 1, 0, 0),
    0.0
);
tc = L/Uinit;


// create output file
fileName outFilePath = runTime.path()/("TaylorGreenVortexProperties.dat_"+runTime.timeName());

if(Pstream::parRun() && Pstream::master())
{
    outFilePath = runTime.path()/"../"/("TaylorGreenVortexProperties.dat_"+runTime.timeName());
}
// Info << outFilePath << endl;

OFstream globalPropertiesFile(outFilePath);

// write to file
if(TaylorGreenVortex && Pstream::master())
{
   if(TaylorGreenDimension =="2D" )
    {
        globalPropertiesFile << "time" << "\t" << "Ek" << "\t" << "epsilon\t\t"
        << "L2(U)\t" << "L2(p)\t" << "Linf(U)\t" << "Linf(p)" << endl;
    }
    else if (TaylorGreenDimension =="3D" && Pstream::master())
    {
        globalPropertiesFile << "time" << "\t" << "Ek" << "\t" << "epsilon" << endl;
    }
    else {
        globalPropertiesFile << "TG Dimension mismatch: check for dimension in the transportproperties file" << endl;
    }
} else
{
    rm(outFilePath);
}


