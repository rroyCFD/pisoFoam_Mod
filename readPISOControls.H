const dictionary& pisoDict = mesh.solutionDict().subDict("PISO");

// const int nOuterCorr =
//     pisoDict.lookupOrDefault<int>("nOuterCorrectors", 1);

// default minimum execution 1
const int nCorr =  pisoDict.lookupOrDefault<int>("nCorrectors", 1);
// Info<< "nCorrectors: " << nCorr << endl;

// default minimum execution 0
const int nNonOrthCorr = pisoDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
// Info << "nNonOrthogonalCorrectors: " << nNonOrthCorr << endl;

// default minimum execution 0
const int nPrimeIterations = pisoDict.lookupOrDefault<int>("nPrimeIterations", 0);
// Info<< "nPrimeIterations: " << nPrimeIterations << endl;

// const bool momentumPredictor =
//     pisoDict.lookupOrDefault("momentumPredictor", true);

// const bool transonic =
//     pisoDict.lookupOrDefault("transonic", false);

// const bool consistent =
//     pisoDict.lookupOrDefault("consistent", false);

// declared in createFields.H
// const bool pRefOn = pisoDict.lookupOrDefault("pRefOn", false);
// label  pRefCell =  pisoDict.lookupOrDefault<scalar>("pRefCell", 0);
// scalar pRefValue = pisoDict.lookupOrDefault<scalar>("pRefValue", 0);