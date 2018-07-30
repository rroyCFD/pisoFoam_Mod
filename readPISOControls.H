const dictionary& pisoDict = mesh.solutionDict().subDict("PISO");

// default value 1, minimum execution 1
const int nCorr =  pisoDict.lookupOrDefault<int>("nCorrectors", 1);

// default value 1; minimum execution 1
const int nPrimeIterations = pisoDict.lookupOrDefault<int>("nPrimeIterations", 1);

// default value 0; minimum execution 1
const int nNonOrthCorr = pisoDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

if(runTime.value() == (runTime.deltaTValue() + runTime.startTime().value()) )
{
    Info << "nCorrectors: " << nCorr << endl;
    Info << "nPrimeIterations: " << nPrimeIterations << endl;
    Info << "nNonOrthogonalCorrectors: " << nNonOrthCorr  << nl << endl;
}


// const int nOuterCorr =
//     pisoDict.lookupOrDefault<int>("nOuterCorrectors", 1);

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