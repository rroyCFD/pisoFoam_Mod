if(TaylorGreenVortex && TaylorGreenDimension =="2D")
{
    // Velocity error norms--------------------------------------------------------------//
    Info << "Calculating Uerror" << endl;
    Uerror = Ustart*Foam::exp(-2.0*nu.value()*runTime.value()) - U;

    ULinfErr = max(mag(Uerror));   // infinity norm
    UL2err = sqrt(magSqr(Uerror)().weightedAverage(mesh.V())); // 2nd norm

    // Pressure error norms ----------------------------------------------------------//
    Info << "Calculating perror" << endl;

    if(pRefOn) {

        // Reference analytical pressure at pRefCell
        const scalar xPRefCell = mesh.C()[pRefCell].x();
        const scalar yPRefCell = mesh.C()[pRefCell].y();

        pRefValueAna  = sqr(Uinit)* 0.25*(
                          Foam::cos(2.*xPRefCell/L.value())
                        + Foam::cos(2.*yPRefCell/L.value()) )
                        * Foam::exp(-4.*nu.value()*runTime.value());

        Info << "Analytical value at pRefCell: "<< pRefValueAna.value() << endl;
        perror = pstart*Foam::exp(-4.0*nu.value()*runTime.value()) - pRefValueAna - p;
    } else {
        perror = pstart*Foam::exp(-4.0*nu.value()*runTime.value()) - p;
    }


    pLinfErr = max(mag(perror));    // infinity norm
    pL2err = sqrt(magSqr(perror)().weightedAverage(mesh.V()));  // 2nd norm

    Info << "L2   error norm | Velocity: " << UL2err.value()   << "\tPressure: " << pL2err.value()   << endl;
    Info << "Linf error norm | Velocity: " << ULinfErr.value() << "\tPressure: " << pLinfErr.value() << endl;

    Info << "Pressure decay coeffs: " << setprecision(10) << Foam::exp(-4.0*nu.value()*runTime.value()) << endl;
    Info << "Velocity decay coeffs: " << setprecision(10) << Foam::exp(-2.0*nu.value()*runTime.value()) << endl;

    // write the error induced during simulation
    if(runTime.timeName() == "0" || runTime.write())
    {
        Uerror.write();
        perror.write();
    }
}