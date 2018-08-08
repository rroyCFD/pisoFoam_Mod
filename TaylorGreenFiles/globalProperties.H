if(TaylorGreenVortex)
{
    // Dissipation rate normalized
    epsilon = 2*nu*magSqr(symm(fvc::grad(U)))().weightedAverage(mesh.V())/(Uinit*Uinit*Uinit/L);

    // Kinetic enery normalized
    Ek = 0.5*magSqr(U)().weightedAverage(mesh.V())/ (Uinit*Uinit);

    if(Pstream::master())
    {
        if(TaylorGreenDimension =="2D")
        {
            globalPropertiesFile
            << runTime.value()/tc.value() << "\t" << Ek.value() << "\t" << epsilon.value()
            << "\t" << UL2err.value() << "\t" << pL2err.value()
            << "\t" << ULinfErr.value() << "\t" << pLinfErr.value()<< endl;
        }
        else if (TaylorGreenDimension =="3D")
        {
            globalPropertiesFile << runTime.value()/tc.value() << "\t" << Ek.value() << "\t" << epsilon.value()<< endl;
        }
    }
}

