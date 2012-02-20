
    volScalarField pTot
    (
        IOobject
        (
            "pTot",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p + 0.5*rho*(U&U)
    );

    volScalarField NOx
    (
        IOobject
        (
            "NOx",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

