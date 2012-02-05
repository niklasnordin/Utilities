
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


    volScalarField v1
    (
        IOobject
        (
            "v1",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField v2
    (
        IOobject
        (
            "v2",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField v3
    (
        IOobject
        (
            "v3",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField v4
    (
        IOobject
        (
            "v4",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    volScalarField v5
    (
        IOobject
        (
            "v5",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField v6
    (
        IOobject
        (
            "v6",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField v7
    (
        IOobject
        (
            "v7",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField v8
    (
        IOobject
        (
            "v8",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    wordList htcBCTypes
    (
        U.boundaryField().size(),
        zeroGradientFvPatchScalarField::typeName
    );


    volScalarField mix
    (
        IOobject
        (
            "mix",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0.0),
	htcBCTypes
    );

    volScalarField addMix
    (
        IOobject
        (
            "addMix",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0.0),
	htcBCTypes
    );
