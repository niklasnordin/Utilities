    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo.rho()
    );

    volScalarField& p = thermo.p();
    volScalarField& h = thermo.h();
    const volScalarField& psi = thermo.psi();
    const volScalarField& T = thermo.T();

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

    #include "compressibleCreatePhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

    dimensionedScalar rhoMax
    (
        mesh.solutionDict().subDict("SIMPLE").lookup("rhoMax")
    );

    dimensionedScalar rhoMin
    (
        mesh.solutionDict().subDict("SIMPLE").lookup("rhoMin")
    );

    dimensionedScalar pMin
    (
        mesh.solutionDict().subDict("SIMPLE").lookup("pMin")
    );

    dimensionedScalar pMax
    (
        mesh.solutionDict().subDict("SIMPLE").lookup("pMax")
    );

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::RASModel> turbulence
    (
        compressible::RASModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    dimensionedScalar initialMass = fvc::domainIntegrate(rho);

    wordList htcBCTypes
    (
        U.boundaryField().size(),
        calculatedFvPatchScalarField::typeName
    );

    for (label i=0; i<htcBCTypes.size(); i++)
    {
        if (isA<wallFvPatch>(mesh.boundary()[i]))
        {
	  //htcBCTypes[i] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField htc
    (
        IOobject
        (
            "htc",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0.0),
	htcBCTypes
    );


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
