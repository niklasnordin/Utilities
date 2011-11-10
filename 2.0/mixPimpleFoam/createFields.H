    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<basicPsiThermo> pThermo
    (
        basicPsiThermo::New(mesh)
    );
    basicPsiThermo& thermo = pThermo();

    volScalarField& p = thermo.p();
    volScalarField& h = thermo.h();
    const volScalarField& psi = thermo.psi();
    const volScalarField& T = thermo.T();

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

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    wordList BCTypes
    (
        U.boundaryField().size(),
        calculatedFvPatchScalarField::typeName
    );


    volScalarField magGradU
    (
        IOobject
        (
            "magGradU",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,1,-2,0,0), 0.0),
        BCTypes
    );

    volVectorField Up
    (
        IOobject
        (
            "Up",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0), vector::zero),
        BCTypes
    );

    volVectorField Uav
    (
        IOobject
        (
            "Uav",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("zero", dimensionSet(0,1,-1,0,0), vector::zero),
        BCTypes
    );

    scalar startAveraging
    (
        readScalar(mesh.solutionDict().subDict("PIMPLE").lookup("startAveraging"))
    );

    #include "compressibleCreatePhi.H"

    dimensionedScalar rhoMax
    (
        mesh.solutionDict().subDict("PIMPLE").lookup("rhoMax")
    );

    dimensionedScalar rhoMin
    (
        mesh.solutionDict().subDict("PIMPLE").lookup("rhoMin")
    );

    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );

    Info<< "Creating field DpDt\n" << endl;
    volScalarField DpDt
    (
        fvc::DDt(surfaceScalarField("phiU", phi/fvc::interpolate(rho)), p)
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
