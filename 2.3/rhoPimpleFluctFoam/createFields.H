    Info<< "Reading thermophysical properties\n" << endl;

    autoPtr<psiThermo> pThermo
    (
        psiThermo::New(mesh)
    );
    psiThermo& thermo = pThermo();
    thermo.validate(args.executable(), "h", "e");

    volScalarField& p = thermo.p();
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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "compressibleCreatePhi.H"

    dimensionedScalar rhoMax(pimple.dict().lookup("rhoMax"));
    dimensionedScalar rhoMin(pimple.dict().lookup("rhoMin"));

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

    Info<< "Creating field dpdt\n" << endl;
    volScalarField dpdt
    (
        IOobject
        (
            "dpdt",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
    );

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));


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

    volScalarField Uproc
    (
        IOobject
        (
	 "Uproc",
	 runTime.timeName(),
	 mesh,
	 IOobject::READ_IF_PRESENT,
	 IOobject::AUTO_WRITE
	 ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0.0),
        BCTypes
     );

    volScalarField UprocAv
    (
        IOobject
        (
	 "UprocAv",
	 runTime.timeName(),
	 mesh,
	 IOobject::READ_IF_PRESENT,
	 IOobject::AUTO_WRITE
	 ),
        mesh,
        dimensionedScalar("zero", dimensionSet(0,0,0,0,0), 0.0),
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

    volVectorField UpAv
    (
        IOobject
        (
	 "UpAv",
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
