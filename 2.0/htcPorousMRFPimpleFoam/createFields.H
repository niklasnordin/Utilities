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
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
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


    wordList htcBCTypes
    (
     U.boundaryField().size(),
     calculatedFvPatchScalarField::typeName
     );
/*
    for (label i=0; i<htcBCTypes.size(); i++)
    {
        if (isA<wallFvPatch>(mesh.boundary()[i]))
        {
	//htcBCTypes[i] = fixedValueFvPatchScalarField::typeName;                                                                                         
        }
    }
*/
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

