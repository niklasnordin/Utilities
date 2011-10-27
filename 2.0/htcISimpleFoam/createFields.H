    Info<< "Reading field p\n" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
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

    #include "createPhi.H"


    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell(p, mesh.solutionDict().subDict("SIMPLE"), pRefCell, pRefValue);

    singlePhaseTransportModel laminarTransport(U, phi);

    autoPtr<incompressible::RASModel> turbulence
    (
        incompressible::RASModel::New(U, phi, laminarTransport)
    );

    wordList htcBCTypes
    (
        U.boundaryField().size(),
	calculatedFvPatchScalarField::typeName
    );

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
        p + 0.5*(U&U)
    );


    const dictionary& transpDict = U.db().lookupObject<dictionary>("transportProperties");

    const dimensionedScalar rhoDim(transpDict.lookup("rho"));
    const scalar rho = rhoDim.value();
    const scalar Cp(readScalar(transpDict.lookup("Cp")));
