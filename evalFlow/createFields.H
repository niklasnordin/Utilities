Info<< nl << "Reading thermophysicalProperties" << endl;

autoPtr<psiChemistryModel> pChemistry
(
    psiChemistryModel::New(mesh)
);
psiChemistryModel& chemistry = pChemistry();

hsCombustionThermo& thermo = chemistry.thermo();

basicMultiComponentMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

word inertSpecie(thermo.lookup("inertSpecie"));

if (!composition.contains(inertSpecie))
{
    FatalErrorIn(args.executable())
        << "Specified inert specie '" << inertSpecie << "' not found in "
        << "species list. Available species:" << composition.species()
        << exit(FatalError);
}

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh
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


volScalarField& p = thermo.p();
//const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();
//volScalarField& hs = thermo.hs();

#include "compressibleCreatePhi.H"

Info << "Creating turbulence model.\n" << nl;
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


Info << "Constructing Spray" << endl;

PtrList<gasThermoPhysics> gasProperties(Y.size());
forAll(gasProperties, i)
{
    gasProperties.set
    (
        i,
        new gasThermoPhysics
        (
            dynamic_cast<const reactingMixture<gasThermoPhysics>&>
                (thermo).speciesData()[i]
        )
    );
}
