    typedef regionModels::surfaceFilmModels::surfaceFilmModel filmModelType;

    // Initialise fluid field pointer lists
    PtrList<hsCombustionThermo> thermoFluid(fluidRegions.size());
    PtrList<volScalarField> rhoFluid(fluidRegions.size());
    PtrList<volScalarField> KFluid(fluidRegions.size());
    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
    PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
    PtrList<compressible::turbulenceModel> turbulence(fluidRegions.size());
    PtrList<volScalarField> p_rghFluid(fluidRegions.size());
    PtrList<volScalarField> ghFluid(fluidRegions.size());
    PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
    //PtrList<radiation::radiationModel> radiation(fluidRegions.size());
    PtrList<volScalarField> DpDtFluid(fluidRegions.size());

    List<scalar> initialMassFluid(fluidRegions.size());

    PtrList<psiChemistryModel> pChemistry(fluidRegions.size());
    PtrList<basicReactingCloud> parcels(fluidRegions.size());
    PtrList<SLGThermo> slgThermos(fluidRegions.size());
    PtrList<autoPtr<filmModelType> > tsurfaceFilm(fluidRegions.size());
    PtrList<filmModelType> surfaceFilm(fluidRegions.size());

//autoPtr<filmModelType> tsurfaceFilm(filmModelType::New(mesh, g));
//filmModelType& surfaceFilm = tsurfaceFilm();


    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermoFluid\n" << endl;
	pChemistry.set
	(
	    i,
	    psiChemistryModel::New(fluidRegions[i])
        );

	psiChemistryModel& chemistry = pChemistry[i];
	hsCombustionThermo& thermo = chemistry.thermo();

	
        thermoFluid.set
        (
            i,
            &thermo
        );

	slgThermos.set
	(
	    i,
	    new SLGThermo(fluidRegions[i], thermo)
        );
	
        Info<< "    Adding to rhoFluid\n" << endl;
        rhoFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                thermo.rho()
            )
        );

        Info<< "    Adding to KFluid\n" << endl;
        KFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "K",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                thermo.Cp()*thermo.alpha()
            )
        );

        Info<< "    Adding to UFluid\n" << endl;
        UFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to phiFluid\n" << endl;
        phiFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(rhoFluid[i]*UFluid[i])
                    & fluidRegions[i].Sf()
            )
        );

        Info<< "    Adding to gFluid\n" << endl;
        gFluid.set
        (
            i,
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );

        Info<< "    Adding to turbulence\n" << endl;
        turbulence.set
        (
            i,
            autoPtr<compressible::turbulenceModel>
            (
                compressible::turbulenceModel::New
                (
                    rhoFluid[i],
                    UFluid[i],
                    phiFluid[i],
                    thermo
                )
            ).ptr()
        );

        Info<< "    Adding to ghFluid\n" << endl;
        ghFluid.set
        (
            i,
            new volScalarField("gh", gFluid[i] & fluidRegions[i].C())
        );

        Info<< "    Adding to ghfFluid\n" << endl;
        ghfFluid.set
        (
            i,
            new surfaceScalarField("ghf", gFluid[i] & fluidRegions[i].Cf())
        );

        p_rghFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "p_rgh",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        // Force p_rgh to be consistent with p
        p_rghFluid[i] = thermo.p() - rhoFluid[i]*ghFluid[i];
	/*
        radiation.set
        (
            i,
            radiation::radiationModel::New(thermo.T())
        );
	*/
        initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

        Info<< "    Adding to DpDtFluid\n" << endl;
        DpDtFluid.set
        (
            i,
            new volScalarField
            (
                "DpDt",
                fvc::DDt
                (
                    surfaceScalarField
                    (
                        "phiU",
                        phiFluid[i]/fvc::interpolate(rhoFluid[i])
                    ),
                    thermo.p()
                )
            )
        );

	parcels.set
	(
	    i,
	    new basicReactingCloud
	    (
	        "reactingCloud",
		rhoFluid[i],
		UFluid[i],
		gFluid[i],
		slgThermos[i]
	    )
	);

        //autoPtr<filmModelType> tsurfaceFilmi(filmModelType::New(fluidRegions[i], gFluid[i]));
	//tsurfaceFilm[i] = filmModelType::New(fluidRegions[i], gFluid[i]);
	//surfaceFilm[i] = (tsurfaceFilm[i])();

	
	surfaceFilm.set
	(
	    i,
	    filmModelType::New(fluidRegions[i], gFluid[i])
        );
	
	//surfaceFilm.set(i, (tsurfaceFilm[i])());
	//filmModelType& surfaceFilm = tsurfaceFilm();


    }
