label nPassives(readLabel(simple.dict().lookup("nPassives")));

Info << "number of passive scalars are " << nPassives << endl;

PtrList<volScalarField> vs(nPassives);

for(label i=0; i<nPassives; i++)
{
    std::ostringstream oss;
    oss << "v" << i+1;
    word name( oss.str() );

    vs.set
    (
        i,
	new volScalarField
	(
	    IOobject
	    (
	        name,
		mesh.time().timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	    ),
	    mesh
        )
    );

}

volScalarField errMix
(
    IOobject
    (
        "errMix",
       	mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
       	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("errMix", dimensionSet(0,0,0,0,0), 0.0)
);

volScalarField bland
(
    IOobject
    (
        "bland",
       	mesh.time().timeName(),
        mesh,
        IOobject::NO_READ,
       	IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("bland", dimensionSet(0,0,0,0,0), 0.0)
);
