
basicMultiComponentMixture& composition = thermoFluid[i].composition();
PtrList<volScalarField>& Y = composition.Y();

const word inertSpecie(thermo.lookup("inertSpecie"));

if (!composition.contains(inertSpecie))
{
    FatalErrorIn(args.executable())
      << "Specified inert specie '" << inertSpecie << "' not found in "
      << "species list. Available species:" << composition.species()
      << exit(FatalError);
}

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, ii)
{
  fields.add(Y[ii]);
}
fields.add(hs);

tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);



{
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, j)
    {
        if (Y[j].name() != inertSpecie)
        {
            volScalarField& Yi = Y[j];
            solve
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turb.muEff(), Yi)
              ==
                parcels[i].SYi(j, Yi),
		//+ chemistry.RR(i)().dimensionedInternalField(),
                mesh.solver("Yi")
            );

            Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = j;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
