// RAS
{
  //    const compressible::RASModel& rasModel = U.db().lookupObject<compressible::RASModel>("RASProperties");
    scalar Pr = 1.0;
    scalar sigma_t = 1.0;
    scalar Cmu = 0.09;
    scalar Cmu4 = ::pow(Cmu, 0.25);
    scalar K = 0.41;
    scalar E = 9.0;
    scalar Y = 9.24*(::pow(Pr/sigma_t, 0.75) - 1.0)*(1.0 + 0.28*::exp(-0.007*Pr/sigma_t));

    surfaceScalarField snGradT = fvc::snGrad(T);

    for (label i=0; i<htcBCTypes.size(); i++)
    {
        if (isA<wallFvPatch>(mesh.boundary()[i]))
        {
	    scalarField Tp = T.boundaryField()[i].patchInternalField();
            const scalarField& rhop = rho.boundaryField()[i];
	    scalarField kp = turbulence->k()().boundaryField()[i].patchInternalField();

            scalarField yp = turbulence->yPlus(i, Cmu);
	    scalarField Cpp = thermo.Cp(Tp, i);

            scalarField Tstar = Cmu4*(log(E*yp)/K + Y);

	    scalarField qw = Cmu4*pow(kp, 0.5)*rhop*Cpp/Tstar;

	    htc.boundaryField()[i] == qw;

        }
    }
}

