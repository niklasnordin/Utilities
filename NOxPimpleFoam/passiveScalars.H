    label nInlets = 8;
    label nInlet = nInlets*rndGen.sample01<scalar>();

    for(label i=0; i<nInlets; i++)
    {

        OStringStream num;
	label j = i+1;
	num << j;
        word patchName = "STL_inlet" + num.str();
        label patchID = mesh.boundaryMesh().findPatchID(patchName);
	if (patchID == -1)
	{
	    Info << "warning!!! Could not find patch " << patchName << endl;
	}

	if (i == nInlet)
	{
	    // set patch value to max
	    NOx.boundaryField()[patchID] == 160.0;
	}
	else
	{
	    // set patch value to min
	    NOx.boundaryField()[patchID] == 0.0;
	}
    }

    fvScalarMatrix NOxEqn
    (
        fvm::ddt(rho, NOx)
      + fvm::div(phi, NOx)
      - fvm::laplacian(turbulence->muEff(), NOx)
    );

    NOxEqn.relax();
    NOxEqn.solve();
