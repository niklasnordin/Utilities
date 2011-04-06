    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho.dimensionedInternalField()*g
      + sprayParcels.SU()
    );

    if (momentumPredictor)
    {
        solve(UEqn == -fvc::grad(p));
    }

