{
    solve
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
     ==
       DpDt
     + sprayParcels.Sh()
     + chemistrySh
    );

    thermo.correct();
}
