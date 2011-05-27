{
    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(turbulence->alphaEff(), hs)
     ==
        DpDt
     + dieselSpray.heatTransferSource()().dimensionedInternalField()
     + chemistrySh

    );

    if (oCorr == nOuterCorr-1)
    {
        hsEqn.relax();
        hsEqn.solve(mesh.solver("hsFinal"));
    }
    else
    {
        hsEqn.relax();
        hsEqn.solve();
    }

    thermo.correct();
}
