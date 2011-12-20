{

    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + fvm::div(phi, hs)
      - fvm::laplacian(turb.alphaEff(), hs)
     ==
        DpDt
      + parcels[i].Sh(hs)
	//+ rad.Sh(thermo)
    );

    hsEqn.relax();
    hsEqn.solve(mesh.solver(hs.select(finalIter)));

    thermo.correct();
    //rad.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
