
    fvScalarMatrix v1Eqn
    (
        fvm::ddt(rho, v1)
      + fvm::div(phi, v1)
      - fvm::laplacian(turbulence->muEff(), v1)
    );

    v1Eqn.relax();
    v1Eqn.solve();

    fvScalarMatrix v2Eqn
    (
        fvm::ddt(rho, v2)
      + fvm::div(phi, v2)
      - fvm::laplacian(turbulence->muEff(), v2)
    );

    v2Eqn.relax();
    v2Eqn.solve();

    fvScalarMatrix v3Eqn
    (
        fvm::ddt(rho, v3)
      + fvm::div(phi, v3)
      - fvm::laplacian(turbulence->muEff(), v3)
    );

    v3Eqn.relax();
    v3Eqn.solve();

    fvScalarMatrix v4Eqn
    (
        fvm::ddt(rho, v4)
      + fvm::div(phi, v4)
      - fvm::laplacian(turbulence->muEff(), v4)
    );

    v4Eqn.relax();
    v4Eqn.solve();

    fvScalarMatrix v5Eqn
    (
        fvm::ddt(rho, v5)
      + fvm::div(phi, v5)
      - fvm::laplacian(turbulence->muEff(), v5)
    );

    v5Eqn.relax();
    v5Eqn.solve();

    fvScalarMatrix v6Eqn
    (
        fvm::ddt(rho, v6)
      + fvm::div(phi, v6)
      - fvm::laplacian(turbulence->muEff(), v6)
    );

    v6Eqn.relax();
    v6Eqn.solve();

    fvScalarMatrix v7Eqn
    (
        fvm::ddt(rho, v7)
      + fvm::div(phi, v7)
      - fvm::laplacian(turbulence->muEff(), v7)
    );

    v7Eqn.relax();
    v7Eqn.solve();

    fvScalarMatrix v8Eqn
    (
        fvm::ddt(rho, v8)
      + fvm::div(phi, v8)
      - fvm::laplacian(turbulence->muEff(), v8)
    );

    v8Eqn.relax();
    v8Eqn.solve();

    scalar fr = 4.0;
    volScalarField c1=fr*v1*(1.0-v1);
    volScalarField c2=fr*v2*(1.0-v2);
    volScalarField c3=fr*v3*(1.0-v3);
    volScalarField c4=fr*v4*(1.0-v4);
    volScalarField c5=fr*v5*(1.0-v5);
    volScalarField c6=fr*v6*(1.0-v6);
    volScalarField c7=fr*v7*(1.0-v7);
    volScalarField c8=fr*v8*(1.0-v8);
    mix = c1*c2*c3*c4*c5*c6*c7*c8;
    addMix = c1 + c2 + c3 + c4 + c5 + c6 + c7 + c8;

    Info << "max(mix) = " << max(mix).value() 
    << ", max(addMix) = " << max(addMix).value()
    << endl;
