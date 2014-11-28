for(label i=0; i<nPassives; i++)
{

  volScalarField& v = vs[i];

    fvScalarMatrix VEqn
    (
     //fvm::ddt(rho, v)
        fvm::div(phi, v)
      - fvm::laplacian(turbulence->muEff(), v)
    );

    VEqn.relax();
    VEqn.solve(mesh.solver("vs"));

}

//scalar cNorm = 2.0*(nPassives-1.0)/nPassives;
// if the inlet values are set to 1, cNorm = 0.5;
// otherwise cNorm = 0.5/(nPassives * inletValue/nPassives)
scalar inletValue = 1.0;
scalar cNorm = 0.5/inletValue;
errMix = dimensionedScalar("errMix", dimensionSet(0,0,0,0,0), 0.0);
dimensionedScalar average(dimensionedScalar("average", dimensionSet(0,0,0,0,0), inletValue/nPassives));
for(label i=0; i<nPassives; i++)
{
    volScalarField dev = mag(vs[i] - average);
    errMix += cNorm*dev;
}

bland = 1.0 - errMix;
