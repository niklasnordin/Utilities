// Solve the Momentum equation

tmp<fvVectorMatrix> UEqn
(
    pZones.ddt(rho, U)
  + fvm::div(phi, U)
  + turbulence->divDevRhoReff(U)
  ==
    rho*g
  + dieselSpray.momentumSource()

);

if (oCorr == nOuterCorr-1)
{
    UEqn().relax(1);
}
else
{
    UEqn().relax();
}

pZones.addResistance(UEqn());

volScalarField rUA = 1.0/UEqn().A();

if (momentumPredictor)
{
    if (oCorr == nOuterCorr-1)
    {
        solve(UEqn() == -fvc::grad(p), mesh.solver("UFinal"));
    }
    else
    {
        solve(UEqn() == -fvc::grad(p));
    }
}
else
{
    U = rUA*(UEqn().H() - fvc::grad(p));
    U.correctBoundaryConditions();
}
