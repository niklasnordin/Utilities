Info << "Constructing Spray" << endl;

/*
scalar gasMass0 = fvc::domainIntegrate(rho).value();

if (dieselSpray.twoD())
{
    gasMass0 *= 2.0*mathematicalConstant::pi/dieselSpray.angleOfWedge();
}

gasMass0 -=
    dieselSpray.injectedMass(runTime.value()) - dieselSpray.liquidMass();
*/

thermoSprayCloud sprayParcels
(
    "sprayCloud",
    rho,
    U,
    g,
    thermo
);
