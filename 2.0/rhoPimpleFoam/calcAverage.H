{
scalar time = runTime.value();
scalar dt = runTime.deltaT().value();

if (time > startAveraging )
{
    if ( (time - dt) < startAveraging )
    {
        Uav = U;
    }
    else
    {
        Uav = ( Uav*(time-startAveraging) + U*dt ) / (time-startAveraging + dt);
    }
}

dimensionedScalar Usmall("Usmall", Uav.dimensions(), SMALL);

Up = U - Uav;
Uproc = mag(Up)/(mag(Uav)+Usmall);

const polyBoundaryMesh& patches = mesh.boundaryMesh();
forAll(patches, i)
{
    magGradU.boundaryField()[i] = mag(U.boundaryField()[i].snGrad());
}
 

}
