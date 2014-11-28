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

if (time > startAveraging )
{
    if ( (time - dt) < startAveraging )
    {
        UpAv = Up;
	UprocAv = Uproc;
    }
    else
    {
        UpAv = ( UpAv*(time-startAveraging) + Up*dt ) / (time-startAveraging + dt);
        UprocAv = ( UprocAv*(time-startAveraging) + Uproc*dt ) / (time-startAveraging + dt);
    }
}

/*
const polyBoundaryMesh& patches = mesh.boundaryMesh();
forAll(patches, i)
{
    magGradU.boundaryField()[i] = mag(U.boundaryField()[i].snGrad());
}
*/

}
