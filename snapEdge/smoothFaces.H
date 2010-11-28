if (nSmootherIterations>0)
{
    Info << "Smoothing patch faces...";
    flush(Info);
}

for(label ix=0; ix<nSmootherIterations; ix++)
{
    forAll(snapPatches, index)
    {
        label id = mesh.boundaryMesh().findPatchID(snapPatches[index]);
	const faceList& bFaces = patches[id];
	//const labelList& addr = patches[id].meshPoints();
	const vectorField& bCf = patches[id].faceCentres();

	forAll(bFaces, j)
	{
  	    const face& f = bFaces[j];
	    label n = f.size();

	    forAll(f, i)
	    {
	        label i0 = f[i];
		label i1 = f[(i+1) % n];
		label i2 = f[(i+2) % n];

		const vector& p0 = newPoints[i0];
		vector& p1 = newPoints[i1];
		const vector& p2 = newPoints[i2];

		vector v0 = p0 - p1;
		v0 /= mag(v0);
		vector v1 = p2 - p1;
		v1 /= mag(v1);
		scalar dot = mag(v1 & v0);
		if (dot > 0.999999)
		{
		    const vector& Cf = bCf[j];
		    p1 = 0.5*(p0 + p2);
		    vector n = p1 - Cf;
		    p1 += straightEdgePerturbation*n;
		    //Info << "Straight line, dot = " << dot << ", i1 = " << i1 << endl;
		}

	    }
	}
    }
}
