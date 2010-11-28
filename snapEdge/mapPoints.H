

    Info << "Mapping to stl...";
    flush(Info);
    // map points to stl
    forAll(patches, i)
    {
        const labelList& addr = patches[i].meshPoints();
	const vectorField& bPoints = patches[i].localPoints();

	// dont use the localPoint, use the newPoints
	forAll(bPoints, j)
	{
	    vector& p = newPoints[addr[j]];
	    
	    // find the stl triangle it should be in
	    forAll(stlSurface, k)
	    {
	        vector Sf = stlSf[k];
		Sf /= mag(Sf);
		const face& f = stlSurface[k];
		const vector& p0 = stlGlobalPoints[f[0]];
		const vector& p1 = stlGlobalPoints[f[1]];
		const vector& p2 = stlGlobalPoints[f[2]];
		scalar Area = mag((p1-p0)^(p2-p0));
		scalar dist = (p-p0)&Sf;

		if (mag(dist) < mapTolerance*::sqrt(Area))
		{
		    vector px = p - dist*Sf;
		    
		    if (insideTriangle(p0, p1, p2, px))
		    {
			p = px;
		    }
		}
	    }
	}
    }

    Info << "Done!" << endl;

