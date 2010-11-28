
for(label is=0; is<nSTL; is++)
{

    label nMeshOutline = meshOutline.size();
    edgeDetails dummy;
    dummy.d0.featurei = -1;
    dummy.d1.featurei = -1;
    List<edgeDetails> edgeBelongsToFeature(nMeshOutline, dummy);

    forAll(stlFeatures[is], j)
    {

        const edge& feat = stlFeatures[is][j];
	vector n1 = globalStlPoints[is][feat[1]] - globalStlPoints[is][feat[0]];
	n1 /= mag(n1);
	
	// first collect all the edges close enough to be mapped to this feature line
	forAll(meshOutline, i)
	{
	    
	    const edge& e = meshOutline[i];
	    label i0 = addr[e[0]];
	    label i1 = addr[e[1]];
	    vector p0 = newPoints[i0];
	    vector p1 = newPoints[i1];
	    scalar scale = mag(p1 - p0);
	    
	    vector n2 = (p1 - p0)/scale;
	    scalar ac = mag(n1 & n2);
	    
	    // only align features which aren't close to perpendicular
	    if (ac > cosExclude)
	    {
		scalar lam0 = GREAT;
		vector px0 = vector::zero;
		scalar d0 = distance(globalStlPoints[is], feat, p0, lam0, px0);
		scalar lam1 = GREAT;
		vector px1 = vector::zero;
		scalar d1 = distance(globalStlPoints[is], feat, p1, lam1, px1);
		
		details det0, det1;
		det0.pointi = i0;
		det0.pointo = i1;
		det0.featurei = j;
		det0.distance = d0;
                det0.distanceo = d1;
		det0.lambda = lam0;
		det0.p = px0;
		det0.po = px1;
		det0.lambdaOther = lam1;
		det0.cos = ac;
		
		det1.pointi = i1;
		det1.pointo = i0;
		det1.featurei = j;
		det1.distance = d1;
                det1.distanceo = d0;
		det1.lambda = lam1;
		det1.p = px1;
		det1.po = px0;
		det1.lambdaOther = lam0;
		det1.cos = ac;
		
		// chose the feature that has the smallest angle to the edge
		// instead of letting the distance alone dictate to which feature the edge belongs.
		// this allows points to belong to multiple features and makes it
		// possible to capture the corners as well

		if ( (det0.lambda >= minFit) && (det0.lambda <= maxFit) && ( d0 < tolerance*scale) )
		{
		    if (edgeBelongsToFeature[i].d0.featurei == -1)
		    {
			edgeBelongsToFeature[i].d0 = det0;
		    }
		    else
		    {
			scalar ar0 = edgeBelongsToFeature[i].d0.cos;

			// if both edges are close to parallel, use distance, otherwise use angle
			if ( (ar0 > cosParAngle) && (ac > cosParAngle) )
			{
			    if (d0 < edgeBelongsToFeature[i].d0.distance)
			    {
				edgeBelongsToFeature[i].d0 = det0;
			    }
			}
			else if ( ac > ar0 )
			{
			    edgeBelongsToFeature[i].d0 = det0;
			}
		    }
		}
		
		if ( (det1.lambda >= minFit) && (det1.lambda <= maxFit) && ( d1 < tolerance*scale) )
		{
		    if (edgeBelongsToFeature[i].d1.featurei == -1)
		    {
			edgeBelongsToFeature[i].d1 = det1;
		    }
		    else
		    {
			scalar ar0 = edgeBelongsToFeature[i].d1.cos;

			// if both edges are parallel, use distance, otherwise use angle
			if ( (ar0 > cosParAngle) && (ac > cosParAngle) )
			{
			    if (d1 < edgeBelongsToFeature[i].d1.distance)
			    {
				edgeBelongsToFeature[i].d1 = det1;
			    }
			}
			else if ( ac > ar0 )
			{
			    edgeBelongsToFeature[i].d1 = det1;
			}
			
		    }
		}
	    }
	    
	} // end forAll(meshOutline, i)

    } // end forAll(stlFeatures[is]...
    
    // now every edge is assigned to 1-2 feature line(s)
    // check for overlaps
    label nEdgesToMap = 0;

#   include "checkEdges.H"

    scalar sumMovement = 0.0;
    forAll(meshOutline, i)
    {
        bool map = false;
	
	if ( edgeBelongsToFeature[i].d0.featurei != -1 )
	{
	    map = true;
	    label i0 = edgeBelongsToFeature[i].d0.pointi;
	    vector dx = relaxation*(edgeBelongsToFeature[i].d0.p - newPoints[i0]);
	    vector np0 = newPoints[i0] + dx;

	    if (edgeBelongsToFeature[i].d1.featurei == -1)
	    {
		// contract line
		label i1 = edgeBelongsToFeature[i].d0.pointo;
		newPoints[i1] += 0.01*relaxation*(np0 - newPoints[i1]);
	    }

	    //if (!checkOverlap(newPoints, isEdgePoint, np0, i0, overlap))
	    {
		newPoints[i0] = np0;
		sumMovement += mag(dx);
	    }
	}
	
	if (edgeBelongsToFeature[i].d1.featurei != -1)
	{
	    map = true;
	    label i1 = edgeBelongsToFeature[i].d1.pointi;
	    vector dx = relaxation*(edgeBelongsToFeature[i].d1.p - newPoints[i1]);
	    vector np1 = newPoints[i1] + dx;
	    if ( edgeBelongsToFeature[i].d0.featurei == -1 )
	    {
		// contract line
		//label i0 = edgeBelongsToFeature[i].d1.pointo;
		//newPoints[i0] += 0.01*relaxation*(np1 - newPoints[i0]);
	    }

	    //if (!checkOverlap(newPoints, isEdgePoint, np1, i1, overlap))
	    {
		newPoints[i1] = np1;
		sumMovement += mag(dx);
	    }
	}
	if (map) nEdgesToMap++;
    }
    
    // move the corner points for the last iteration
    /*
    if (ni == (nIterations-1))
    {
	List<label> closestID(cornerPoints[is].size(), 0);
	List<scalar> closestDist(cornerPoints[is].size(), GREAT);
	forAll(newPoints, i)
	{
	    forAll(closestID, j)
	    {
	        scalar dist = mag(cornerPoints[is][j] - newPoints[i]);
		if (dist < closestDist[j])
		{
		    closestDist[j] = dist;
		    closestID[j] = i;
		}
	    }
	}
	forAll(closestID, i)
	{
	    if (closestDist[i] < smallestEdgeLength)
	    {
	      //newPoints[closestID[i]] = cornerPoints[is][i];
	    }
	}
    }
    */
    Info << " ( " << nEdgesToMap << " / " << meshOutline.size() << " ). sum(move) = " << sumMovement << ". ";
}  // end for(label is=0...

