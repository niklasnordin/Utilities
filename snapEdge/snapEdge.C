/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    snapEdge

Description
    snapEdges to feature lines

Usage
    Options are:

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "ReadFields.H"
#include "pointFields.H"
#include "transformField.H"
#include "transformGeometricField.H"
#include "IStringStream.H"
#include "Random.H"
#include "triSurface.H"

using namespace Foam;
using namespace Foam::mathematicalConstant;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "types.H"
#include "functions.H"

//  Main program:
int main(int argc, char *argv[])
{
    argList::validOptions.insert("overwrite", "");

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    const word oldInstance = mesh.pointsInstance();
    const bool overwrite = args.optionFound("overwrite");

#   include "readSnapEdgeDict.H"

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const faceZoneMesh& zones = mesh.faceZones();
    const pointField& points = mesh.points();
    vectorField newPoints(points);
    List<bool> isEdgePoint(newPoints.size(), false);
    scalar deg2rad = 3.1415926535/180.0;

    scalar cosFeature = ::cos(featureAngle*deg2rad);
    scalar cosFeature2 = ::cos(featureAngle*0.5*deg2rad);
    scalar cosExclude = ::cos(excludeEdgeAngle*deg2rad);
    scalar cosParAngle = ::cos(parallelAngle*deg2rad);

    labelListList edgePoints(0);

    label nSTL = stlFileNames.size();
    List<edgeList> stlFeatures(nSTL);
    List<vectorField> globalStlPoints(nSTL);
    List< List<vector> > cornerPoints(nSTL);
    List< List<label> > nEdgePoints(nSTL);

    forAll(stlFileNames, is)
    {
        word stlFileName(stlFileNames[is]);
        triSurface stlSurface(stlFileName);
        
        const edgeList& stlEdges = stlSurface.edges();    
        const labelListList& stlEdgeFaces = stlSurface.edgeFaces();
        const vectorField& stlSf = stlSurface.faceNormals();
        const vectorField& stlPoints = stlSurface.localPoints();
	addToList(globalStlPoints[is], stlPoints);
	nEdgePoints[is].setSize(stlPoints.size());

        Info << "Finding features for stl : " << stlFileName << endl;
        forAll(stlEdgeFaces, i)
        {
            if(!stlSurface.isInternalEdge(i))
            {
                addToList(stlFeatures[is], stlEdges[i]);
		nEdgePoints[is][stlEdges[i][0]]++;
		nEdgePoints[is][stlEdges[i][1]]++;
            }
            else
            {
	        label fi0 = stlEdgeFaces[i][0];
		label fi1 = stlEdgeFaces[i][1];
		
		scalar magF0 = mag(stlSf[fi0]);
		scalar magF1 = mag(stlSf[fi1]);
		if ( (magF0 < 1.0e-5) || (magF1 < 1.0e-5))
		{
		    Info << "Warning:: face normals are zero!?!" << endl;
		}
		vector f0 = stlSf[fi0]/magF0;
		vector f1 = stlSf[fi1]/magF1;
		
		scalar cosa = mag(f0 & f1);
		if (cosa < cosFeature)
		{
		    addToList(stlFeatures[is], stlEdges[i]);
		    nEdgePoints[is][stlEdges[i][0]]++;
		    nEdgePoints[is][stlEdges[i][1]]++;
		}
            }
        }

        if (stlFeatures[is].size() == 0)
	{
	    Info << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
	    Info << "WARNING!!! Your stl does not contain any feature lines." << endl;
	    Info << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
	}

    }

    // collect all corner points
    /*
    forAll(stlFileNames, is)
    {
        forAll(globalStlPoints[is], j)
	{
	    if (nEdgePoints[is][j] > 2)
	    {
	        addToList(cornerPoints[is], globalStlPoints[is][j]);
	    }
	}
	Info << "nCornerPoints = " << cornerPoints[is].size() << endl;
    }
    */

    Info << "Done!" << endl;
    isEdgePoint = false;

    scalar smallestEdgeLength = GREAT;
    // check that all patches are available and get the smallest length scale of the mesh
    forAll(snapPatches, index)
    {

        label id = mesh.boundaryMesh().findPatchID(snapPatches[index]);
	if (id == -1)
	{
            Info << "Could not find patch: " << snapPatches[index] <<endl;
	    Info << "Available patches are:" << endl;
	    forAll(mesh.boundaryMesh(), i)
	    {
	      Info << mesh.boundaryMesh()[i].name() << endl;
	    }
	    FatalError << abort(FatalError);
	}
	else
	{
	    const labelList& addr = patches[id].meshPoints();
	    const edgeList& edges = patches[id].edges();
	    
	    forAll(edges, i)
	    {
	        const edge& e = edges[i];
		label i0 = addr[e[0]];
		label i1 = addr[e[1]];
		isEdgePoint[i0] = true;
		isEdgePoint[i1] = true;
		scalar scale = mag(points[i0] - points[i1]);
		smallestEdgeLength = min(smallestEdgeLength, scale);
	    }
	}
    }

    // simple copy/paste/adapt for faceZones
    forAll(snapZones, index)
    {
        label id = mesh.faceZones().findZoneID(snapZones[index]);
	if (id == -1)
	{
            Info << "Could not find zone: " << snapZones[index] <<endl;
	    Info << "Available zones are:" << endl;
            Info << zones.names() << endl;
	    FatalError << abort(FatalError);
	}
	else
	{
	    const labelList& addr = zones[id]().meshPoints();
	    const edgeList& edges = zones[id]().edges();
	    
	    forAll(edges, i)
	    {
	        const edge& e = edges[i];
		label i0 = addr[e[0]];
		label i1 = addr[e[1]];
		isEdgePoint[i0] = true;
		isEdgePoint[i1] = true;
		scalar scale = mag(points[i0] - points[i1]);
		smallestEdgeLength = min(smallestEdgeLength, scale);
	    }
	}
    }
    // allow for overlap, minFit=0 and maxFit=1 implies no overlap
    scalar minFit = -mag(fitFactor);
    scalar maxFit = 1.0 + mag(fitFactor);
    scalar overlap = mag(overlapFactor*smallestEdgeLength);

    Info << "smallestEdgeLength = " << smallestEdgeLength << endl;
    for(label ni=0; ni<nIterations; ni++)
    {

        Info << "(" << ni+1 << "/" << nIterations << ") Matching edges...";
        flush(Info);

        forAll(snapPatches, index)
	{

	    label id = mesh.boundaryMesh().findPatchID(snapPatches[index]);

	    const labelList& addr = patches[id].meshPoints();
	    const vectorField& bSf = patches[id].faceNormals();
	    const edgeList& edges = patches[id].edges();
	    const labelListList& ef = patches[id].edgeFaces();

	    edgeList meshOutline(0);
	    
	    // construct the list of all the mesh edges to be used for matching
	    forAll(ef, i)
	    {
	        if (!patches[id].isInternalEdge(i))
		{
		    addToList(meshOutline, edges[i]);
		}
		else
		{
		    if (includeInterior)
		    {
			label fi0 = ef[i][0];
			label fi1 = ef[i][1];
			
			scalar mag0 = mag(bSf[fi0]);
			scalar mag1 = mag(bSf[fi1]);
			if ( (mag0<SMALL) || (mag1<SMALL) )
			{
			    Info << "Warning:: zero area faces in the mesh!!!" << endl;
			}
			vector f0 = bSf[fi0]/(mag0 + VSMALL);
			vector f1 = bSf[fi1]/(mag1 + VSMALL);
			
			// check angle between faces (use a smaller feature angle than the one used for the stl)
			scalar cosa = mag(f0 & f1);
			if (cosa < cosFeature2)
			{
			    addToList(meshOutline, edges[i]);
			}
		    }
		}
	    }
	    
#           include "newMethod.H"

	}

        // again: simple copy/paste/adapt for faceZones
        forAll(snapZones, index)
	{

	    label id = mesh.faceZones().findZoneID(snapZones[index]);

	    const labelList& addr = zones[id]().meshPoints();
	    const vectorField& bSf = zones[id]().faceNormals();
	    const edgeList& edges = zones[id]().edges();
	    const labelListList& ef = zones[id]().edgeFaces();

	    edgeList meshOutline(0);
	    
	    // construct the list of all the mesh edges to be used for matching
	    forAll(ef, i)
	    {
	        if (!zones[id]().isInternalEdge(i))
		{
		    addToList(meshOutline, edges[i]);
		}
		else
		{
		    if (includeInterior)
		    {
			label fi0 = ef[i][0];
			label fi1 = ef[i][1];
			
			scalar mag0 = mag(bSf[fi0]);
			scalar mag1 = mag(bSf[fi1]);
			if ( (mag0<SMALL) || (mag1<SMALL) )
			{
			    Info << "Warning:: zero area faces in the mesh!!!" << endl;
			}
			vector f0 = bSf[fi0]/(mag0 + VSMALL);
			vector f1 = bSf[fi1]/(mag1 + VSMALL);
			
			// check angle between faces (use a smaller feature angle than the one used for the stl)
			scalar cosa = mag(f0 & f1);
			if (cosa < cosFeature2)
			{
			    addToList(meshOutline, edges[i]);
			}
		    }
		}
	    }
	    
#           include "newMethod.H"

	}

	Info << "Moving points...";
	flush(Info);
	
	mesh.movePoints(newPoints);
	//mesh.checkMesh();
	Info << "Done!" << endl;
    }
    
#   include "smoothFaces.H"
    //#include "mapPoints.H"

    if (nSmootherIterations>0)
    {
	mesh.movePoints(newPoints);
    }
    if(!overwrite) {
        runTime++;
    } else {
        mesh.setInstance(oldInstance);
    }

 
    Info << "Writing points to Time = " << runTime.value() << endl;
    mesh.write();

    return 0;
}


// ************************************************************************* //
