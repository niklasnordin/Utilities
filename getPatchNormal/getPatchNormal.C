/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    rhoPimpleFoam

Description
    Transient solver for laminar or turbulent flow of compressible fluids
    for HVAC and similar applications.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    argList::addOption
    (
        "patchName",
        "word",
        "specify the name of the patch to calculate the average normal for"
    );

    #include "setRootCase.H"

    bool patchFound = args.optionFound("patchName");
    if (!patchFound)
    {
        Info << "please specify the name of the patch" << endl;
        return 0;
    }

    #include "createTime.H"
    #include "createMesh.H"

    word patchName = args.optionLookupOrDefault<word>("patchName", "inlet");
    label patchI = mesh.boundaryMesh().findPatchID(patchName);

    if (patchI == -1)
    {
        Info << "could not find patch " << patchName << endl;
	return 0;
    }
    //const vectorField& Sf = mesh.boundaryMesh()[patchI].faceNormals();
    const vectorField& Sf = mesh.boundaryMesh()[patchI].faceAreas();

    vector vSum = vector::zero;
    scalar aSum = 0.0;
    forAll(Sf, i)
    {
        vSum += Sf[i];
        aSum += mag(Sf[i]);
    }
    Info << "NB! vector is pointing outward" << endl;
    Info << "vector = " << vSum/aSum << endl;
    Info << "area = " << aSum << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
