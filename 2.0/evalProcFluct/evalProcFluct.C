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
    evalProcFluct

Description

\*---------------------------------------------------------------------------*/
#include "calc.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    bool writeResults = !args.optionFound("noWrite");

    IOobject Uavheader
    (
        "Uav",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    IOobject Upheader
    (
        "Up",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );


    // Check Uav and Up exists
    if (Uavheader.headerOk() && Upheader.headerOk())
    {
        volVectorField Uav(Uavheader, mesh);
        volVectorField Up(Upheader, mesh);
	dimensionedScalar Usmall("Usmall", Uav.dimensions(), SMALL);

        autoPtr<volScalarField> UprocPtr;

	UprocPtr.set
	(
	    new volScalarField
	    (
	        IOobject
		(
		    "Uproc",
		    runTime.timeName(),
		    mesh
		),
		mag(Up)/(mag(Uav) + Usmall)
	    )
	);

        if (writeResults)
	{
	    UprocPtr().write();
	}

    }

    Info<< "\nEnd\n" << endl;

}


// ************************************************************************* //
