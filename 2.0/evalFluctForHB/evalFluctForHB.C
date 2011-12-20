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
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:


void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{

    OFstream logSummaryFile
    (
       runTime.path()/("Udev_" + runTime.timeName() + ".dat")
    );

    //bool writeResults = !args.optionFound("noWrite");

    vector normal(-1, 0, 0);
    IOobject Uheader
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check Uav and Up exists
    if (Uheader.headerOk())
    {
        volVectorField U(Uheader, mesh);
	dimensionedScalar Usmall("Usmall", U.dimensions(), SMALL);
	word key("vulkanOut");

	const labelList zoneID = mesh.cellZones().findIndices(key);

	if (zoneID.size() == 0)
	{
	    FatalError << "no zone found" << abort(FatalError);
	}
        const labelList& cells = mesh.cellZones()[zoneID[0]];

	vector Uav = vector::zero;
	scalar vol = 0.0;
	scalar Udev = 0.0;
	scalar torque = 0.0;
	vector center = vector::zero;

        forAll(cells, i)
        {
	    label celli = cells[i];
	    scalar voli = mesh.V()[celli];
	    vector C = mesh.C()[celli];
	    vol += voli;
	    center += C*voli;
        }
	reduce(center, sumOp<vector>());
	reduce(vol, sumOp<scalar>());
	center /= vol;

        forAll(cells, i)
        {
	    label celli = cells[i];
	    scalar voli = mesh.V()[celli];
	    vector C = mesh.C()[celli];
	    Uav += U[celli]*voli;
	    vector Utau = U[celli] - ((U[celli] & normal)*normal);
	    scalar Udevi = mag(Utau);
	    Udev += Udevi*voli;
	    vector radie = C - center;
	    vector torqi = Utau ^ radie;
	    torque += mag(torqi)*voli;
        }
	reduce(Uav, sumOp<vector>());
	reduce(Udev, sumOp<scalar>());
	reduce(torque, sumOp<scalar>());

	Uav /= vol;
	Udev /= vol;
	torque /= vol;

	Info << "Udev = " << Udev << ", Uav = " << Uav << ", torque = " << torque << endl;
	logSummaryFile << Udev << " " << torque << endl;
    }

    Info<< "\nEnd\n" << endl;

}


// ************************************************************************* //
