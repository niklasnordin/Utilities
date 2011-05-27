/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    evalFlow

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "OFstream.H"
#include "hCombustionThermo.H"
#include "turbulenceModel.H"
#include "spray.H"
#include "psiChemistryModel.H"
#include "chemistrySolver.H"
#include "multivariateScheme.H"
#include "Switch.H"
#include "OFstream.H"
#include "volPointInterpolation.H"
#include "thermoPhysicsTypes.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"

    instantList timeDirs = timeSelector::select0(runTime, args);
    
    OFstream dataFile("wMass.dat");

    forAll(timeDirs, timeI)
    {
        word timeName = timeDirs[timeI].name();
	runTime.setTime(timeDirs[timeI], timeI);

        volScalarField wMass
	(
	    IOobject
	    (
	        "wallFilmMass",
		runTime.timeName(),
		mesh,
		IOobject::READ_IF_PRESENT,
		IOobject::NO_WRITE
	    ),
	    mesh,
	    dimensionedScalar("zero", dimMass, 0.0)
	);

        spray dieselSpray
        (
            U,
            rho,
            p,
            T,
            composition,
            gasProperties,
            thermo,
            g
        );

	scalar mTot = sum(wMass).value();
	scalar mLiq = dieselSpray.liquidMass();

	//reduce(mTot,sumOp<scalar>());
	//reduce(mLiq,sumOp<scalar>());

	dataFile << timeName << " " << mTot << " " << mLiq << endl;
	
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
