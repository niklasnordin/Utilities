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
    rhoPorousSimpleFoam

Description
    Steady-state solver for turbulent flow of compressible fluids with
    RANS turbulence modelling, implicit or explicit porosity treatment
    and MRF for HVAC and similar applications.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "dynamicRefineFvMesh.H"
#include "basicPsiThermo.H"
#include "RASModel.H"
#include "MRFZones.H"
#include "thermalPorousZones.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
    #include "createZones.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    bool hasChanged = false;
    dimensionedScalar Usmall("small", dimensionSet(0,0,-1,0,0,0,0), SMALL);
    dimensionedScalar epsSmall("small", dimensionSet(0,2,-3,0,0,0,0), SMALL);

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Indicators for refinement. Note: before runTime++
        // only for postprocessing reasons.
	const scalarField& V = mesh.V();
	scalarField length = pow(V, 1.0/3.0);
	const volScalarField& k = turbulence->k();

	volScalarField epsilon = turbulence->epsilon() + epsSmall;
	volScalarField lt = 0.09*pow(k, 1.5)/epsilon;
	lt.writeOpt() = IOobject::AUTO_WRITE;
        tmp<volScalarField> tscaleQ = mag(lt);

        volScalarField scaleQ
        (
            "scaleQ",
            tscaleQ()
        );

	forAll(scaleQ, i)
	{
            scaleQ[i] = length[i]/lt[i];
	}
	Info << "scaleQ min/max = " << min(scaleQ).value() << ", " << max(scaleQ).value() << endl;
        scaleQ.writeOpt() = IOobject::AUTO_WRITE;
        tscaleQ.clear();

        {
            // Make the fluxes absolute
            fvc::makeAbsolute(phi, rho, U);

            // Test : disable refinement for some cells
            PackedBoolList& protectedCell =
                refCast<dynamicRefineFvMesh>(mesh).protectedCell();

            if (protectedCell.empty())
            {
                protectedCell.setSize(mesh.nCells());
                protectedCell = 0;
            }

            // Flux estimate for introduced faces.
            volVectorField rhoU("rhoU", rho*U);

            // Do any mesh changes
            bool meshChanged = mesh.update();

            if (meshChanged)
            {
                hasChanged = true;
            }

            if (runTime.write() && hasChanged)
            {
                hasChanged = false;
            }

            // Make the fluxes relative to the mesh motion
            fvc::makeRelative(phi, rho, U);
        }

        // Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "hEqn.H"
            #include "pEqn.H"
        }

        turbulence->correct();

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
