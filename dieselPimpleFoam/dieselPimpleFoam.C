/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    dieselPimpleFoam

Description
    Solver for diesel spray and combustion, based on rhoPimpleFoam

    Transient solver for laminar or turbulent flow of compressible fluids
    for HVAC and similar applications.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "hsCombustionThermo.H"
#include "turbulenceModel.H"
#include "bound.H"
#include "spray.H"
#include "psiChemistryModel.H"
#include "chemistrySolver.H"
#include "multivariateScheme.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readGravitationalAcceleration.H"
    #include "readCombustionProperties.H"
    #include "createSpray.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensionedScalar muZero("muZero", dimensionSet(0,0,2,0,0,0,0), 0.0);
    dimensionedScalar epsZero("muZero", dimensionSet(0,2,-3,0,0,0,0), 1.0e-6);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "readPIMPLEControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (nOuterCorr != 1)
        {
            p.storePrevIter();
            rho.storePrevIter();
        }

        Info<< "Evolving Spray" << endl;

        dieselSpray.evolve();

        Info<< "Solving chemistry" << endl;

        chemistry.solve
        (
            runTime.value() - runTime.deltaT().value(),
            runTime.deltaT().value()
        );

        // turbulent time scale
        {

  	    volScalarField eps = max(epsZero, turbulence->epsilon());
            volScalarField tc = chemistry.tc();
	    volScalarField muTemp = max(muZero, turbulence->muEff()/rho/eps);
            volScalarField tk = Cmix*sqrt(muTemp);

            // Chalmers PaSR model
            kappa = (runTime.deltaT() + tc)/(runTime.deltaT()+tc+tk);
        }

        chemistrySh = kappa*chemistry.Sh()();

        volScalarField Cht = mag(kappa*chemistry.Sh()());
        Cht.max(1.0e-10);
        volScalarField tHR = mag(rho*hs)/Cht;

        #include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        for (int oCorr=0; oCorr<nOuterCorr; oCorr++)
        {
            #include "UEqn.H"
            #include "YEqn.H"
            #include "hsEqn.H"

            // --- PISO loop
            for (int corr=0; corr<nCorr; corr++)
            {
                #include "pEqn.H"
            }

            turbulence->correct();
        }

        #include "spraySummary.H"
        if (runTime.write())
        {
            chemistry.dQ()().write();
        }

	Info << "min/max T = " << min(T).value() << ", " << max(T).value() << endl;
	Info << "min/max U = " << min(mag(U)).value() << ", " << max(mag(U)).value() << endl;
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
