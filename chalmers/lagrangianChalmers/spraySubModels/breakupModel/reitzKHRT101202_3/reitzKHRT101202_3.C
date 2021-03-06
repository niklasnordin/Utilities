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

\*---------------------------------------------------------------------------*/

#include "reitzKHRT101202_3.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(reitzKHRT101202_3, 0);

addToRunTimeSelectionTable
(
    breakupModel,
    reitzKHRT101202_3,
    dictionary
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
reitzKHRT101202_3::reitzKHRT101202_3
(
    const dictionary& dict,
    spray& sm
)
:
    breakupModel(dict, sm),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    g_(sm.g()),
    b0_(readScalar(coeffsDict_.lookup("B0"))),
    b1_(readScalar(coeffsDict_.lookup("B1"))),
    cTau_(readScalar(coeffsDict_.lookup("Ctau"))),
    cRT_(readScalar(coeffsDict_.lookup("CRT"))),
    msLimit_(readScalar(coeffsDict_.lookup("msLimit"))),
    weberLimit_(readScalar(coeffsDict_.lookup("WeberLimit")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

reitzKHRT101202_3::~reitzKHRT101202_3()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void reitzKHRT101202_3::breakupParcel
(
    parcel& p,
    const scalar deltaT,
    const vector& vel,
    const liquidMixture& fuels
) const
{

    label celli = p.cell();
    scalar T = p.T();
    scalar r = 0.5*p.d();
    scalar pc = spray_.p()[celli];

    scalar sigma = fuels.sigma(pc, T, p.X());
    scalar rhoLiquid = fuels.rho(pc, T, p.X());
    scalar muLiquid = fuels.mu(pc, T, p.X());
    scalar rhoGas = spray_.rho()[celli];
    scalar Np = p.N(rhoLiquid);
    scalar semiMass = Np*pow(p.d(), 3.0);

    scalar weGas      = p.We(vel, rhoGas, sigma);
    scalar weLiquid   = p.We(vel, rhoLiquid, sigma);
    // correct the Reynolds number. Reitz is using radius instead of diameter
    scalar reLiquid   = 0.5*p.Re(rhoLiquid, vel, muLiquid);
    scalar ohnesorge  = sqrt(weLiquid)/(reLiquid + VSMALL);
    scalar taylor     = ohnesorge*sqrt(weGas);

    vector acceleration = p.Urel(vel)/p.tMom();
    vector trajectory = p.U()/mag(p.U());
    scalar gt = (g_ + acceleration) & trajectory;

    // frequency of the fastest growing KH-wave
    scalar omegaKH =
        (0.34 + 0.38*pow(weGas, 1.5))
       /((1 + ohnesorge)*(1 + 1.4*pow(taylor, 0.6)))
       *sqrt(sigma/(rhoLiquid*pow(r, 3)));

    // corresponding KH wave-length.
    scalar lambdaKH =
        9.02
       *r
       *(1.0 + 0.45*sqrt(ohnesorge))
       *(1.0 + 0.4*pow(taylor, 0.7))
       /pow(1.0 + 0.865*pow(weGas, 1.67), 0.6);

    // characteristic Kelvin-Helmholtz breakup time
    scalar tauKH = 3.726*b1_*r/(omegaKH*lambdaKH);

    // stable KH diameter
    scalar dc = 2.0*b0_*lambdaKH;

    // the frequency of the fastest growing RT wavelength.
    scalar helpVariable = mag(gt*(rhoLiquid - rhoGas));
    scalar omegaRT = sqrt
    (
        2.0*pow(helpVariable, 1.5)
       /(3.0*sqrt(3.0*sigma)*(rhoGas + rhoLiquid))
    );

    // RT wave number
    scalar KRT = sqrt(helpVariable/(3.0*sigma + VSMALL));

    // wavelength of the fastest growing RT frequency
    scalar lambdaRT = 2.0*mathematicalConstant::pi*cRT_/(KRT + VSMALL);

    // if lambdaRT < diameter, then RT waves are growing on the surface
    // and we start to keep track of how long they have been growing
    if ((p.ct() > 0) || (lambdaRT < p.d()))
    {
        p.ct() += deltaT;
    }

    // characteristic RT breakup time
    scalar tauRT = cTau_/(omegaRT + VSMALL);

    // check if we have RT breakup
    if ((p.ct() > tauRT) && (lambdaRT < p.d()))
    {
        // the RT breakup creates diameter/lambdaRT new droplets
        p.ct() = -GREAT;

//--------------------------------AL_____101201------------------------------//
      scalar multiplier = p.d()/lambdaRT;
      scalar nDrops = multiplier*Np;
      p.d() = cbrt(semiMass/nDrops);
//------------------------------------END------------------------------------//
    }
    // otherwise check for KH breakup
    else if (dc < p.d())
    {
        // no breakup below Weber = 12
        if (weGas > weberLimit_)
        {

            label injector = label(p.injector());
            scalar fraction = deltaT/tauKH;

            // reduce the diameter according to the rate-equation
            p.d() = (fraction*dc + p.d())/(1.0 + fraction);

//--------------------------------AL_____101201------------------------------//
//          2. Correct evaluation of the stripped mass
//             In this case, msLimit should be equal to 0.03
//            scalar ms = rhoLiquid*Np*pow3(dc)*mathematicalConstant::pi/6.0;
//            p.ms() += ms;		
	    p.ms() = p.m()*(1.0-pow(p.d()/p.d0(),3.0));
//------------------------------------END------------------------------------//

            // Total number of parcels for the whole injection event
            label nParcels =
                spray_.injectors()[injector].properties()->nParcelsToInject
                (
                    spray_.injectors()[injector].properties()->tsoi(),
                    spray_.injectors()[injector].properties()->teoi()
                );

            scalar averageParcelMass =
                spray_.injectors()[injector].properties()->mass()/nParcels;

            if (p.ms()/averageParcelMass > msLimit_)
            {
                // set the initial ms value to -GREAT. This prevents
                // new droplets from being formed from the child droplet
                // from the KH instability

//--------------------------------AL_____101201------------------------------//
//          2. Correct evaluation of the number of child droplets and the diameter of parcel droplets after breaukp
//		Solution of cubic equation for the diameter of the parent drops after breakup, see Eq. 18 in Patterson & Reitz, SAE 980131
		bool br3 = true;
		scalar ae3 = 1.;
		scalar be3 = -dc;
		scalar ce3 = 0.;
		scalar de3 = p.d()*p.d()*(dc-p.d());
		scalar qe3 = pow3(be3/(3.*ae3))-be3*ce3/(6.*ae3*ae3)+de3/(2.*ae3);
		scalar pe3 = (3.*ae3*ce3-be3*be3)/(9.*ae3*ae3);
		scalar D3 = qe3*qe3 + pe3*pe3*pe3; 
		if (D3 < 0) br3 = false;
		if (br3) 
		{
			D3 = sqrt(D3);
			scalar ue3 = cbrt(-qe3+D3);	
                        scalar ve3 = cbrt(-qe3-D3);
			scalar dParenDrops = ue3 + ve3 - be3/3.;
                        scalar mc = Np*(pow3(p.d())-pow3(dParenDrops));
			scalar nChildDrops = mc/pow3(dc);
			if (nChildDrops >= Np)
			{
		                p.m() -= mc;
                                p.d() = dParenDrops;
                		p.ms() = 0.;

		                spray_.addParticle
                		(
                    			new parcel
                    			(
//	Add the following lines
//      					dc,
//						1
//	in order to set p.d0() and pIndKH()
                        			spray_,
                        			p.position(),
                        			p.cell(),
                        			p.n(),
                        			dc,
                                                dc,
                        			p.T(),
                        			mc,
                                                1,
                        			0.0,
                        			0.0,
                        			0.0,
                        			-GREAT,
                        			p.tTurb(),
                        			0.0,
                        			p.injector(),
                        			p.U(),
                        			p.Uturb(),
                        			p.X(),
                        			p.fuelNames()
                    			)		
                		);
			}
		}

                // mass of stripped child parcel
//              scalar mc = p.ms();
                // Prevent child parcel from taking too much mass
//              if (mc > 0.5*p.m())
//              {
//                  mc = 0.5*p.m();
//              }
//
//              spray_.addParticle
//              (
//                  new parcel
//                  (
//                      spray_,
//                      p.position(),
//                      p.cell(),
//                      p.n(),
//                      dc,
//                      p.T(),
//                      mc,
//                      0.0,
//                      0.0,
//                      0.0,
//                      -GREAT,
//                      p.tTurb(),
//                      0.0,
//                      p.injector(),
//                      p.U(),
//                      p.Uturb(),
//                      p.X(),
//                      p.fuelNames()
//                  )
//              );
//
//              p.m() -= mc;
//              p.ms() = 0.0;
//------------------------------------END------------------------------------//
            }
        }
    }
//--------------------------------AL_____101202------------------------------//
//	3. Case of larger drops after breakup (Reitz, Atomization & Spray Technology 3 (1987) 309-337, p.322)
//         pIndKH() should be introduced
    else if (p.pIndKH() == 0)
    {
	scalar lengthScale = min(lambdaKH,2.0*mathematicalConstant::pi*mag(p.Urel(vel))/omegaKH);
	scalar diameterLargerDrop = cbrt(1.5*p.d()*p.d()*lengthScale);
	p.d() = diameterLargerDrop;
        p.ms() = 0.;
	p.pIndKH() = 1;
    }
//------------------------------------END------------------------------------//
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
