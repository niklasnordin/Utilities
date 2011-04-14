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
    ANY WARRANTY; without even the implied warranty of MERCHANTAB_RR3ILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "TAB_RR3.H"
#include "addToRunTimeSelectionTable.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(TAB_RR3, 0);

addToRunTimeSelectionTable
(
    breakupModel,
    TAB_RR3,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
TAB_RR3::TAB_RR3
(
    const dictionary& dict,
    spray& sm
)
:
    breakupModel(dict, sm),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    Cmu_(readScalar(coeffsDict_.lookup("Cmu"))),
    Comega_(readScalar(coeffsDict_.lookup("Comega"))),
    WeCrit_(readScalar(coeffsDict_.lookup("WeCrit"))),
    K_(readScalar(coeffsDict_.lookup("K")))
{

//--------------------------------AL_____101110------------------------------//
/*
    // calculate the inverse function of the Rossin-Rammler Distribution
//    const scalar xx0 = 12.0;
//    const scalar rrd100 = 1.0/(1.0-exp(-xx0)*(1+xx0+pow(xx0, 2)/2+pow(xx0, 3)/6));

    for(label n=0; n<100; n++)
    {
        scalar xx = 0.12*(n+1);
        rrd_[n] = (1-exp(-xx)*(1 + xx + pow(xx, 2)/2 + pow(xx, 3)/6))*rrd100;
    }
*/
//    const scalar xx0 = 18.0;
//    const scalar rrd100 = 1.0/(1.0-exp(-xx0)*(1+xx0+pow(xx0, 2)/2+pow(xx0, 3)/6));
//    for(label n=0; n<100; n++)
//    {
//        scalar xx = 0.12*(n+1);
//        rrd_[n] = (1-exp(-xx)*(1 + xx + pow(xx, 2)/2 + pow(xx, 3)/6))*rrd100;
//    }
//-----------------------------------END-------------------------------------//
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

TAB_RR3::~TAB_RR3()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void TAB_RR3::breakupParcel
(
    parcel& p,
    const scalar deltaT,
    const vector& Ug,
    const liquidMixture& fuels
) const
{

    scalar T = p.T();
    scalar pc = spray_.p()[p.cell()];
    scalar r = 0.5*p.d();
    scalar r2 = r*r;
    scalar r3 = r*r2;

    scalar rho = fuels.rho(pc, T, p.X());
    scalar sigma = fuels.sigma(pc, T, p.X());
    scalar mu = fuels.mu(pc, T, p.X());

    // inverse of characteristic viscous damping time
    scalar rtd = 0.5*Cmu_*mu/(rho*r2);
    
    // oscillation frequency (squared)
    scalar omega2 = Comega_*sigma/(rho*r3) - rtd*rtd;

    if (omega2 > 0)
    {
        scalar omega = sqrt(omega2);
        scalar rhog = spray_.rho()[p.cell()];
        scalar We = p.We(Ug, rhog, sigma);
        scalar Wetmp = We/WeCrit_;

        scalar y1 = p.dev() - Wetmp;
        scalar y2 = p.ddev()/omega;

        scalar a = sqrt(y1*y1 + y2*y2);

//-----------------------------------CH 101111----------------------------------------//
        scalar cosomt = 0.0;
        scalar sinomt = 0.0;
        scalar expttd = 0.0;
//--------------------------------------END------------------------------------------//

        // scotty we may have break-up
        if (a+Wetmp > 1.0)
        {
            scalar phic = y1/a;

            // constrain phic within -1 to 1
            phic = max(min(phic, 1), -1);

            scalar phit = acos(phic);
            scalar phi = phit;
            scalar quad = -y2/a;
            if (quad < 0)
            {
                phi = 2*mathematicalConstant::pi - phit;
            }
            
            scalar tb = 0;
          
//-----------------------------------CH 101111----------------------------------------//
            bool breakupYes = false;
//--------------------------------------END------------------------------------------//

            if (mag(p.dev()) < 1.0)
            {
                scalar coste = 1.0;
                if
                (
                    (Wetmp - a < -1)
                 && (p.ddev() < 0)
                )
                {
                    coste = -1.0;
                }
                
                scalar theta = acos((coste-Wetmp)/a);
                
                if (theta < phi)
                {
                    if (2*mathematicalConstant::pi-theta >= phi)
                    {
                        theta = -theta;
                    }
                    theta += 2*mathematicalConstant::pi;
                }
                tb = (theta-phi)/omega;

                // breakup occurs
                if (deltaT >= tb)
                {
                    p.dev() = 1.0;
                    p.ddev() = -a*omega*sin(omega*tb + phi);
                    breakupYes = true;
                }
//--------------------------------AL_____101109------------------------------//
                else
                {
                    breakupYes = false;
                    cosomt = cos(omega*deltaT);
                    sinomt = sin(omega*deltaT);
                    expttd = exp(-deltaT*rtd);
                    y2 = (p.ddev()+y1*rtd)/omega;
                    p.dev() = Wetmp + expttd*(y1*cosomt+y2*sinomt);
                    p.ddev() = (Wetmp-p.dev())*rtd+expttd*omega*(y2*cosomt-y1*sinomt);
                }
            }
            else
            {
            breakupYes = true;
            }
            // update droplet size
//            if (deltaT > tb)
            if (breakupYes)
//-----------------------------------END-------------------------------------//
            {
                scalar rs = r/
                (
                    1 
                  + (8.0/20.0)*pow(p.dev(), 2)*K_
                  + rho*r3/sigma*pow(p.ddev(), 2)*(6.0*K_-5.0)/120.0
                );

//--------------------------------AL_____101012------------------------------//
// Calculation of the mean radius based on SMR rs. Coefficient factorGamma depends on nExp.
//              scalar factorGamma = 0.75*sqrt(mathematicalConstant::pi);       //nExp=2
                scalar factorGamma = 1.;
                scalar delta = rs/factorGamma;
//              scalar delta = rs/(0.75*sqrt(mathematicalConstant::pi));
//-----------------------------------END-------------------------------------//
//                Info << "deltaDTab = " << delta << endl;
//--------------------------------AL_____101104------------------------------//
                scalar minValue = min (p.d()/2.0, 0.04*rs);
//              scalar minValue = min (p.d()/2.0, delta/40.0);
                scalar maxValue = rs*4.0;
//              scalar maxValue = delta*4.0;
//-----------------------------------END-------------------------------------//
                scalar range = maxValue - minValue;

                if(maxValue - minValue < SMALL)
                {
//--------------------------------AL_____101015------------------------------//
                     minValue = p.d()/20.0;
                     maxValue = p.d();
//                   minValue = maxValue/10.0;
//-----------------------------------END-------------------------------------//
                }
                    
                scalar nExp = 3.5;
                scalar rrd_[100];
//--------------------------------AL_____101012------------------------------//
                scalar probFactorMin = exp(-pow(minValue/delta,nExp));
                scalar probFactorMax = exp(-pow(maxValue/delta,nExp));
                scalar probFactor = 1./(probFactorMin - probFactorMax);
//-----------------------------------END-------------------------------------//
                for(label n=0; n<100; n++)
                {
                    scalar xx = minValue + range*n/100;
//-------------------------------A-L_____101012------------------------------//
//                  rrd_[n] = 1 - exp(-pow(xx/delta,nExp));
                    rrd_[n] = (probFactorMin - exp(-pow(xx/delta,nExp)))*probFactor;
//-----------------------------------END-------------------------------------//
                }
                
                label n = 0;
                bool found = false;
                scalar random = rndGen_.scalar01();
                scalar rNew = 0;
                while (!found && (n<100))
                {
                    if (rrd_[n]>random)
                    {
                        found = true;
                    }
                    n++;

                }
//--------------------------------AL_____101012------------------------------//
//              rNew = minValue + range*n/100;
                rNew = minValue + range*(n-0.5)/100.0;
//------------------------------------END------------------------------------//

                if (rNew > r) rNew=r;
                    p.d() = 2*rNew;
                    p.dev() = 0;
                    p.ddev() = 0;

            }

        }

        else
        {
        cosomt = cos(omega*deltaT);
        sinomt = sin(omega*deltaT);
        expttd = exp(-deltaT*rtd);
        y2 = (p.ddev()+y1*rtd)/omega;
        p.dev() = Wetmp + expttd*(y1*cosomt+y2*sinomt);
        p.ddev() = (Wetmp-p.dev())*rtd+expttd*omega*(y2*cosomt-y1*sinomt);
        }
//------------------------------------END------------------------------------//       
    }
    else
    {
        // reset droplet distortion parameters
        p.dev() = 0;
        p.ddev() = 0;
    }

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
