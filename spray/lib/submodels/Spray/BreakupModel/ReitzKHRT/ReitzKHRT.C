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

\*---------------------------------------------------------------------------*/

#include "ReitzKHRT.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ReitzKHRT<CloudType>::ReitzKHRT
(
    const dictionary& dict,
    CloudType& owner
)
:
    BreakupModel<CloudType>(owner),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    b0_(readScalar(coeffsDict_.lookup("B0"))),
    b1_(readScalar(coeffsDict_.lookup("B1"))),
    cTau_(readScalar(coeffsDict_.lookup("Ctau"))),
    cRT_(readScalar(coeffsDict_.lookup("CRT"))),
    msLimit_(readScalar(coeffsDict_.lookup("msLimit"))),
    weberLimit_(readScalar(coeffsDict_.lookup("WeberLimit")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ReitzKHRT<CloudType>::~ReitzKHRT()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ReitzKHRT<CloudType>::update
(
    const scalar& dt,
    scalar& d,
    scalar& tc,
    scalar& ms,
    scalar& nParticle,
    const scalar& rho,
    const scalar& mu,
    const scalar& sigma,
    const vector& U,
    const scalar& rhoc,
    const scalar& muc,
    const vector& Urel,
    const scalar& Urmag,
    const scalar& tMom,
    const scalar& averageParcelMass,
    scalar& dChild,
    scalar& massChild
) const
{

    bool addParcel = false;

  /*
    scalar r = 0.5*d;
    scalar d3 = pow(d, 3.0);
    scalar mass = nParticle*rho*mathematicalConstant::pi*d3/6.0;
    scalar Urmag = mag(Urel);
    scalar We = 0.5*rhoc*pow(Urmag, 2)*d/sigma;
    scalar Re = Urel*d/nuc;

    scalar weGas      = 0.5*rhoc*pow(Urmag, 2)*d/sigma;
    scalar weLiquid   = 0.5*rho*pow(Urmag, 2)*d/sigma;
    // Reitz is using radius instead of diameter for Re-number
    scalar reLiquid   = rho*Urmag*r/mu;
    scalar ohnesorge  = sqrt(weLiquid)/(reLiquid + VSMALL);
    scalar taylor     = ohnesorge*sqrt(weGas);

    vector acceleration = Urel/tMom;
    vector trajectory = U/mag(U);
    scalar gt = (g + acceleration) & trajectory;

    // frequency of the fastest growing KH-wave
    scalar omegaKH =
        (0.34 + 0.38*pow(weGas, 1.5))
       /((1.0 + ohnesorge)*(1.0 + 1.4*pow(taylor, 0.6)))
       *sqrt(sigma/(rho*pow(r, 3.0)));

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
    scalar helpVariable = mag(gt*(rho - rhoc));
    scalar omegaRT = sqrt
    (
        2.0*pow(helpVariable, 1.5)
       /(3.0*sqrt(3.0*sigma)*(rhoc + rho))
    );

    // RT wave number
    scalar KRT = sqrt(helpVariable/(3.0*sigma + VSMALL));

    // wavelength of the fastest growing RT frequency
    scalar lambdaRT = 2.0*mathematicalConstant::pi*cRT_/(KRT + VSMALL);

    // if lambdaRT < diameter, then RT waves are growing on the surface
    // and we start to keep track of how long they have been growing
    if ((tc > 0) || (lambdaRT < d) )
    {
        tc += dt;
    }

    // characteristic RT breakup time
    scalar tauRT = cTau_/(omegaRT + VSMALL);

    // check if we have RT breakup
    if ((tc > tauRT) && (lambdaRT < d) )
    {
        // the RT breakup creates diameter/lambdaRT new droplets
        tc = -GREAT;
        scalar nDrops = d/lambdaRT;
        d = cbrt(d3/nDrops);
    }
    // otherwise check for KH breakup
    else if (dc < d)
    {
        // no breakup below Weber = 12
        if (weGas > weberLimit_)
        {
            scalar fraction = dt/tauKH;

            // reduce the diameter according to the rate-equation
            d = (fraction*dc + d)/(1.0 + fraction);

            scalar ms0 = rho*pow3(dc)*mathematicalConstant::pi/6.0;
            ms += ms0;

            if (ms/averageParcelMass > msLimit_)
            {
		addParcel = true;
                // set the initial ms value to -GREAT. This prevents
                // new droplets from being formed from the child droplet
                // from the KH instability

                // mass of stripped child parcel
                // Prevent child parcel from taking too much mass
                scalar mc = min(ms, 0.5*mass);

		// reduce the parcel mass by reducing nParticle
		scalar mass1 = mass - mc;
		scalar massDrop = rho*mathematicalConstant::pi*d3/6.0;
		nParticle = mass/massDrop;
                ms = 0.0;
            }
        }
    }
  */
    return addParcel;
}


// ************************************************************************* //
