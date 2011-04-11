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
void Foam::ReitzKHRT<CloudType>::breakup
(
    const scalar& dt,
    scalar& d,
    const scalar& rho,
    const scalar& mu,
    const scalar& sigma,
    const scalar& rhoc,
    const scalar& muc,
    const scalar& Urel
) const
{
  /*
    scalar r = 0.5*d;
    scalar d3 = pow(d, 3.0);

    scalar We = 0.5*rhoc*pow(Urel, 2)*d/sigma;
    scalar Re = Urel*d/nuc;

    scalar weGas      = 0.5*rhoc*pow(Urel, 2)*d/sigma;
    scalar weLiquid   = 0.5*rho*pow(Urel, 2)*d/sigma;
    // Reitz is using radius instead of diameter for Re-number
    scalar reLiquid   = rho*Urel*r/mu;
    scalar ohnesorge  = sqrt(weLiquid)/(reLiquid + VSMALL);
    scalar taylor     = ohnesorge*sqrt(weGas);

    vector acceleration = p.Urel(vel)/p.tMom();
    vector trajectory = p.U()/mag(p.U());
    scalar gt = (g + acceleration) & trajectory;

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
        scalar multiplier = p.d()/lambdaRT;
        scalar nDrops = multiplier*Np;
        p.d() = cbrt(semiMass/nDrops);
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

            scalar ms = rhoLiquid*Np*pow3(dc)*mathematicalConstant::pi/6.0;
            p.ms() += ms;

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

                // mass of stripped child parcel
                scalar mc = p.ms();
                // Prevent child parcel from taking too much mass
                if (mc > 0.5*p.m())
                {
                    mc = 0.5*p.m();
                }

                spray_.addParticle
                (
                    new parcel
                    (
                        spray_,
                        p.position(),
                        p.cell(),
                        p.n(),
                        dc,
                        p.T(),
                        mc,
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

                p.m() -= mc;
                p.ms() = 0.0;
            }
        }
    }
  */
}


// ************************************************************************* //
