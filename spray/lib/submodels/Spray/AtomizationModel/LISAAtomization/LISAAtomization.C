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

#include "LISAAtomization.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::LISAAtomization<CloudType>::LISAAtomization
(
    const dictionary& dict,
    CloudType& owner
)
:
    AtomizationModel<CloudType>(owner),
    coeffsDict_(dict.subDict(typeName+"Coeffs")),
    Cl_(readScalar(coeffsDict_.lookup("Cl"))),
    cTau_(readScalar(coeffsDict_.lookup("cTau"))),
    Q_(readScalar(coeffsDict_.lookup("Q"))),
    lisaExp_(readScalar(coeffsDict_.lookup("lisaExp"))),
    injectorDirection_(coeffsDict_.lookup("injectorDirection"))
{
    // NN. Would be good if this could be picked up from the injector
    injectorDirection_ /= mag(injectorDirection_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::LISAAtomization<CloudType>::~LISAAtomization()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::LISAAtomization<CloudType>::initLiquidCore() const
{
    return 1.0;
}

template<class CloudType>
Foam::scalar Foam::LISAAtomization<CloudType>::Taverage
(
    const scalar& Tliq,
    const scalar& Tc
) const
{
    return (2.0*Tliq + Tc)/3.0;
}


template<class CloudType>
void Foam::LISAAtomization<CloudType>::update
(
    const scalar& dt,
    scalar& d,
    scalar& liquidCore,
    scalar& tc,
    const scalar& rho,
    const scalar& mu,
    const scalar& sigma,
    const scalar& massflowRate,
    const scalar& rhoAv,
    const scalar& Urel,
    const vector& pos,
    const vector& injectionPos,
    const scalar& pAmbient,
    const scalar& chi,
    Random& rndGen
) const
{

    scalar tau = 0.0;
    scalar dL = 0.0;
    scalar k = 0.0;

    // update atomization characteristic time
    tc += dt;

    scalar We = 0.5*rhoAv*pow(Urel, 2)*d/sigma;
    scalar nu = mu/rho;

    scalar Q = rhoAv/rho;

    vector diff = pos - injectionPos;
    scalar pWalk = mag(diff);

    scalar h = diff & injectorDirection_;
    scalar delta = sqrt(sqr(pWalk) - sqr(h));

    scalar hSheet = massflowRate/(mathematicalConstant::pi*delta*rho*Urel);

    // update drop diameter
    d = min(d, hSheet);


    if(We > 27.0/16.0)
    {

        scalar kPos = 0.0;
        scalar kNeg = Q*pow(Urel, 2.0)*rho/sigma;

        scalar derivativePos = sqrt
        (
            Q*pow(Urel, 2.0)
        );

        scalar derivativeNeg =
        (
            8.0*pow(nu, 2.0)*pow(kNeg, 3.0)
            + Q*pow(Urel, 2.0)*kNeg
            - 3.0*sigma/2.0/rho*pow(kNeg, 2.0)
        )
        /
        sqrt
        (
            4.0*pow(nu, 2.0)*pow(kNeg, 4.0)
            + Q*pow(Urel, 2.0)*pow(kNeg, 2.0)
            - sigma*pow(kNeg, 3.0)/rho
        )
        -
        4.0*nu*kNeg;

        scalar kOld = 0.0;

        for(label i=0; i<40; i++)
        {

            k = kPos - (derivativePos/((derivativeNeg-derivativePos)/(kNeg-kPos)));

            scalar derivativek =
            (
                8.0*pow(nu, 2.0)*pow(k, 3.0)
                + Q*pow(Urel, 2.0)*k
                - 3.0*sigma/2.0/rho*pow(k, 2.0)
            )
            /
            sqrt
            (
                4.0*pow(nu, 2.0)*pow(k, 4.0)
                + Q*pow(Urel, 2.0)*pow(k, 2.0)
                - sigma*pow(k, 3.0)/rho
            )
            -
            4.0*nu*k;

            if(derivativek > 0)
            {
                derivativePos = derivativek;
                kPos = k;
            }
            else
            {
                derivativeNeg = derivativek;
                kNeg = k;
            }

            if(mag(k-kOld)/k < 1e-4)
            {
                break;
            }

            kOld = k;

        }

        scalar omegaS =
        - 2.0 * nu * pow(k, 2.0)
        + sqrt
        (
                4.0*pow(nu, 2.0)*pow(k, 4.0)
            +   Q*pow(Urel, 2.0)*pow(k, 2.0)
            -   sigma*pow(k, 3.0)/rho
        );

        tau = cTau_/omegaS;

        dL = sqrt(8.0*d/k);

    }
    else
    {

        k =
        rhoAv*pow(Urel, 2.0)
        /
        2.0*sigma;

        scalar J = pWalk*d/2.0;

        tau = pow(3.0*cTau_,2.0/3.0)*cbrt(J*sigma/(sqr(Q)*pow(Urel,4.0)*rho));

        dL = sqrt(4.0*d/k);
    }

    scalar kL =
        1.0
        /
        (
            dL *
            pow(0.5 + 1.5 * mu/pow((rho*sigma*dL), 0.5), 0.5)
        );

    scalar dD = cbrt(3.0*mathematicalConstant::pi*pow(dL, 2.0)/kL);

    scalar atmPressure = 1.0e+5;

    scalar pRatio = pAmbient/atmPressure;

    dD = dD*pow(pRatio, lisaExp_);

    scalar pExp = 0.135;

    /*
//  modifications to take account of the flash boiling on primary breakUp

    scalar chi = 0.0;

    label Nf = fuels.components().size();

    scalar Td = p.T();

    for(label i = 0; i < Nf ; i++)
    {

        if(fuels.properties()[i].pv(spray_.ambientPressure(), Td) >= 0.999*spray_.ambientPressure())
        {

//          The fuel is boiling.....
//          Calculation of the boiling temperature

            scalar tBoilingSurface = Td;

            label Niter = 200;

            for(label k=0; k< Niter ; k++)
            {
                scalar pBoil = fuels.properties()[i].pv(pressure, tBoilingSurface);

                if(pBoil > pressure)
                {
                    tBoilingSurface = tBoilingSurface - (Td-temperature)/Niter;
                }
                else
                {
                    break;
                }
            }

            scalar hl = fuels.properties()[i].hl(spray_.ambientPressure(), tBoilingSurface);
            scalar iTp = fuels.properties()[i].h(spray_.ambientPressure(), Td) - spray_.ambientPressure()/fuels.properties()[i].rho(spray_.ambientPressure(), Td);
            scalar iTb = fuels.properties()[i].h(spray_.ambientPressure(), tBoilingSurface) - spray_.ambientPressure()/fuels.properties()[i].rho(spray_.ambientPressure(), tBoilingSurface);

            chi += p.X()[i]*(iTp-iTb)/hl;

        }
    }

    //  bounding chi

    chi = max(chi, 0.0);
    chi = min(chi, 1.0);
    */

    //  modifing dD to take account of flash boiling

    dD = dD*(1.0 - chi*pow(pRatio, -pExp));

    scalar lBU = Cl_ * mag(Urel)*tau;

    if(pWalk > lBU)
    {

        liquidCore = 0.0;

//      calculate the new diameter with a Rosin Rammler distribution

        scalar minValue = min(d, dD/10.0);

        scalar maxValue = dD;

        if(maxValue - minValue < SMALL)
        {
            minValue = d/10.0;
        }

        scalar range = maxValue - minValue;

        scalar y = 0;
        scalar x = 0;

        bool success = false;

        while(!success)
        {

            x = minValue + range*rndGen.scalar01();
            y = rndGen.scalar01();

            scalar p = 0.0;

            scalar nExp = 1;

            scalar xx = pow(x/dD, nExp);

            p = xx*exp(-xx);

            if (y<p)
            {
                success = true;
            }

        }

        //  New droplet properties
        d = x;
        tc = 0.0;

    }

}

// ************************************************************************* //
