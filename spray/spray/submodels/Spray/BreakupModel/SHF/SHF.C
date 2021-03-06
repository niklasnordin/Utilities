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

#include "SHF.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::SHF<CloudType>::SHF
(
    const dictionary& dict,
    CloudType& owner
)
:
    BreakupModel<CloudType>(dict, owner, typeName),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    weCorrCoeff_(readScalar(coeffsDict_.lookup("weCorrCoeff"))),
    weBuCrit_(readScalar(coeffsDict_.lookup("weBuCrit"))),
    weBuBag_(readScalar(coeffsDict_.lookup("weBuBag"))),
    weBuMM_(readScalar(coeffsDict_.lookup("weBuMM"))),
    ohnCoeffCrit_(readScalar(coeffsDict_.lookup("ohnCoeffCrit"))),
    ohnCoeffBag_(readScalar(coeffsDict_.lookup("ohnCoeffBag"))),
    ohnCoeffMM_(readScalar(coeffsDict_.lookup("ohnCoeffMM"))),
    ohnExpCrit_(readScalar(coeffsDict_.lookup("ohnExpCrit"))),
    ohnExpBag_(readScalar(coeffsDict_.lookup("ohnExpBag"))),
    ohnExpMM_(readScalar(coeffsDict_.lookup("ohnExpMM"))),
    cInit_(readScalar(coeffsDict_.lookup("Cinit"))),
    c1_(readScalar(coeffsDict_.lookup("C1"))),
    c2_(readScalar(coeffsDict_.lookup("C2"))),
    c3_(readScalar(coeffsDict_.lookup("C3"))),
    cExp1_(readScalar(coeffsDict_.lookup("Cexp1"))),
    cExp2_(readScalar(coeffsDict_.lookup("Cexp2"))),
    cExp3_(readScalar(coeffsDict_.lookup("Cexp3"))),
    weConst_(readScalar(coeffsDict_.lookup("Weconst"))),
    weCrit1_(readScalar(coeffsDict_.lookup("Wecrit1"))),
    weCrit2_(readScalar(coeffsDict_.lookup("Wecrit2"))),
    coeffD_(readScalar(coeffsDict_.lookup("CoeffD"))),
    onExpD_(readScalar(coeffsDict_.lookup("OnExpD"))),
    weExpD_(readScalar(coeffsDict_.lookup("WeExpD"))),
    mu_(readScalar(coeffsDict_.lookup("mu"))),
    sigma_(readScalar(coeffsDict_.lookup("sigma"))),
    d32Coeff_(readScalar(coeffsDict_.lookup("d32Coeff"))),
    cDmaxBM_(readScalar(coeffsDict_.lookup("cDmaxBM"))),
    cDmaxS_(readScalar(coeffsDict_.lookup("cDmaxS"))),
    corePerc_(readScalar(coeffsDict_.lookup("corePerc")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::SHF<CloudType>::~SHF()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::SHF<CloudType>::update
(
    const scalar& dt,
    const vector& g,
    scalar& d,
    scalar& tc,
    scalar& ms,
    scalar& nParticle,
    scalar& KHindex,
    scalar& y,
    scalar& yDot,
    const scalar& d0,
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
    scalar& massChild,
    Random& rndGen
) const
{
    bool addChild = false;

    scalar d03 = pow(d, 3);
    scalar rhopi6 = rho*mathematicalConstant::pi/6.0;
    scalar mass0 = nParticle*rhopi6*d03;
    scalar mass = mass0;

    scalar weGas      = 0.5*rhoc*pow(Urmag, 2.0)*d/sigma;
    scalar weLiquid   = 0.5*rho*pow(Urmag, 2.0)*d/sigma;

    // correct the Reynolds number. Reitz is using radius instead of diameter
    scalar reLiquid   = 0.5*Urmag*d/mu;
    scalar ohnesorge  = sqrt(weLiquid)/(reLiquid + VSMALL);

    vector acceleration = Urel/tMom;
    vector trajectory = U/mag(U);

    scalar weGasCorr = weGas/(1.0 + weCorrCoeff_ * ohnesorge);

    // droplet deformation characteristic time

    scalar tChar = d/Urmag*sqrt(rho/rhoc);

    scalar tFirst = cInit_ * tChar;

    scalar tSecond = 0;
    scalar tCharSecond = 0;

    bool bag = false;
    bool multimode = false;
    bool shear = false;
    bool success = false;


    //  updating the droplet characteristic time
    tc += dt;

    if(weGas > weConst_)
    {
        if(weGas < weCrit1_)
        {
            tCharSecond = c1_*pow((weGas - weConst_),cExp1_);
        }
        else if(weGas >= weCrit1_ && weGas <= weCrit2_)
        {
            tCharSecond = c2_*pow((weGas - weConst_),cExp2_);
        }
        else
        {
            tCharSecond = c3_*pow((weGas - weConst_),cExp3_);
        }
    }

    scalar weC = weBuCrit_*(1.0+ohnCoeffCrit_*pow(ohnesorge,ohnExpCrit_));
    scalar weB = weBuBag_*(1.0+ohnCoeffBag_*pow(ohnesorge, ohnExpBag_));
    scalar weMM = weBuMM_*(1.0+ohnCoeffMM_*pow(ohnesorge, ohnExpMM_));

    if(weGas > weC && weGas < weB)
    {
        bag = true;
    }

    if(weGas >= weB && weGas <= weMM)
    {
        multimode = true;
    }

    if(weGas > weMM)
    {
        shear = true;
    }

    tSecond = tCharSecond * tChar;

    scalar tBreakUP = tFirst + tSecond;
    if(tc > tBreakUP)
    {

        scalar d32 = coeffD_*d*pow(ohnesorge,onExpD_)*pow(weGasCorr,weExpD_);

        if(bag || multimode)
        {

            scalar d05 = d32Coeff_ * d32;

            scalar x = 0.0;
            scalar yGuess = 0.0;
            scalar dGuess = 0.0;

            while(!success)
            {
                x = cDmaxBM_*rndGen.scalar01();
                dGuess = sqr(x)*d05;
                yGuess = rndGen.scalar01();

                scalar p = x/(2.0*sqrt(2.0*mathematicalConstant::pi)*sigma_)*exp(-0.5*sqr((x-mu_)/sigma_));

                if (yGuess < p)
                {
                    success = true;
                }
            }

            d = dGuess;
            tc = 0.0;
        }

        if(shear)
        {
            scalar dC = weConst_*sigma/(rhoc*sqr(Urmag));
            scalar d32Red = 4.0*(d32 * dC)/(5.0 * dC - d32);

            scalar d05 = d32Coeff_ * d32Red;

            scalar x = 0.0;
            scalar yGuess = 0.0;
            scalar dGuess = 0.0;

            while(!success)
            {

                x = cDmaxS_*rndGen.scalar01();
                dGuess = sqr(x)*d05;
                yGuess = rndGen.scalar01();
                
                scalar p = x/(2.0*sqrt(2.0*mathematicalConstant::pi)*sigma_)*exp(-0.5*sqr((x-mu_)/sigma_));
                
                if (yGuess<p)
                {
                    success = true;
                }
            }
            
            d = dC;
            dChild = dGuess;
            massChild = corePerc_ * mass0;
            mass -= massChild;

            addChild = true;
            // reset timer
            tc = 0.0;
        }

        // correct nParticle to conserve mass
        nParticle = mass/(rhopi6*pow(d, 3.0));
    }

    return addChild;
}


// ************************************************************************* //
