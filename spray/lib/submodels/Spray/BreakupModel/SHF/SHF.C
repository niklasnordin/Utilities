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
    BreakupModel<CloudType>(owner),
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
    // Do nothing
    return false;
}


// ************************************************************************* //
