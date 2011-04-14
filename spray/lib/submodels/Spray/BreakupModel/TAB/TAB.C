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

#include "TAB.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::TAB<CloudType>::TAB
(
    const dictionary& dict,
    CloudType& owner
)
:
    BreakupModel<CloudType>(owner),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    Cmu_(readScalar(coeffsDict_.lookup("Cmu"))),
    Comega_(readScalar(coeffsDict_.lookup("Comega"))),
    WeCrit_(readScalar(coeffsDict_.lookup("WeCrit")))
{

    // calculate the inverse function of the Rossin-Rammler Distribution
    const scalar xx0 = 12.0;
    const scalar rrd100 = 1.0/(1.0-exp(-xx0)*(1.0+xx0+pow(xx0, 2.0)/2.0 + pow(xx0, 3.0)/6.0));

    for(label n=0; n<100; n++)
    {
        scalar xx = 0.12*(n+1);
        rrd_[n] = (1.0 - exp(-xx)*(1.0 + xx + pow(xx, 2.0)/2.0 + pow(xx, 3.0)/6.0))*rrd100;
    }
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::TAB<CloudType>::~TAB()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::TAB<CloudType>::update
(
    const scalar& dt,
    const vector& g,
    scalar& d,
    scalar& tc,
    scalar& ms,
    scalar& nParticle,
    scalar& KHindex,
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
    scalar& massChild
) const
{
    // Do nothing
    return false;
}


// ************************************************************************* //
