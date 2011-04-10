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

#include "SprayParcel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ParcelType>
Foam::SprayParcel<ParcelType>::SprayParcel
(
    const SprayParcel<ParcelType>& p
)
:
    ReactingParcel<ParcelType>(p),
    d0_(p.d0_),
    position0_(p.position0_),
    liquidCore_(p.liquidCore_),
    KHindex_(p.KHindex_),
    y_(p.y_),
    yDot_(p.yDot_),
    tc_(p.tc_),
    ms_(p.ms_),
    injector_(p.injector_),
    tMom_(p.tMom_),
    user_(p.user_)
{}


// * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * * //

// NN. Dont think all these functions are needed, but I'm adding them in case 
//     one might have to add anything later

template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::setCellValues
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ReactingParcel<ParcelType>::setCellValues(td, dt, cellI);
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::cellValueSourceCorrection
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
    ReactingParcel<ParcelType>::cellValueSourceCorrection(td, dt, cellI);
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::correctSurfaceValues
(
    TrackData& td,
    const label cellI,
    const scalar T,
    const scalarField& Cs,
    scalar& rhos,
    scalar& mus,
    scalar& Pr,
    scalar& kappa
)
{
    ReactingParcel<ParcelType>::correctSurfaceValues
    (
        td,
        cellI,
        T,
        Cs,
        rhos,
        mus,
        Pr,
        kappa
    );
}

template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::calc
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{

    // check if parcel belongs to liquid core
    if (liquidCore() > 0.5)
    {
        // liquid core parcels should not interact with the gas
        if (td.cloud().coupled())
	{
	    td.cloud().coupled() = false;
	}
    }

    ReactingParcel<ParcelType>::calc(td, dt, cellI);

    if (liquidCore() > 0.5)
    {
        calcAtomization(td, dt, cellI);
    }
    else
    {
        calcBreakup(td, dt, cellI);
    }

    // check if parcel is uncoupled and does not belong to liquid core
    if ( (!td.cloud().coupled()) && (liquidCore() < 0.5) )
    {
        td.cloud().coupled() = true;
    }

}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::calcAtomization
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
  /*
    // cell state info is updated in ReactingParcel calc

    const scalarField& Y(this->Y());
    scalarField X(td.cloud().composition().liquids().X(Y));

    scalar rho = td.cloud().composition().liquids().rho(this->pc_, this->T(), X);
    scalar mu = td.cloud().composition().liquids().mu(this->pc_, this->T(), X);
    scalar sigma = td.cloud().composition().liquids().sigma(this->pc_, this->T(), X);

    // Average molecular weight of carrier mix - assumes perfect gas
    scalar Wc = this->rhoc_*specie::RR*this->Tc_/this->pc_;

    scalar t0 = this->cloud().db().time().value();
    scalar t1 = t0 + dt;
    Info << "t0 = " << t0 << endl;
    scalar massflowRate = rho*td.cloud().injection().volumeToInject(t0, t1)/dt;
  */

    /*
    td.cloud().atomization().update
    (
        this->d(),
        this->liquidCore(),
	this->tc(),
	rho,
	mu,
	sigma,
	massflowRate
    );
    */
}


template<class ParcelType>
template<class TrackData>
void Foam::SprayParcel<ParcelType>::calcBreakup
(
    TrackData& td,
    const scalar dt,
    const label cellI
)
{
}

// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "SprayParcelIO.C"


// ************************************************************************* //
