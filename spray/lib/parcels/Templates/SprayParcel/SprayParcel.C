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


// * * * * * * * * * * *  Member Functions * * * * * * * * * * * * //

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

    // store the parcel properties
    const scalarField& Y(this->Y());
    scalarField X(td.cloud().composition().liquids().X(Y));

    scalar T0 = this->T();
    this->cp() = td.cloud().composition().liquids().cp(this->pc_, T0, X);
    scalar rho0 = td.cloud().composition().liquids().rho(this->pc_, T0, X);
    this->rho() = rho0;

    ReactingParcel<ParcelType>::calc(td, dt, cellI);

    // update drop cp, diameter and density because of change in temperature/composition
    scalar T1 = this->T();
    const scalarField& Y1(this->Y());
    scalarField X1(td.cloud().composition().liquids().X(Y1));

    this->cp() = td.cloud().composition().liquids().cp(this->pc_, T1, X1);

    scalar rho1 = td.cloud().composition().liquids().rho(this->pc_, T1, X1);
    this->rho() = rho1;
    scalar d1 = this->d()*pow(rho0/rho1, 1.0/3.0);
    this->d() = d1;

    if (liquidCore() > 0.5)
    {
        calcAtomization(td, dt, cellI);

        // preserve the total mass/volume, by increasing the number of particles in parcels due to breakup
        scalar d2 = this->d();
        this->nParticle() *= pow(d1/d2, 3.0);
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

    // cell state info is updated in ReactingParcel calc
    const scalarField& Y(this->Y());
    scalarField X(td.cloud().composition().liquids().X(Y));

    scalar rho = td.cloud().composition().liquids().rho(this->pc_, this->T(), X);
    scalar mu = td.cloud().composition().liquids().mu(this->pc_, this->T(), X);
    scalar sigma = td.cloud().composition().liquids().sigma(this->pc_, this->T(), X);

    // Average molecular weight of carrier mix - assumes perfect gas
    scalar Wc = this->rhoc_*specie::RR*this->Tc_/this->pc_;
    scalar R = specie::RR/Wc;
    scalar Tav = td.cloud().atomization().Taverage(this->T(), this->Tc_);

    // calculate average gas density based on average temperature
    scalar rhoAv = this->pc_/(R*Tav);

    scalar soi = td.cloud().injection().timeStart();
    scalar t1 = this->cloud().db().time().value() - soi;
    scalar t0 = t1 - dt;

    scalar massflowRate = rho*td.cloud().injection().volumeToInject(t0, t1)/dt;
    const vector& pos = this->position();
    const vector& injectionPos = this->position0();

    // disregard the continous phase when calculating the relative velocity
    scalar Urel = mag(this->U());

    scalar chi = 0.0;
    if (td.cloud().atomization().calcChi())
    {
	chi = this->chi(td, X);
    }

    td.cloud().atomization().update
    (
        dt,
        this->d(),
        this->liquidCore(),
	this->tc(),
	rho,
	mu,
	sigma,
	massflowRate,
	rhoAv,
	Urel,
	pos,
	injectionPos,
        td.cloud().pAmbient(),
	chi,
	td.cloud().rndGen()
    );

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

    // cell state info is updated in ReactingParcel calc
    const scalarField& Y(this->Y());
    scalarField X(td.cloud().composition().liquids().X(Y));

    scalar rho = td.cloud().composition().liquids().rho(this->pc_, this->T(), X);
    scalar mu = td.cloud().composition().liquids().mu(this->pc_, this->T(), X);
    scalar sigma = td.cloud().composition().liquids().sigma(this->pc_, this->T(), X);

    // Average molecular weight of carrier mix - assumes perfect gas
    scalar Wc = this->rhoc_*specie::RR*this->Tc_/this->pc_;
    scalar R = specie::RR/Wc;
    scalar Tav = td.cloud().atomization().Taverage(this->T(), this->Tc_);

    // calculate average gas density based on average temperature
    scalar rhoAv = this->pc_/(R*Tav);
    scalar muAv = this->muc_;
    vector Urel = this->U() - this->Uc_;
    scalar Urmag = mag(Urel);
    scalar As = this->areaS(this->d());
    scalar Re = rhoAv*Urmag*this->d()/muAv;

    scalar utc = td.cloud().drag().utc(Re, this->d(), muAv) + ROOTVSMALL;
    scalar tMom = 1.0/(As*utc);

    scalar averageParcelMass = 1.0;
    const vector g = td.cloud().g().value();

    scalar massChild = 0.0;
    scalar dChild = 0.0;
    if
    (
        td.cloud().breakup().update
        (
            dt,
	    g,
            this->d(),
	    this->tc(),
	    this->ms(),
	    this->nParticle(),
	    this->KHindex(),
	    this->d0(),
            rho,
            mu,
            sigma,
	    this->U(),
            rhoAv,
            muAv,
	    Urel,
            Urmag,
	    tMom,
	    averageParcelMass,
	    dChild,
	    massChild
        )
    )
    {
        // add child parcel. most properties will be identical to the parent
        ParcelType* child = new ParcelType(td.cloud(), this->position(), cellI);
	child->d() = dChild;
	child->tc() = -GREAT;
	td.cloud().addParticle(child);
    }
}


template<class ParcelType>
template<class TrackData>
Foam::scalar Foam::SprayParcel<ParcelType>::chi
(
    TrackData& td,
    const scalarField& X
) const
{


//  modifications to take account of the flash boiling on primary breakUp

    scalar chi = 0.0;
    label Nf = td.cloud().composition().liquids().components().size();

    scalar Td = this->T();
    scalar pAmb = td.cloud().pAmbient();

    for(label i = 0; i < Nf ; i++)
    {
        scalar pv = td.cloud().composition().liquids().sigma(this->pc_, this->T(), X);

        if(pv >= 0.999*pAmb)
        {

//          The fuel is boiling.....
//          Calculation of the boiling temperature

            scalar tBoilingSurface = Td;

            label Niter = 200;

            for(label k=0; k< Niter ; k++)
            {
                scalar pBoil = td.cloud().composition().liquids().properties()[i].pv(this->pc_, tBoilingSurface);

                if(pBoil > this->pc_)
                {
                    tBoilingSurface = tBoilingSurface - (Td-this->Tc_)/Niter;
                }
                else
                {
                    break;
                }
            }

            scalar hl = td.cloud().composition().liquids().properties()[i].hl(pAmb, tBoilingSurface);
            scalar iTp = td.cloud().composition().liquids().properties()[i].h(pAmb, Td) - pAmb/td.cloud().composition().liquids().properties()[i].rho(pAmb, Td);
            scalar iTb = td.cloud().composition().liquids().properties()[i].h(pAmb, tBoilingSurface) - pAmb/td.cloud().composition().liquids().properties()[i].rho(pAmb, tBoilingSurface);

            chi += X[i]*(iTp-iTb)/hl;

        }
    }

    //  bound chi
    chi = max(chi, 0.0);
    chi = min(chi, 1.0);

  return chi;
}

// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "SprayParcelIO.C"


// ************************************************************************* //
