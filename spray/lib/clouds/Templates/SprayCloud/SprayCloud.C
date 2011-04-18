/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "SprayCloud.H"
#include "AtomizationModel.H"
#include "BreakupModel.H"
#include "CollisionModel.H"
#include "PtrList.H"

template<class ParcelType>
void Foam::SprayCloud<ParcelType>::preEvolve()
{
    ReactingCloud<ParcelType>::preEvolve();
    
    pAmbient_ = this->carrierThermo().p().average().value();
}


template<class ParcelType>
void Foam::SprayCloud<ParcelType>::evolveCloud()
{
    const volScalarField& T = this->carrierThermo().T();
    const volScalarField cp = this->carrierThermo().Cp();
    const volScalarField& p = this->carrierThermo().p();

    autoPtr<interpolation<scalar> > rhoInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->rho()
    );

    autoPtr<interpolation<vector> > UInterp = interpolation<vector>::New
    (
        this->interpolationSchemes(),
        this->U()
    );

    autoPtr<interpolation<scalar> > muInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        this->mu()
    );

    autoPtr<interpolation<scalar> > TInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        T
    );

    autoPtr<interpolation<scalar> > cpInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        cp
    );

    autoPtr<interpolation<scalar> > pInterp = interpolation<scalar>::New
    (
        this->interpolationSchemes(),
        p
    );

    typename ParcelType::trackData td
    (
        *this,
        constProps_,
        rhoInterp(),
        UInterp(),
        muInterp(),
        TInterp(),
        cpInterp(),
        pInterp(),
        this->g().value()
    );

    if (this->coupled())
    {
        resetSourceTerms();
    }
    
    if (collision().active())
    {

        label i = 0;
        scalar dt = this->db().time().deltaTValue();
        forAllIter(typename Cloud<ParcelType>, *this, iter)
        {
            ParcelType& p = iter();
            scalar Vi = this->mesh().V()[p.cell()];
            scalarField X1(this->composition().liquids().X(p.Y()));
            scalar sigma1 = this->composition().liquids().sigma(p.pc(), p.T(), X1);
            
            label j = 0;
            forAllIter(typename Cloud<ParcelType>, *this, jter)
            {
                if (j > i)
                {
                    ParcelType& q = jter();
                    scalar Vj = this->mesh().V()[q.cell()];
                    scalarField X2(this->composition().liquids().X(q.Y()));
                    scalar sigma2 = this->composition().liquids().sigma(q.pc(), q.T(), X2);
                    bool updateRho = collision().update
                    (
                        dt,
                        this->rndGen(),
                        p.position(),
                        p.mass0(),
                        p.d(),
                        p.nParticle(),
                        p.U(),
                        p.rho(),
                        p.T(),
                        p.Y(),
                        sigma1,
                        p.cell(),
                        Vi,
                        q.position(),
                        q.mass0(),
                        q.d(),
                        q.nParticle(),
                        q.U(),
                        q.rho(),
                        q.T(),
                        q.Y(),
                        sigma2,
                        q.cell(),
                        Vj
                    );
                    
                    // for coalescence we need to update the density and slightly correct
                    // the statistical number of particles
                    if (updateRho)
                    {
                        if (p.mass0() > VSMALL)
                        {
                            scalarField Xp(this->composition().liquids().X(p.Y()));
                            p.rho() = this->composition().liquids().rho(p.pc(), p.T(), Xp);
                            scalar md = p.rho()*mathematicalConstant::pi*pow(p.d(), 3.0)/6.0;
                            p.nParticle() = p.mass0()/md;
                        }

                        if (q.mass0() > VSMALL)
                        {
                            scalarField Xq(this->composition().liquids().X(q.Y()));
                            q.rho() = this->composition().liquids().rho(q.pc(), q.T(), Xq);
                            scalar md = q.rho()*mathematicalConstant::pi*pow(q.d(), 3.0)/6.0;
                            q.nParticle() = q.mass0()/md;
                        }

                    }
                }
                j++;
            }
            
            i++;
        }
        // remove coalesced particles (diameter set to 0)
        forAllIter(typename Cloud<ParcelType>, *this, iter)
        {
            ParcelType& p = iter();
            if (p.mass0() < VSMALL)
            {
                deleteParticle(p);
            }
        }
    }

    Cloud<ParcelType>::move(td);
    this->injection().inject(td);
}


template<class ParcelType>
void Foam::SprayCloud<ParcelType>::postEvolve()
{
    ReactingCloud<ParcelType>::postEvolve();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::SprayCloud<ParcelType>::SprayCloud
(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    basicThermo& thermo,
    bool readFields
)
:
    ReactingCloud<ParcelType>(cloudName, rho, U, g, thermo, false),
    sprayCloud(),
    pAmbient_(0.0),
    averageParcelMass_(this->injection().averageParcelMass()),
    constProps_(this->particleProperties()),
    atomizationModel_
    (
        AtomizationModel<SprayCloud<ParcelType> >::New
        (
	     this->particleProperties(),
	     *this
        )
    ),
    breakupModel_
    (
        BreakupModel<SprayCloud<ParcelType> >::New
        (
	     this->particleProperties(),
	     *this
        )
    ),
    collisionModel_
    (
        CollisionModel<SprayCloud<ParcelType> >::New
        (
	     this->particleProperties(),
	     *this
        )
    )
{
    if (readFields)
    {
        ParcelType::readFields(*this);
    }
    Info << "Average parcel mass: " << averageParcelMass_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::SprayCloud<ParcelType>::~SprayCloud()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
void Foam::SprayCloud<ParcelType>::checkParcelProperties
(
    ParcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed
)
{
    ReactingCloud<ParcelType>::checkParcelProperties
    (
        parcel,
        lagrangianDt,
        fullyDescribed
    );

    const scalarField& Y(parcel.Y());
    scalarField X(this->composition().liquids().X(Y));

    // override rho and cp from constantProperties
    parcel.cp() = this->composition().liquids().cp(parcel.pc(), parcel.T(), X);
    parcel.rho() = this->composition().liquids().rho(parcel.pc(), parcel.T(), X);

    // store the injection position and initial drop size
    parcel.position0() = parcel.position();
    parcel.d0() = parcel.d();

    parcel.y() = breakup().y0();
    parcel.yDot() = breakup().yDot0();
}


template<class ParcelType>
void Foam::SprayCloud<ParcelType>::resetSourceTerms()
{
    ReactingCloud<ParcelType>::resetSourceTerms();
}


template<class ParcelType>
void Foam::SprayCloud<ParcelType>::evolve()
{
    if (this->active())
    {
        preEvolve();

        evolveCloud();

        postEvolve();

        info();
        Info<< endl;
    }
}


template<class ParcelType>
void Foam::SprayCloud<ParcelType>::info() const
{
    ReactingCloud<ParcelType>::info();
}


// ************************************************************************* //
