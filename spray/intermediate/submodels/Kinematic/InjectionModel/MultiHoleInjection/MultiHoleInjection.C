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

#include "MultiHoleInjection.H"
#include "DataEntry.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::MultiHoleInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return round((time1 - time0)*parcelsPerSecond_);
    }
    else
    {
        return 0;
    }
}


template<class CloudType>
Foam::scalar Foam::MultiHoleInjection<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
) const
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        return volumeFlowRate_().integrate(time0, time1);
    }
    else
    {
        return 0.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MultiHoleInjection<CloudType>::MultiHoleInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    centerPosition_(this->coeffDict().lookup("centerPosition")),
    xyAngle_(readScalar(this->coeffDict().lookup("xyAngle"))),
    zAngle_(readScalar(this->coeffDict().lookup("zAngle"))),
    nHoles_(readLabel(this->coeffDict().lookup("nHoles"))),
    umbrellaAngle_(readScalar(this->coeffDict().lookup("umbrellaAngle"))),
    nozzleTipDiameter_(readScalar(this->coeffDict().lookup("nozzleTipDiameter"))),
    angleSpacing_(this->coeffDict().lookup("angleSpacing")),
    nozzleDiameter_(readScalar(this->coeffDict().lookup("nozzleDiameter"))),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    position_(nHoles_),
    injectorCell_(-1),
    direction_(nHoles_),
    parcelsPerSecond_
    (
        readScalar(this->coeffDict().lookup("parcelsPerSecond"))
    ),
    volumeFlowRate_
    (
        DataEntry<scalar>::New
        (
            "volumeFlowRate",
            this->coeffDict()
        )
    ),
    Cd_
    (
        DataEntry<scalar>::New
        (
            "Cd",
            this->coeffDict()
        )
    ),
    thetaInner_
    (
        DataEntry<scalar>::New
        (
            "thetaInner",
            this->coeffDict()
        )
    ),
    thetaOuter_
    (
        DataEntry<scalar>::New
        (
            "thetaOuter",
            this->coeffDict()
        )
    ),
    parcelPDF_
    (
        pdfs::pdf::New
        (
            this->coeffDict().subDict("parcelPDF"),
            owner.rndGen()
        )
    ),
    tanVec1_(nHoles_),
    tanVec2_(nHoles_),
    injectorI_(-1)
{
    scalar pi180 = mathematicalConstant::pi/180.0;
    scalar alpha = xyAngle_*pi180;
    scalar phi = zAngle_*pi180;

    // define the local coordiante system
    vector xp(cos(alpha), sin(alpha), 0.0);
    vector zp(cos(alpha)*sin(phi), sin(alpha)*sin(phi), cos(phi));

    if (mag(zp-xp) < 1.0e-15)
    {
        xp = vector(0.0, 0.0, -1.0);
        xp -= (xp & zp)*zp;
        xp /= mag(xp);
    }
    vector yp = zp^xp;

    scalar angle = 0.0;
    scalar u = umbrellaAngle_*pi180/2.0;
    for(label i=0; i<nHoles_; i++)
    {
        angle += angleSpacing_[i];
        scalar v = angle*pi180;
        direction_[i] = cos(v)*sin(u)*xp + sin(v)*sin(u)*yp + cos(u)*zp;
        vector dp = direction_[i] - (direction_[i] & zp)*direction_[i];
        if (mag(dp) > SMALL)
        {
            dp /= mag(dp);
        }
        position_[i] = centerPosition_ + 0.5*nozzleTipDiameter_*dp;
    }

    for(label i=0; i<nHoles_; i++)
    {
        vector tangent(vector::zero);
        scalar magV = 0;
        while (magV < SMALL)
        {
            vector testThis = this->owner().rndGen().vector01();
            tangent = testThis - (testThis & direction_[i])*direction_[i];
            magV = mag(tangent);
        }

        tanVec1_[i] = tangent/magV;
        tanVec2_[i] = direction_[i] ^ tanVec1_[i];

    } 

    // Set total volume to inject
    this->volumeTotal_ = volumeFlowRate_().integrate(0.0, duration_);

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MultiHoleInjection<CloudType>::~MultiHoleInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::MultiHoleInjection<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::MultiHoleInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
void Foam::MultiHoleInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector& position,
    label& cellOwner
)
{
    injectorI_ = int(nHoles_*this->owner().rndGen().scalar01());
 
    position = position_[injectorI_];
    this->findCellAtPosition(injectorCell_, position);
    cellOwner = injectorCell_;
}


template<class CloudType>
void Foam::MultiHoleInjection<CloudType>::setProperties
(
    const label parcelI,
    const label,
    const scalar time,
    const scalar pressure,
    typename CloudType::parcelType& parcel
)
{
    // set particle velocity
    const scalar deg2Rad = mathematicalConstant::pi/180.0;

    scalar t = time - this->SOI_;
    scalar ti = thetaInner_().value(t);
    scalar to = thetaOuter_().value(t);
    scalar coneAngle = this->owner().rndGen().scalar01()*(to - ti) + ti;

    coneAngle *= deg2Rad;
    scalar alpha = sin(coneAngle);
    scalar dcorr = cos(coneAngle);
    scalar beta = mathematicalConstant::twoPi*this->owner().rndGen().scalar01();

    vector normal = alpha*(tanVec1_[injectorI_]*cos(beta) + tanVec2_[injectorI_]*sin(beta));
    vector dirVec = dcorr*direction_[injectorI_];
    dirVec += normal;
    dirVec /= mag(dirVec);

    scalar A = nHoles_*0.25*mathematicalConstant::pi*nozzleDiameter_*nozzleDiameter_;
    scalar massFlowRate = this->massTotal()*volumeFlowRate_().value(t)/this->volumeTotal();

    scalar Umag = massFlowRate/(parcel.rho()*Cd_().value(t)*A);
    parcel.U() = Umag*dirVec;

    // set particle diameter
    parcel.d() = parcelPDF_().sample();
}


template<class CloudType>
bool Foam::MultiHoleInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::MultiHoleInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
