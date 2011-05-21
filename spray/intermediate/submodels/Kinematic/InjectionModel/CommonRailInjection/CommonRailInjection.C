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

#include "CommonRailInjection.H"
#include "DataEntry.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::CommonRailInjection<CloudType>::parcelsToInject
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
Foam::scalar Foam::CommonRailInjection<CloudType>::volumeToInject
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
Foam::CommonRailInjection<CloudType>::CommonRailInjection
(
    const dictionary& dict,
    CloudType& owner
)
:
    InjectionModel<CloudType>(dict, owner, typeName),
    outerNozzleDiameter_(readScalar(this->coeffDict().lookup("outerNozzleDiameter"))),
    innerNozzleDiameter_(readScalar(this->coeffDict().lookup("innerNozzleDiameter"))),
    duration_(readScalar(this->coeffDict().lookup("duration"))),
    position_(this->coeffDict().lookup("position")),
    injectorCell_(-1),
    direction_(this->coeffDict().lookup("direction")),
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
    Pinj_
    (
        DataEntry<scalar>::New
        (
            "Pinj",
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
    tanVec1_(vector::zero),
    tanVec2_(vector::zero),
    normal_(vector::zero)
{

    if (innerNozzleDiameter_ >= outerNozzleDiameter_)
    {
        FatalErrorIn
        (
            "Foam::CommonRailInjection<CloudType>::CommonRailInjection"
            "("
                "..."
            ")"
        )<< "innerNozzleDiameter >= outerNozzleDiameter" << nl
         << exit(FatalError); 
    }

    word injectionMethodType = this->coeffDict().lookup("injectionMethod");

    if (injectionMethodType == "disc")
    {
        this->injectionMethod_ = InjectionModel<CloudType>::imDisc;
    }
    else if (injectionMethodType == "point")
    {
        this->injectionMethod_ = InjectionModel<CloudType>::imPoint;
    }
    else
    {
        FatalErrorIn
        (
            "Foam::InjectionModel<CloudType>::InjectionModel"
            "("
                "const dictionary&, "
                "CloudType&, "
                "const word&"
            ")"
        )<< "injectionMethodType must be either 'point' or 'disc'" << nl
         << exit(FatalError);
    }

    // Normalise direction vector
    direction_ /= mag(direction_);

    // Determine direction vectors tangential to direction
    vector tangent = vector::zero;
    scalar magTangent = 0.0;

    while (magTangent < SMALL)
    {
        vector v = this->owner().rndGen().vector01();

        tangent = v - (v & direction_)*direction_;
        magTangent = mag(tangent);
    }

    tanVec1_ = tangent/magTangent;
    tanVec2_ = direction_^tanVec1_;

    // Set total volume to inject
    this->volumeTotal_ = volumeFlowRate_().integrate(0.0, duration_);

    // Set/cache the injector cell
    //this->findCellAtPosition(injectorCell_, position_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CommonRailInjection<CloudType>::~CommonRailInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::CommonRailInjection<CloudType>::active() const
{
    return true;
}


template<class CloudType>
Foam::scalar Foam::CommonRailInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}


template<class CloudType>
void Foam::CommonRailInjection<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector& position,
    label& cellOwner
)
{
    scalar beta = mathematicalConstant::twoPi*this->owner().rndGen().scalar01();
    normal_ = tanVec1_*cos(beta) + tanVec2_*sin(beta);

    switch (InjectionModel<CloudType>::injectionMethod_)
    {
        case InjectionModel<CloudType>::imPoint:
        {
            position = position_;
            break;
        }
        case InjectionModel<CloudType>::imDisc:
        {
            scalar frac = this->owner().rndGen().scalar01();
            scalar r = 0.5*(innerNozzleDiameter_ + frac*(outerNozzleDiameter_ - innerNozzleDiameter_));
            position = position_ + r*normal_;
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "void Foam::CommonRailInjection<CloudType>::setPositionAndCell"
                "("
                    "..."
                ")"
            )<< "Unknown injectionMethod type" << nl
             << exit(FatalError);
        }
    }
    this->findCellAtPosition(injectorCell_, position);
    cellOwner = injectorCell_;
}


template<class CloudType>
void Foam::CommonRailInjection<CloudType>::setProperties
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
    // moved to setPositionAndCell
    //scalar beta = mathematicalConstant::twoPi*this->owner().rndGen().scalar01();
    //vector normal = alpha*(tanVec1_*cos(beta) + tanVec2_*sin(beta));
    vector normal = alpha*normal_;
    vector dirVec = dcorr*direction_;
    dirVec += normal;
    dirVec /= mag(dirVec);

    scalar Umag = ::sqrt(2.0*(Pinj_().value(t) - pressure)/parcel.rho());

    parcel.U() = Umag*dirVec;

    // set particle diameter
    parcel.d() = parcelPDF_().sample();
}


template<class CloudType>
bool Foam::CommonRailInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::CommonRailInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
