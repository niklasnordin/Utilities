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

Class
    Foam::MultiHoleInjection

Description
    Cone injection

    - User specifies
      - time of start of injection
      - injector position
      - direction (along injection axis)
      - parcel flow rate
      - Discharge coefficiend, Cd
      - inner and outer cone angles
    - Parcel diameters obtained by PDF model

    Parcel velocity is calculated as:
        U = V_dot/(A * Cd), where V_dot is the volume flow rate
SourceFiles
    MultiHoleInjection.C

\*---------------------------------------------------------------------------*/

#ifndef MultiHoleInjection_H
#define MultiHoleInjection_H

#include "InjectionModel.H"
#include "pdf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class Type>
class DataEntry;

/*---------------------------------------------------------------------------*\
                       Class MultiHoleInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MultiHoleInjection
:
    public InjectionModel<CloudType>
{
    // Private data

        //- position from which all the holes are centered
        const vector centerPosition_;

        const scalar xyAngle_;

        const scalar zAngle_;
        const label nHoles_;
        const scalar umbrellaAngle_;
        const scalar nozzleTipDiameter_;

        //- angle between the different holes
        const List<scalar> angleSpacing_;

        //- Nozzle diameter [m]
        const scalar nozzleDiameter_;

        //- Injection duration [s]
        const scalar duration_;

        //- Injector position [m]
        List<vector> position_;

        //- Cell containing injector position []
        label injectorCell_;

        //- Injector direction []
        List<vector> direction_;

        //- Number of parcels to introduce per second []
        const label parcelsPerSecond_;

        //- Volume flow rate of parcels to introduce relative to SOI [m^3/s]
        const autoPtr<DataEntry<scalar> > volumeFlowRate_;

        //- Discharge coefficient, relative to SOI [m/s]
        const autoPtr<DataEntry<scalar> > Cd_;

        //- Inner cone angle relative to SOI [deg]
        const autoPtr<DataEntry<scalar> > thetaInner_;

        //- Outer cone angle relative to SOI [deg]
        const autoPtr<DataEntry<scalar> > thetaOuter_;

        //- Parcel size PDF model
        const autoPtr<pdfs::pdf> parcelPDF_;


        // Tangential vectors to the direction vector

            //- First tangential vector
            List<vector> tanVec1_;

            //- Second tangential vector
            List<vector> tanVec2_;

            label injectorI_;

protected:

    // Protected member functions

        //- Number of parcels to introduce over the time step relative to SOI
        label parcelsToInject
        (
            const scalar time0,
            const scalar time1
        ) const;


public:

    //- Runtime type information
    TypeName("MultiHoleInjection");


    // Constructors

        //- Construct from dictionary
        MultiHoleInjection
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~MultiHoleInjection();


    // Member Functions

        //- Flag to indicate whether model activates injection model
        bool active() const;

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Volume of parcels to introduce over the time step relative to SOI
        scalar volumeToInject
        (
            const scalar time0,
            const scalar time1
        ) const;


        // Injection geometry

            //- Set the injection position and owner cell
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
		const scalar pressure,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "MultiHoleInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
