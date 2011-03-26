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

Class
    Foam::SprayCloud

Description

\*---------------------------------------------------------------------------*/

#ifndef SprayCloud_H
#define SprayCloud_H

#include "ReactingCloud.H"
#include "thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                      Class SprayCloud Declaration
\*---------------------------------------------------------------------------*/

namespace Foam
{

template<class ParcelType>
class SprayCloud
:
    public ReactingCloud<ParcelType>
{
public:

    //- Type of thermodynamics the cloud was instantiated for
    typedef typename ParcelType::thermoType thermoType;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        SprayCloud(const SprayCloud&);

        //- Disallow default bitwise assignment
        void operator=(const SprayCloud&);


protected:

        //- Parcel constant properties
        typename ParcelType::constantProperties constProps_;

        // Cloud evolution functions

            //- Pre-evolve
            void preEvolve();

            //- Evolve the cloud
            void evolveCloud();

            //- Post-evolve
            void postEvolve();


public:

    // Constructors

        //- Construct given carrier gas fields
        SprayCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const dimensionedVector& g,
            basicThermo& thermo,
            bool readFields = true
        );


    //- Destructor
    virtual ~SprayCloud();


    //- Type of parcel the cloud was instantiated for
    typedef ParcelType parcelType;

    // Member Functions

        // Access

            //- Return the constant properties
            inline const typename ParcelType::constantProperties&
                constProps() const;

        // Check

            //- Print cloud information
            void info() const;

        // Cloud evolution functions

            //- Check parcel properties
            void checkParcelProperties
            (
                ParcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Reset the spray source terms
            void resetSourceTerms();

            //- Evolve the spray (inject, move)
            void evolve();
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SprayCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SprayCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
