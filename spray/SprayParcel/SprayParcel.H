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
    Foam::SprayParcel

Description
    Derived from reacingParcel, with added functionality 
    for atomization and breakup

\*---------------------------------------------------------------------------*/

#ifndef SprayParcel_H
#define SprayParcel_H

#include "ReactingParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class SprayParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class SprayParcel
:
    public ReactingParcel<SprayParcel<ThermoType> >
{
    // Private data

        //- Description of data_
        dataType data_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        SprayParcel(const SprayParcel&);

        //- Disallow default bitwise assignment
        void operator=(const SprayParcel&);


public:

    //- The type of thermodynamics this parcel was instantiated for
    typedef ThermoType thermoType;

    //- Run-time type information
    TypeName("SprayParcel");

    // Constructors

        //- Construct null
        SprayParcel();

        //- Construct from components
        SprayParcel(const dataType& data);

        //- Construct from Istream
        SprayParcel(Istream&);

        //- Construct as copy
        SprayParcel(const SprayParcel&);


    // Selectors

        //- Select null constructed
        static autoPtr<SprayParcel> New();


    //- Destructor
    ~SprayParcel();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const SprayParcel&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, SprayParcel&);
        friend Ostream& operator<<(Ostream&, const SprayParcel&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SprayParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
