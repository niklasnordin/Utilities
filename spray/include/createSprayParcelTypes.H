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
    Foam::createSprayParcelTypes

Description

SourceFiles
    createSprayParcelTypesI.H
    createSprayParcelTypes.C
    createSprayParcelTypesIO.C

\*---------------------------------------------------------------------------*/

#ifndef createSprayParcelTypes_H
#define createSprayParcelTypes_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class createSprayParcelTypes;
Istream& operator>>(Istream&, createSprayParcelTypes&);
Ostream& operator<<(Ostream&, const createSprayParcelTypes&);


/*---------------------------------------------------------------------------*\
                          Class createSprayParcelTypes Declaration
\*---------------------------------------------------------------------------*/

class createSprayParcelTypes
:
    public baseClassName
{
    // Private data

        //- Description of data_
        dataType data_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        createSprayParcelTypes(const createSprayParcelTypes&);

        //- Disallow default bitwise assignment
        void operator=(const createSprayParcelTypes&);


public:

    // Static data members

        //- Static data staticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        createSprayParcelTypes();

        //- Construct from components
        createSprayParcelTypes(const dataType& data);

        //- Construct from Istream
        createSprayParcelTypes(Istream&);

        //- Construct as copy
        createSprayParcelTypes(const createSprayParcelTypes&);


    // Selectors

        //- Select null constructed
        static autoPtr<createSprayParcelTypes> New();


    //- Destructor
    ~createSprayParcelTypes();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const createSprayParcelTypes&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, createSprayParcelTypes&);
        friend Ostream& operator<<(Ostream&, const createSprayParcelTypes&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "createSprayParcelTypesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
