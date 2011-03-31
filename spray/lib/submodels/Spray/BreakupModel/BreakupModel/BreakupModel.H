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
    Foam::BreakupModel

Description

SourceFiles
    BreakupModelI.H
    BreakupModel.C
    BreakupModelIO.C

\*---------------------------------------------------------------------------*/

#ifndef BreakupModel_H
#define BreakupModel_H

#include ".H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class BreakupModel;
Istream& operator>>(Istream&, BreakupModel&);
Ostream& operator<<(Ostream&, const BreakupModel&);


/*---------------------------------------------------------------------------*\
                          Class BreakupModel Declaration
\*---------------------------------------------------------------------------*/

class BreakupModel
:
    public baseClassName
{
    // Private data

        //- Description of data_
        dataType data_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        BreakupModel(const BreakupModel&);

        //- Disallow default bitwise assignment
        void operator=(const BreakupModel&);


public:

    // Static data members

        //- Static data staticData
        static const dataType staticData;


    // Constructors

        //- Construct null
        BreakupModel();

        //- Construct from components
        BreakupModel(const dataType& data);

        //- Construct from Istream
        BreakupModel(Istream&);

        //- Construct as copy
        BreakupModel(const BreakupModel&);


    // Selectors

        //- Select null constructed
        static autoPtr<BreakupModel> New();


    //- Destructor
    ~BreakupModel();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Member Operators

        void operator=(const BreakupModel&);


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, BreakupModel&);
        friend Ostream& operator<<(Ostream&, const BreakupModel&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "BreakupModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
