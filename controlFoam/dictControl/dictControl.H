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
    Foam::dictControl

Description

SourceFiles
    dictControlI.H
    dictControl.C
    dictControlIO.C

\*---------------------------------------------------------------------------*/

#ifndef dictControl_H
#define dictControl_H

#include "IOdictionary.H"
#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class dictControl;
Istream& operator>>(Istream&, dictControl&);
Ostream& operator<<(Ostream&, const dictControl&);


/*---------------------------------------------------------------------------*\
                          Class dictControl Declaration
\*---------------------------------------------------------------------------*/

class dictControl
{
    // Private data

        //- Description of data_
        word name_;
        dictionary dict_;
        string options_;
        string script_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        dictControl(const dictControl&);

        //- Disallow default bitwise assignment
        void operator=(const dictControl&);


public:

    // Constructors

        //- Construct null
        dictControl
        (
            const IOdictionary& dict,
            const word& name
        );
  /*
    // Selectors

        //- Select null constructed
        static autoPtr<dictControl> New();
  */

    //- Destructor
    ~dictControl();


    // Member Functions

        // Access

        inline const word& name() const
        {
	    return name_;
        }

        inline const dictionary& dict() const
        {
	    return dict_;
        }

        inline const string& options() const
        {
	    return options_;
	}

        inline const string& script() const
        {
	    return script_;
	}

        // Check

        // Edit

        // Write

    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators

        friend Istream& operator>>(Istream&, dictControl&);
        friend Ostream& operator<<(Ostream&, const dictControl&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "dictControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
