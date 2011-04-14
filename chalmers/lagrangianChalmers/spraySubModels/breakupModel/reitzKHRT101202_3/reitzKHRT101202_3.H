/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::reitzKHRT101202_3

Description
    secondary breakup model which uses the Kelvin-Helmholtz
    instability theory to predict the 'stripped' droplets... and
    the Raleigh-Taylor instability as well.

\*---------------------------------------------------------------------------*/

#ifndef reitzKHRT101202_3_H
#define reitzKHRT101202_3_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class reitzKHRT101202_3 Declaration
\*---------------------------------------------------------------------------*/

class reitzKHRT101202_3
:
    public breakupModel
{

private:

    // Private data

        dictionary coeffsDict_;

        // reference to gravity
        const vector& g_;

        // model constants
        scalar b0_;
        scalar b1_;
        scalar cTau_;
        scalar cRT_;
        scalar msLimit_;
        scalar weberLimit_;

public:

    //- Runtime type information
    TypeName("ReitzKHRT101202_3");


    // Constructors

        //- Construct from components
        reitzKHRT101202_3
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        ~reitzKHRT101202_3();


    // Member Operators

        void breakupParcel
        (
            parcel& parcel,
            const scalar deltaT,
            const vector& vel,
            const liquidMixture& fuels
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
