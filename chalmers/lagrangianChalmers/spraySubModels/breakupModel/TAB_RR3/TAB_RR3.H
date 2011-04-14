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
    ANY WARRANTY; without even the implied warranty of MERCHANTAB_RR3ILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::TAB_RR3

Description
    The TAB_RR3 Method for Numerical Calculation of Spray Droplet Breakup.

    @verbatim
        O'Rourke, P.J. and Amsden, A.A.,
        "The TAB_RR3 Method for Numerical Calculation of Spray Droplet Breakup,"
        1987 SAE International Fuels and Lubricants Meeting and Exposition,
        Toronto, Ontario, November 2-5, 1987,
        Los Alamos National Laboratory document LA-UR-87-2105;
        SAE Technical Paper Series, Paper 872089.
    @endverbatim

    This implementation follows the kiva version.

See Also
    The Enhanced %TAB_RR3 model - ETAB_RR3

SourceFiles
    TAB_RR3I.H
    TAB_RR3.C
    TAB_RR3IO.C

\*---------------------------------------------------------------------------*/

#ifndef TAB_RR3_H
#define TAB_RR3_H

#include "breakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class TAB_RR3 Declaration
\*---------------------------------------------------------------------------*/

class TAB_RR3
:
    public breakupModel
{
private:

    // Private data

        dictionary coeffsDict_;

    // inverse function approximation
    // of the Rossin-Rammler Distribution
    // used when calculating the droplet size after breakup
//        scalar rrd_[100];

        // model constants

        scalar Cmu_;
        scalar Comega_;
        scalar WeCrit_;
        scalar K_;

public:


    //- Runtime type information
    TypeName("TAB_RR3");


    // Constructors

        //- Construct from components
        TAB_RR3
        (
            const dictionary& dict,
            spray& sm
        );


    // Destructor

        ~TAB_RR3();


    // Member Functions

        void breakupParcel
        (
            parcel& parcel,
            const scalar deltaT,
            const vector& Ug,
            const liquidMixture& fuels
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
