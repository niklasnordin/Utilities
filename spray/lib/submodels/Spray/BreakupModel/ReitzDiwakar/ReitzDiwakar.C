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

\*---------------------------------------------------------------------------*/

#include "ReitzDiwakar.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ReitzDiwakar<CloudType>::ReitzDiwakar
(
    const dictionary& dict,
    CloudType& owner
)
:
    BreakupModel<CloudType>(owner),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    Cbag_(readScalar(coeffsDict_.lookup("Cbag"))),
    Cb_(readScalar(coeffsDict_.lookup("Cb"))),
    Cstrip_(readScalar(coeffsDict_.lookup("Cstrip"))),
    Cs_(readScalar(coeffsDict_.lookup("Cs")))
{}

    /*
        These are the default values for this model...
        static const scalar Cbag   = 6.0;
        static const scalar Cb     = 0.785;
        static const scalar Cstrip = 0.5;
        static const scalar Cs     = 10.0;
    */


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template <class CloudType>
Foam::ReitzDiwakar<CloudType>::~ReitzDiwakar()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ReitzDiwakar<CloudType>::breakup() const
{
  // Do nothing
}


// ************************************************************************* //
