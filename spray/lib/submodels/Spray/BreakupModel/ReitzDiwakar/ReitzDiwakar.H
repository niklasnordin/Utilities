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
    Foam::ReitzDiwakar

Description
    secondary breakup model

\*---------------------------------------------------------------------------*/

#ifndef ReitzDiwakar_H
#define ReitzDiwakar_H

#include "BreakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class ReitzDiwakar Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ReitzDiwakar
:
    public BreakupModel<CloudType>
{
private:

    // Private data

        dictionary coeffsDict_;
        scalar Cbag_;
        scalar Cb_;
        scalar Cstrip_;
        scalar Cs_;

public:

    //- Runtime type information
    TypeName("ReitzDiwakar");


    // Constructors

        //- Construct from dictionary
        ReitzDiwakar(const dictionary&, CloudType&);


    //- Destructor
    virtual ~ReitzDiwakar();


    // Member Functions

        //- update the parcel properties
        virtual bool update
        (
	    const scalar& dt,
	    scalar& d,
	    scalar& tc,
	    scalar& ms,
	    scalar& nParticle,
	    const scalar& rho,
	    const scalar& mu,
	    const scalar& sigma,
	    const vector& U,
	    const scalar& rhoc,
	    const scalar& muc,
	    const vector& Urel,
	    const scalar& Urmag,
	    const scalar& tMom,
	    const scalar& averageParcelMass,
	    scalar& dChild,
	    scalar& massChild
	) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ReitzDiwakar.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
