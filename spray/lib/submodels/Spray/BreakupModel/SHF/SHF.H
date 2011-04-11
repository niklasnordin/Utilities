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
    Foam::SHF

Description
    Secondary Breakup Model to take account of the different breakup regimes,
    bag, molutimode, shear....

    Accurate description in
    @verbatim
    R. Schmehl, G. Maier, S. Witting
    "CFD Analysis of Fuel Atomization, Secondary Droplet Breakup and Spray
    Dispersion in the Premix Duct of a LPP Combustor".
    Eight International Conference on Liquid Atomization and Spray Systems, 2000
    @endverbatim

\*---------------------------------------------------------------------------*/

#ifndef SHF_H
#define SHF_H

#include "BreakupModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class SHF Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SHF
:
    public BreakupModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("SHF");


    // Constructors

        //- Construct from dictionary
        SHF(const dictionary&, CloudType&);


    //- Destructor
    virtual ~SHF();


    // Member Functions

        //- update the parcel diameter
        virtual void breakup
        (
            const scalar& dt,
	    scalar& d,
	    const scalar& rho,
	    const scalar& mu,
	    const scalar& sigma,
	    const scalar& rhoc,
	    const scalar& muc,
	    const scalar& Urel
	) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SHF.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
