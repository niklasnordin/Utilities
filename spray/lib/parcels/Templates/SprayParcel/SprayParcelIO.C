/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "SprayParcel.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template <class ParcelType>
Foam::string Foam::SprayParcel<ParcelType>::propHeader =
    ReactingParcel<ParcelType>::propHeader
+ " niklas";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::SprayParcel<ParcelType>::SprayParcel
(
    const Cloud<ParcelType>& cloud,
    Istream& is,
    bool readFields
)
:
    ReactingParcel<ParcelType>(cloud, is, readFields)
{

    // Check state of Istream
    is.check
    (
        "SprayParcel<ParcelType>::SprayParcel"
        "("
            "const Cloud<ParcelType>&, "
            "Istream&, "
            "bool"
        ")"
    );
}


template<class ParcelType>
void Foam::SprayParcel<ParcelType>::readFields(Cloud<ParcelType>& cIn)
{
    if (!cIn.size())
    {
        return;
    }

    SprayCloud<ParcelType>& c =
        dynamic_cast<SprayCloud<ParcelType>&>(cIn);

    ReactingParcel<ParcelType>::readFields(c);

}


template<class ParcelType>
void Foam::SprayParcel<ParcelType>::writeFields
(
    const Cloud<ParcelType>& cIn
)
{
    const SprayCloud<ParcelType>& c =
        dynamic_cast<const SprayCloud<ParcelType>&>(cIn);

    ReactingParcel<ParcelType>::writeFields(c);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const SprayParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ReactingParcel<ParcelType>&>(p);
    }
    else
    {
        os  << static_cast<const ReactingParcel<ParcelType>&>(p);
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<(Ostream&, const SprayParcel<ParcelType>&)"
    );

    return os;
}


// ************************************************************************* //
