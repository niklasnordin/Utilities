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
    Templated phase change model class

SourceFiles
    BreakupModel.C
    NewBreakupModel.C

\*---------------------------------------------------------------------------*/

#ifndef BreakupModel_H
#define BreakupModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class BreakupModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class BreakupModel
{

protected:

    // Protected data

        //- The cloud dictionary
        const dictionary& dict_;

        //- Reference to the owner cloud class
        CloudType& owner_;

        //- The coefficient dictionary
        const dictionary coeffDict_;


public:

    //- Runtime type information
    TypeName("BreakupModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        BreakupModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        BreakupModel(CloudType& owner);

        //- Construct from dictionary
        BreakupModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );


    //- Destructor
    virtual ~BreakupModel();


    //- Selector
    static autoPtr<BreakupModel<CloudType> > New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Access

        //- Return the owner cloud object
        const CloudType& owner() const;

        //- Return the cloud dictionary
        const dictionary& dict() const;

        //- Return the coefficient dictionary
        const dictionary& coeffDict() const;


    // Member Functions

        //- Flag to indicate whether model activates phase change model
        virtual void breakup() const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeBreakupModel(CloudType)                                       \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(BreakupModel<CloudType>, 0);      \
                                                                              \
    defineTemplateRunTimeSelectionTable                                       \
    (                                                                         \
        BreakupModel<CloudType>,                                          \
        dictionary                                                            \
    );


#define makeBreakupModelThermoType(SS, CloudType, ParcelType, ThermoType) \
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        SS<CloudType<ParcelType<ThermoType> > >,                              \
        0                                                                     \
    );                                                                        \
                                                                              \
    BreakupModel<CloudType<ParcelType<ThermoType> > >::                   \
        adddictionaryConstructorToTable                                       \
            <SS<CloudType<ParcelType<ThermoType> > > >                        \
            add##SS##CloudType##ParcelType##ThermoType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "BreakupModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
