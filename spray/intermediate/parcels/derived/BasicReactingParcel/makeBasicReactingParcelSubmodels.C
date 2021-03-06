/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "BasicReactingParcel.H"

// Kinematic
#include "makeReactingParcelDispersionModels.H"
#include "makeReactingParcelDragModels.H"
#include "makeReactingParcelInjectionModels.H"
#include "makeReactingParcelPatchInteractionModels.H"
#include "makeReactingParcelPostProcessingModels.H"

// Thermodynamic
#include "makeReactingParcelHeatTransferModels.H"

// Reacting
#include "makeReactingParcelCompositionModels.H"
#include "makeReactingParcelPhaseChangeModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Kinematic sub-models
    makeReactingDispersionModels(BasicReactingParcel);
    makeReactingDragModels(BasicReactingParcel);
    makeReactingInjectionModels(BasicReactingParcel);
    makeReactingPatchInteractionModels(BasicReactingParcel);
    makeReactingPostProcessingModels(BasicReactingParcel);

    // Thermo sub-models
    makeReactingHeatTransferModels(BasicReactingParcel);

    // Reacting sub-models
    makeReactingCompositionModels(BasicReactingParcel);
    makeReactingPhaseChangeModels(BasicReactingParcel);
};


// ************************************************************************* //
