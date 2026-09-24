/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "basicFastChemistryModel.H"
#include "NamedEnum.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    const NamedEnum<basicFastChemistryModel::jacobianType, 2>
    basicFastChemistryModel::jacobianTypeNames_
    {
        "fast",
        "exact"
    };
}


namespace Foam
{
    defineTypeNameAndDebug(basicFastChemistryModel, 0);
    defineRunTimeSelectionTable(basicFastChemistryModel, mesh);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::basicFastChemistryModel::correct()
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicFastChemistryModel::basicFastChemistryModel
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "chemistryProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    // thermo_(thermo),
    chemistry_(lookup("chemistry")),
    deltaTChemIni_(lookup<scalar>("initialChemicalTimeStep")),
    deltaTChemMax_(lookupOrDefault("maxChemicalTimeStep", great)),
    deltaTChem_
    (
        IOobject
        (
            "deltaTChem",
            mesh().time().constant(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimTime, deltaTChemIni_)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicFastChemistryModel::~basicFastChemistryModel()
{}


// ************************************************************************* //
