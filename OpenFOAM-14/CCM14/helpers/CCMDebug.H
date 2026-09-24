/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2025 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2021 Shijie Xu, Shenghui Zhong
                             | Copyright (C) 2025 Yuchen Zhou
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM CCM extension.

    OpenFOAM CCM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM CCM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM CCM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Chemistry Coordinate Mapping (CCM) acceleration method for OpenFOAM-10

\*---------------------------------------------------------------------------*/

#ifndef CCMDebug_H
#define CCMDebug_H

#include "scalarField.H"
#include "wordList.H"
#include "label.H"
#include "multicomponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class CCMDebug Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class CCMDebug
{
    // Private data

        //- Reference to number of species
        const label& nSpecie_;

        //- Reference to global maximum mass fraction field
        const scalarField& Ymax_;

        //- Reference to mixture for species information
        const multicomponentMixture<ThermoType>& mixture_;

        //- References to timing arrays
        scalarField& stepTimes_;
        wordList& stepNames_;
        label& currentStepIndex_;

public:

    // Constructors

        //- Construct from references to required data
        CCMDebug
        (
            const label& nSpecie,
            const scalarField& Ymax,
            const multicomponentMixture<ThermoType>& mixture,
            scalarField& stepTimes,
            wordList& stepNames,
            label& currentStepIndex
        );

        //- Disallow default bitwise copy construction
        CCMDebug(const CCMDebug&) = delete;


    //- Destructor
    ~CCMDebug() = default;


    // Member Functions

        //- Output detailed timing analysis
        void outputTimingAnalysis(const scalarField& stepTimes) const;

        //- Get next step index and increment counter
        label getNextStepIndex() const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CCMDebug&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CCMDebug.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //