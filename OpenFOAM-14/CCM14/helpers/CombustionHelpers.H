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

#ifndef CombustionHelpers_H
#define CombustionHelpers_H

#include "fvMesh.H"
#include "volFields.H"
#include "dictionary.H"
#include "Switch.H"
#include "scalarList.H"
#include "hashedWordList.H"
#include "fluidMulticomponentThermo.H"
#include "atomicWeights.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class CombustionHelpers Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class CombustionHelpers
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to thermo
        const fluidMulticomponentThermo& thermo_;

        //- Reference to the multi component mixture
        const multicomponentMixture<ThermoType>& mixture_;

        //- Reference to thermodynamic data of the species
        const PtrList<ThermoType>& specieThermos_;        

        //- Reference to CCM dictionary
        const dictionary& CCMdict_;

        //- number of species
        const label& nSpecie_;

        //- Reference to mass fraction fields
        const PtrList<volScalarField>& Yvf_;

        //- Reference to update variables list
        hashedWordList& updateVars_;

        //- References to J coefficient arrays
        scalarList& JhCoeff_;
        scalarList& JcCoeff_;
        scalarList& JnCoeff_;
        scalarList& JoCoeff_;
        scalarList& Jh_h2oCoeff_;
        scalarList& Jc_co2Coeff_;
        scalarList& Jo_co2_h2oCoeff_;

        //- References to element mass fraction fields
        volScalarField& Jh_;
        volScalarField& Jc_;
        volScalarField& Jn_;
        volScalarField& Jo_;
        volScalarField& Jh_h2o_;
        volScalarField& Jc_co2_;
        volScalarField& Jo_co2_h2o_;
        volScalarField& J_;

        //- References to equivalence ratio and scalar dissipation rate fields
        volScalarField& phieq_;
        volScalarField& chi_;

        //- Reference to style switches
        const Switch& oldStyleChi_;
        const Switch& oldStylePhi_;

        //- Reference to fuel oxygen to carbon ratio
        const scalar& fuelOtoC_;

        //- temporary variables
        hashedWordList& principalVars_;

public:

    // Constructors

        //- Construct from references to required data
        CombustionHelpers
        (
            const fvMesh& mesh,
            const fluidMulticomponentThermo& thermo,
            const multicomponentMixture<ThermoType>& mixture,
            const PtrList<ThermoType>& specieThermos,
            const dictionary& CCMdict,
            const label& nSpecie,
            const PtrList<volScalarField>& Yvf,
            hashedWordList& updateVars,
            scalarList& JhCoeff,
            scalarList& JcCoeff,
            scalarList& JnCoeff,
            scalarList& JoCoeff,
            scalarList& Jh_h2oCoeff,
            scalarList& Jc_co2Coeff,
            scalarList& Jo_co2_h2oCoeff,
            volScalarField& Jh,
            volScalarField& Jc,
            volScalarField& Jn,
            volScalarField& Jo,
            volScalarField& Jh_h2o,
            volScalarField& Jc_co2,
            volScalarField& Jo_co2_h2o,
            volScalarField& J,
            volScalarField& phieq,
            volScalarField& chi,
            const Switch& oldStyleChi,
            const Switch& oldStylePhi,
            const scalar& fuelOtoC,
            hashedWordList& principalVars
        );

        //- Disallow default bitwise copy construction
        CombustionHelpers(const CombustionHelpers&) = delete;


    //- Destructor
    ~CombustionHelpers() = default;


    // Member Functions

        //- Update carbon element mass fraction
        void updateJc();

        //- Update hydrogen element mass fraction
        void updateJh();

        //- Update oxygen element mass fraction
        void updateJo();

        //- Update nitrogen element mass fraction
        void updateJn();

        //- Update hydrogen element mass fraction excluding H2O
        void updateJh_h2o();

        //- Update carbon element mass fraction excluding CO2
        void updateJc_co2();

        //- Update oxygen element mass fraction excluding CO2 and H2O
        void updateJo_co2_h2o();

        //- Update equivalence ratio
        void updatePhi();

        //- Update scalar dissipation rate
        void updateChi();

        //- Update principal variables (excluding species mass fractions and temperature)
        void updatePV();

        //- Append varName to list if not existing
        bool appendIfNotExisting(const word& varName, hashedWordList& list);

        //- Output hashedWordList in good format
        void outputHashedWordList(const hashedWordList& list) const;

        //- Initialize principal variable lists
        void initPvLists();

        //- Initialize J coefficients for element mass fractions
        void initJCoeffs();

    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CombustionHelpers&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CombustionHelpers.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //