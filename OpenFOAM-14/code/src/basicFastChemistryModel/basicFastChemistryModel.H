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

Class
    Foam::basicFastChemistryModel

Description
    Base class for chemistry models

SourceFiles
    basicFastChemistryModelI.H
    basicFastChemistryModel.C

\*---------------------------------------------------------------------------*/

#ifndef basicFastChemistryModel_H
#define basicFastChemistryModel_H

// #include "fluidReactionThermo.H"
#include "volFields.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

/*---------------------------------------------------------------------------*\
                     class basicFastChemistryModel Declaration
\*---------------------------------------------------------------------------*/

class basicFastChemistryModel
:
    public IOdictionary
{
public:

    // Public Enumerations

        //- Enumeration for the type of Jacobian to be calculated by the
        //  chemistry model
        enum class jacobianType
        {
            fast,
            exact
        };

        //- Jacobian type names
        static const NamedEnum<jacobianType, 2> jacobianTypeNames_;


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Chemistry activation switch
        Switch chemistry_;

        //- Initial chemical time step
        const scalar deltaTChemIni_;

        //- Maximum chemical time step
        const scalar deltaTChemMax_;

        //- Latest estimation of integration step
        volScalarField::Internal deltaTChem_;


    // Protected Member Functions

        //- Correct function - updates due to mesh changes
        void correct();


public:

    //- Runtime type information
    TypeName("basicFastChemistryModel");


    //- Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicFastChemistryModel,
        mesh,
        (const fvMesh& mesh),
        (mesh)
    );


    // Constructors

        //- Construct from mesh
        basicFastChemistryModel(const fvMesh& mesh);

        //- Disallow default bitwise copy construction
        basicFastChemistryModel(const basicFastChemistryModel&) = delete;


    // Selectors

        //- Select based on mesh
        static autoPtr<basicFastChemistryModel> New
        (
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~basicFastChemistryModel();


    // Member Functions

        //- Return const access to the mesh database
        inline const fvMesh& mesh() const;

        //- Return const access to the thermo
        // inline const fluidReactionThermo& thermo() const;

        //- Chemistry activation switch
        inline Switch chemistry() const;

        //- The number of species
        virtual label nSpecie() const = 0;

        //- The number of reactions
        virtual label nReaction() const = 0;

        //- Return the latest estimation of integration step
        inline const volScalarField::Internal& deltaTChem() const;

        // Rate evaluation on the basis of cells

            //- Get reaction rates given Y, T, P
            //  \param Y Species mass fractions
            //  \param T Temperature [K]
            //  \param p Pressure [Pa]
            //  \param deltaT Time step [s]
            //  \param deltaTChem Chemical time step [s] (output)
            //  \param rho Density [kg/m^3]
            //  \param rho0 Old time density [kg/m^3] (not sure)
            //  \return Reaction rates for all species [kg/m^3/s] (not sure)
        virtual scalarField getRRGivenYTP
        (
            const scalarField& Y,
            const scalar T,
            const scalar p,
            const scalar deltaT,
            scalar& deltaTChem,
            const scalar& rho,
            const scalar& rho0
        ) const = 0;


        // Functions to be derived in derived classes

            // Fields

                //- Return const access to chemical source terms [kg/m^3/s]
                virtual const volScalarField::Internal& RR
                (
                    const label i
                ) const = 0;

                //- Return access to chemical source terms [kg/m^3/s]
                virtual volScalarField::Internal& RR
                (
                    const label i
                ) = 0;

                //- Return reaction rate of the speciei in reactioni
                virtual tmp<volScalarField::Internal> calculateRR
                (
                    const label reactioni,
                    const label speciei
                ) const = 0;


            // Chemistry solution

                //- Calculates the reaction rates
                virtual void calculate() = 0;

                //- Solve the reaction system for the given time step
                //  and return the characteristic time
                virtual scalar solve(const scalar deltaT) = 0;

                //- Solve the reaction system for the given time step
                //  and return the characteristic time
                virtual scalar solve(const scalarField& deltaT) = 0;

                //- Return the chemical time scale
                virtual tmp<volScalarField> tc() const = 0;

                //- Return the heat release rate [kg/m/s^3]
                virtual tmp<volScalarField> Qdot() const = 0;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const basicFastChemistryModel&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "basicFastChemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
