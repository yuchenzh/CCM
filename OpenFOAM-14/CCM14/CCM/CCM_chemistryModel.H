/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2026 OpenFOAM Foundation
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
    Foam::chemistryModels::CCM

Description
    Extension to Foam::chemistryModels::CCM templated on thermo
    and provides stiff ODE integration functions.

    Integrates a CCM OpenFOAM reaction system and evaluation of chemical
    source terms with optional support for mechanism reduction and tabulation.

See also
    Foam::chemistryModels::CCM

SourceFiles
    CCM_chemistryModelI.H
    CCM_chemistryModel.C

\*---------------------------------------------------------------------------*/
#define  STOP throw "error";
#define SHUTDOWN FatalErrorIn(__FUNCTION__) << " Shut Down By the User To Debug" << abort(FatalError);

// Timing macros for CCM performance analysis with parallel synchronization
#define CCM_SYNC() \
    if (Pstream::parRun()) { scalar sync = 1.0; reduce(sync, sumOp<scalar>()); }

#define CCM_TIMING_START(stepName, debugTimeFlag, timer) \
    if (debugTimeFlag) \
    { \
        CCM_SYNC(); \
        timer.cpuTimeIncrement(); \
    }

#define CCM_TIMING_END(stepName, debugTimeFlag, timer, chemistryModel) \
    if (debugTimeFlag) \
    { \
        CCM_SYNC(); \
        label stepIdx = (chemistryModel).getNextStepIndex(); \
        stepNames_.append(#stepName); \
        (chemistryModel).stepTimes_[stepIdx] = timer.cpuTimeIncrement(); \
    }


#ifndef CCM_chemistryModel_H
#define CCM_chemistryModel_H

#include "standard_chemistryModel.H"
#include "multicomponentMixture.H"
#include "ReactionList.H"
#include "DynamicField.H"
#include "CCMDebug.H"
#include "CombustionHelpers.H"
#include "atomicWeights.H"
#include "basicFastChemistryModel.H"
#include "CCMVars.H"
#include "MPHashTable.H"
#include "TCEntry.H"
#include "rateEntry.H"
#include "CCMHelpers.H"

#include <vector>
#include <mpi.h>
// #include "chemistryReductionMethod.H"
// #include "chemistryTabulationMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace chemistryModels
{

/*---------------------------------------------------------------------------*\
                           Class CCM Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class CCM
:
    public standard
{
    // Private classes

        //- Class to define scope of reaction evaluation. Runs pre-evaluate
        //  hook on all reactions on construction and post-evaluate on
        //  destruction.
        class reactionEvaluationScope
        {
            const CCM<ThermoType>& chemistry_;

        public:

            reactionEvaluationScope
            (
                const CCM<ThermoType>& chemistry
            )
            :
                chemistry_(chemistry)
            {
                forAll(chemistry_.reactions_, i)
                {
                    chemistry_.reactions_[i].preEvaluate();
                }
            }

            ~reactionEvaluationScope()
            {
                forAll(chemistry_.reactions_, i)
                {
                    chemistry_.reactions_[i].postEvaluate();
                }
            }
        };


    // Private data

        //- Switch to select performance logging
        Switch log_;

        //- Switch to enable per-cell CPU load caching for load-balancing
        Switch cpuLoad_;

        //- Type of the Jacobian to be calculated
        const jacobianType jacobianType_;

        //- Reference to the multi component mixture
        const multicomponentMixture<ThermoType>& mixture_;

        //- Thermodynamic data of the species
        const PtrList<ThermoType>& specieThermos_;

        //- Reactions
        const ReactionList<ThermoType> reactions_;

        //- List of reaction rate per specie [kg/m^3/s]
        PtrList<volScalarField::Internal> RR_;

        //- Temporary mass fraction field
        mutable scalarField Y_;

        //- Temporary simplified mechanism mass fraction field
        DynamicField<scalar> sY_;

        //- Temporary concentration field
        mutable scalarField c_;

        //- Temporary simplified mechanism concentration field
        DynamicField<scalar> sc_;

        //- Specie-temperature-pressure workspace fields
        mutable FixedList<scalarField, 5> YTpWork_;

        //- Specie-temperature-pressure workspace matrices
        mutable FixedList<scalarSquareMatrix, 2> YTpYTpWork_;

        //- Log file for average time spent solving the chemistry
        autoPtr<OFstream> cpuSolveFile_;

        //- The ratio of oxygen to carbon (mole) for the fuel
        scalar fuelOtoC_;

        //- Jh coefficient for each specie in the reaction mechanism
        scalarList JhCoeff_;

        //- Jc coefficient for each specie in the reaction mechanism
        scalarList JcCoeff_;

        //- Jn coefficient for each specie in the reaction mechanism
        scalarList JnCoeff_;

        //- Jo coefficient for each specie in the reaction mechanism
        scalarList JoCoeff_;

        //- Jh_h2o Jh coefficient for each specie excluding H2O
        scalarList Jh_h2oCoeff_;

        //- Jc_co2 Jc coefficient for each specie excluding CO2
        scalarList Jc_co2Coeff_;

        //- Jo_co2_h2o Jo coefficient for each specie excluding CO2 and H2O
        scalarList Jo_co2_h2oCoeff_;

        //- CCM dictionary
        dictionary CCMdict_;


        //- CCM Variable management class
        mutable autoPtr<CCMVars<ThermoType>> ccmVars_;

        //- CCM Debug helper class
        mutable autoPtr<CCMDebug<ThermoType>> ccmDebug_;

        //- CCM Combustion helper class
        mutable autoPtr<CombustionHelpers<ThermoType>> combustionHelpers_;

        //- A list of pricipal variables
        hashedWordList principalVars_;

        //- A list of principal variables to be updated
        hashedWordList updateVars_;


        //- The number of slices between min and max values for non-listed principal variables
        label nSlice_;

        //- The threshold to ignore for a min value
        scalar ignoreMin_;

        //- Several element mass fraction fields (Jh, Jc, Jn, Jo ...)
        volScalarField Jh_,Jc_,Jn_,Jo_,Jh_h2o_,Jc_co2_,Jo_co2_h2o_, J_;

        //- Equivalence ratio(phi) and scalar dissipation rate(chi)
        volScalarField phieq_, chi_;

        //- The maximum indices representations
        label maxiRepresentation_;

        //- Inter-core Communication
        label redistributeEvery_;
        bool redistributedThisStep_;
        labelList loads_;
        label localCores_;

        MPI_Comm chemComm_;

        //- Old style scalar dissipation rate
        Switch oldStyleChi_;

        //- Old style equivalence ratio
        Switch oldStylePhi_;

        // //- Acceleration Information output
        // OFstream accelerationInfo_;

        //- Temporary global maximum mass fraction field
        scalarField Ymax_;

        //- Flag to enable detailed timing analysis
        Switch debugTime_;

        //- Contiguous arena holding every cell's zone key (the "zoneIndex"
        //  bytes). Replaces the old List<string> zoneIndex_ so the keys
        //  live in one cache-friendly buffer instead of nCells separate
        //  heap allocations. Cell celli's key occupies
        //  [celli*keyStride_, celli*keyStride_ + keyLen_[celli]).
        List<char> keyArena_;

        //- Byte length of each cell's key inside keyArena_.
        labelList keyLen_;

        //- Fixed per-cell stride into keyArena_ (worst-case key length for
        //  the current step: nEnc*2 + 2 for the optional logP suffix).
        label keyStride_;

        //- remaidner of the reaction index (named zoneRemainder in the program)
        labelList zoneRemainder_;

        //- per-cell FNV hash of each cell's zone key, computed once in the
        //  grouping loop and reused for the mpReceivedRR_ lookup in
        //  Rate_Assignment (and for the mpToCore_ find/insert in grouping).
        List<uint> zoneHash_;

        //- Timing arrays for the solve function steps
        mutable scalarField stepTimes_;

        //- Step names corresponding to timing entries
        mutable wordList stepNames_;

        //- CPU timer for step timing
        mutable cpuTime stepTimer_;

        //- Current step index for automatic timing
        mutable label currentStepIndex_;

        //- EC mode settings
        Switch ecEnabled_;
        //wordList ecInitialVars_;
        label numECVarsToAdd_;
        label numECVarsToRemove_;
        label ecUpdateFreq_;
        label currentStep_;

        //- Initialization mode: manual or default
        word initMode_;

        //- Flag to treat pressure variations in high Mach cases
        Switch highMach_;

        //- MPHashTable-backed grouping hashtables. All six are allocated
        //  once in the ctor body (sizes derived from nCells/nCores and a
        //  per-block byte estimate). Each step clears them (memory pool
        //  re-used) before re-filling.
        mutable PtrList<MPHashTable> mpToCore_;
        mutable PtrList<MPHashTable> mpFromCore_;
        mutable autoPtr<MPHashTable> mpGathered_;
        mutable autoPtr<MPHashTable> mpLocalRR_;
        mutable PtrList<MPHashTable> mpReturnToCore_;
        mutable PtrList<MPHashTable> mpReceivedRR_;

        //- Send/recv staging buffers for the in-line MPI_Alltoallv pipeline.
        //  Sized once in the ctor to the strict upper bound (nCells *
        //  perBlockTC), so per-step we just memcpy into the prefix and
        //  treat the rest as scratch. perBlockRR < perBlockTC, so the same
        //  buffers also cover the Phase 3 return path.
        mutable std::vector<char> mpSendBuf_;
        mutable std::vector<char> mpRecvBuf_;

        //- Optimized ODE solver
        Switch optimizedODE_;
        const autoPtr<basicFastChemistryModel> fastChemistryPtr_;



        typedef scalarField (chemistryModels::CCM<ThermoType>::*GetRRFuncType)
        (
                scalarField& Y,
                scalar& T, 
                scalar& p,
                const scalar& deltaT,
                scalar& deltaTChem,
                const scalar& rho,
                const scalar& rho0
        );

        GetRRFuncType GetRRFunc_;


    // Private Member Functions

        //- Write access to chemical source terms
        //  (e.g. for multi-chemistry model)
        inline PtrList<volScalarField::Internal>& RR();

        //- Stage and exchange raw byte buffers via MPI_Alltoallv.
        //  Performs 2 preliminary MPI_Alltoall calls for per-peer byte
        //  and entry counts, packs the per-peer src into mpSendBuf_, then
        //  dispatches MPI_Alltoallv. On return, outRecvPtrs[p] points
        //  into mpRecvBuf_ at peer p's segment of length outRecvBytes[p]
        //  carrying outRecvEntries[p] entries.
        //  Serial mode (N==1 or not parRun) routes srcPtrs[0] through
        //  to outRecvPtrs[0] without copying.
        void mpAlltoallvBytes
        (
            int N,
            const void* const srcPtrs[],
            const int srcBytes[],
            const int srcEntries[],
            const char* outRecvPtrs[],
            int outRecvBytes[],
            int outRecvEntries[]
        ) const;

    // Private Member Functions

        //- Solve the reaction system for the given time step
        //  of given type and return the characteristic time
        //  Variable number of species added
        template<class DeltaTType>
        scalar solve(const DeltaTType& deltaT);


public:

    //- Runtime type information
    TypeName("CCM");


    // Constructors

        //- Construct from thermo
        CCM(const fluidMulticomponentThermo& thermo);

        //- Disallow default bitwise copy construction
        CCM(const CCM&) = delete;


    //- Destructor
    virtual ~CCM();


    // Member Functions

        // Access

            //- Return reference to the mixture
            inline const multicomponentMixture<ThermoType>& mixture() const;

            //- The reactions
            inline const PtrList<Reaction<ThermoType>>& reactions() const;

            //- Thermodynamic data of the species
            inline const PtrList<ThermoType>& specieThermos() const;


        // Overrides to chemistryModel functions

            //- The number of reactions
            virtual inline label nReaction() const;

            //- Return reaction rates of the species [kg/m^3/s]
            virtual inline const PtrList<volScalarField::Internal>& RR() const;

            //- Return the name of reactioni
            virtual inline const word& reactionName
            (
                const label reactioni
            ) const;

            //- Return the rate of reactioni [kmol/m^3/s]
            virtual tmp<volScalarField::Internal> reactionRR
            (
                const label reactioni
            ) const;

            //- Return reaction rates of the species in reactioni [kg/m^3/s]
            virtual PtrList<volScalarField::Internal> specieReactionRR
            (
                const label reactioni
            ) const;

            //- Calculates the reaction rates
            virtual void calculate();

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalar deltaT);

            //- Solve the reaction system for the given time step
            //  and return the characteristic time
            virtual scalar solve(const scalarField& deltaT);

            //- Return the chemical time scale
            virtual tmp<volScalarField> tc() const;

            //- Return the heat release rate [kg/m/s^3]
            virtual tmp<volScalarField> Qdot() const;

            //- Get reaction rates from Y, T, and p
            scalarField getRRGivenYTP
            (
                scalarField Y,
                scalar& T,
                scalar& p,
                const scalar& deltaT,
                scalar& deltaTChem,
                const scalar& rho,
                const scalar& rho0
            );


        // ODE solution functions

            //- Calculate the ODE derivatives
            virtual void derivatives
            (
                const scalar t,
                const scalarField& YTp,
                const label li,
                scalarField& dYTpdt
            ) const;

            //- Calculate the ODE jacobian
            virtual void jacobian
            (
                const scalar t,
                const scalarField& YTp,
                const label li,
                scalarField& dYTpdt,
                scalarSquareMatrix& J
            ) const;

            //- Solve the ODE system
            virtual void solve
            (
                scalar& p,
                scalar& T,
                scalarField& Y,
                const label li,
                scalar& deltaT,
                scalar& subDeltaT
            ) const;


            // Functions to distribute reactionEntry/solve/update the reaction rates
            void distributeReactionEntry(
                const hashedWordList& encodingVars = hashedWordList(),
                bool debug = false
            );

            void updateReactionRate();

            void distributeReactionRate();

            // //- Calculate the reaction rates give Y,T and P. 
            // // Note, density should be provided!
            // // deltaTChem will be modified!
            // scalarField getRRGivenYTP
            // (
            //     scalarField& Y,
            //     scalar& T, 
            //     scalar& p,
            //     const scalar& deltaT,
            //     scalar& deltaTChem,
            //     const scalar& rho,
            //     const scalar& rho0
            // );

            scalarField getRRGivenYTP_Optimized
            (
                scalarField& Y,
                scalar& T, 
                scalar& p,
                const scalar& deltaT,
                scalar& deltaTChem,
                const scalar& rho,
                const scalar& rho0
            );

            scalarField getRRGivenYTP_Basic
            (
                scalarField& Y,
                scalar& T, 
                scalar& p,
                const scalar& deltaT,
                scalar& deltaTChem,
                const scalar& rho,
                const scalar& rho0
            );
            
            // //- Reset timing step counter
            // void resetTimingSteps() const;

            //- Get next step index and increment counter
            label getNextStepIndex() const;

            //- update zoneIndex_ and zoneRemainder_ if mesh has changed
            void updateCCM4MeshChange();

            //- redistribute loads
            void redistributeLoads();


            // //- Add step name to the dynamic list
            // void addStepName(const word& stepName) const;



            // - calculate maximum standard deviation of all species inside the kinetics
            scalarField maxYstdNorm() const;


            // For encode and decode reaction index
            inline uint encode(const word& s);
            //- FNV-1a over raw key bytes (same hash as encode(word) and
            //  MPHashTable::encode); used by the contiguous key arena.
            inline uint encode(const char* data, label len);
            inline string toBase256Word(uint value);
            inline void listToBase256Word(const labelList& valueList, string& encodedString);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const CCM&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace chemistryModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CCM_chemistryModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CCM_chemistryModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
