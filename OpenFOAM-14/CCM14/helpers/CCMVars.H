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

#ifndef CCMVars_H
#define CCMVars_H

#include "hashedWordList.H"
#include "dictionary.H"
#include "multicomponentMixture.H"
#include "volFields.H"
#include "Switch.H"
#include "scalarList.H"
#include "fvMesh.H"
#include <unordered_set>
#include <string>
#include <unordered_map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                        Class CCMVars Declaration
\*---------------------------------------------------------------------------*/
template<class ThermoType>
class CCMVars
{
    // Private data

        //- Reference to the multi component mixture
        const multicomponentMixture<ThermoType>& mixture_;

        //- Reference to thermodynamic data of the species
        const PtrList<ThermoType>& specieThermos_;

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to CCMdict_
        const dictionary& CCMdict_;

        //- initMode
        const word initMode_;

        //- ecMode enabled or not
        Switch ecMode_;

        //- number of slices
        label nSlice_;

        //- specified principal variables. principalVars_ = listedVars + nonSpVars + adaptSpVars
        hashedWordList principalVars_;

        //- listed variables in the pvInfo dictionnary. Min/max/span should not be updated
        //- listed variables can include species, whose min/max/span are specified in pvInfo
        std::unordered_set<std::string> listedVars_;

        //- non_species variable list. Min/max/span should be updated.
        std::unordered_set<std::string> nonSpVars_;

        //- adaptive species variables that are not listed in the pvInfo dictionary. adaptSpVars_ = usedAdaptSpVars_ + unusedAdaptSpVars_
        //- their min/max/span should be updated.
        std::unordered_set<std::string> adaptSpVars_;

        //- used adaptive species variable list
        std::unordered_set<std::string> usedAdaptSpVars_;

        //- unused adaptive species variable list. 
        //- if ecMode_ is disabled, all adaptSpVars_ will be in unusedAdaptSpVars_ and usedAdaptSpVars_ will be empty.
        std::unordered_set<std::string> unusedAdaptSpVars_;

        //- specieslist (for the ease of implementation)
        const hashedWordList& speciesList_;

        //- all running species (in listedVars_ and adaptSpVars_) index in speciesList_ (the same as the mechanism in OpenFOAM)
        labelList runningSpIndices_;

        //- Pointer cache for every non-species variable in principalVars_,
        //  in principalVars_ order. Non-species vars never migrate in/out of
        //  the running set, so this list is filled once at construction.
        List<volScalarField*> runningNonSpFields_;

        //- Name -> position in principalVars_.
        //  Filled at construction by  varIndex_[principalVars_[i]] = i,
        //  so varIndex_[name] returns the slot of `name` in
        //  minvals_/maxvals_/spans_ (and in any future array sized to
        //  principalVars_.size()).
        std::unordered_map<std::string, std::size_t> varIndex_;

        //- Per-variable min/max/span, all sized to principalVars_.size().
        //  Slot i corresponds to principalVars_[i] (see varIndex_).
        scalarList minvals_;
        scalarList maxvals_;
        scalarList spans_;

        //- ignoreMin for span floor
        scalar ignoreMin_;

        //- Importance-sorted indices into speciesList_ (NOT principalVars_).
        //  Filled by sortSpeciesByStdNorm(); consumed by addN()/removeN().
        //  used: most-important-first; unused: same ordering.
        labelList sortedUsedSpIndices_;
        labelList sortedUnusedSpIndices_;


    public:
    // Constructors

        //- Construct from parameters (only constructor).
        //  Builds principalVars_, partitions it into listed/nonSp/adaptSp,
        //  initializes varIndex_ and pre-fills listed slots from pvInfo.
        CCMVars
        (
            const multicomponentMixture<ThermoType>& mixture,
            const PtrList<ThermoType>& specieThermos,
            const fvMesh& mesh,
            const dictionary& CCMdict
        );


    // Membership update (collective: same args on every rank)

        //- Move vars from unusedAdaptSpVars_ into usedAdaptSpVars_.
        //  Silently ignores vars not in adaptSpVars_; fatal if a var is
        //  already in usedAdaptSpVars_.
        void addVars(const wordList& varsToAdd);

        //- Inverse of addVars(): move vars from used back to unused.
        void removeVars(const wordList& varsToRemove);


    // Min/max/span maintenance

        //- Refresh minvals_/maxvals_/spans_ for adaptSp+nonSp from the
        //  current volScalarField values. Uses one batched Allreduce
        //  packing every variable's local min/max into one scalarList.
        //  Listed vars keep their pvInfo-pinned values and are skipped.
        //  Collective.
        void updateMinMaxSpan();

        //- Print the *current field* gMin/gMax for listed vars (so the
        //  user can sanity-check pvInfo against actual flow ranges).
        //  Diagnostic only; does NOT update minvals_/maxvals_.
        void reportMinMax() const;


    // Adaptive variable selection (collective: caller must ensure
    // Ymaxstdnorm is bit-identical on every rank, typically by passing
    // a value already reduced via maxOp)

        //- Promote up to n unused species into used, in importance order,
        //  stopping early when importance falls below tolerance().
        //  Returns the names that were actually added.
        wordList addN(label n, const scalarField& Ymaxstdnorm);

        //- Demote up to n used species back to unused, from least to
        //  most important, stopping once importance reaches tolerance().
        //  Returns the names that were actually removed (so the caller
        //  can roll back if removal breaks ecOkay).
        wordList removeN(label n, const scalarField& Ymaxstdnorm);

        //- Sort all species by stdvalues (descending) and split the
        //  result into sortedUsedSpIndices_ / sortedUnusedSpIndices_.
        //  stdvalues must be sized to speciesList_.size().
        void sortSpeciesByStdNorm(const scalarList& stdvalues);

        //- Return species-list indices sorted by importance (descending).
        // The importance is maxstd/maxvalue (global) of each species
        //  importance.size() must equal speciesList_.size(); fatal otherwise.
        labelList sort(const scalarList& importance) const;

        //- True iff every unused adaptSp species has Ystdnorm <= tolerance().
        //  Caller must ensure Ystdnorm is bit-identical on every rank
        //  (e.g. already reduced via maxOp). No MPI inside.
        bool ecOkay(const scalarField& Ystdnorm) const;


    // Read-only accessors (used by distribute to do its own indexing
    // over a subset of principalVars_)

        const std::unordered_map<std::string, std::size_t>& varIndex() const
        {
            return varIndex_;
        }
        const scalarList& minvals() const { return minvals_; }
        const scalarList& maxvals() const { return maxvals_; }
        const scalarList& spans()   const { return spans_; }

        const std::unordered_set<std::string>& usedAdaptSpVars() const
        { return usedAdaptSpVars_; }
        const std::unordered_set<std::string>& unusedAdaptSpVars() const
        { return unusedAdaptSpVars_; }


    // Small helpers (defined in CCMVarsI.H)

        //- Debug: stream contents of a string set to stdout.
        void outputUnorderedSet(const std::unordered_set<std::string>& vars) const;

        //- Stream a hashedWordList space-separated, terminated by endl.
        void outputHashedWordList(const hashedWordList& list) const;

        //- O(1) membership check via std::unordered_set::find.
        bool hasVar(const std::unordered_set<std::string>& varSet, const word& var) const;

        //- 1 / nSlice_ ; the cut-off used by addN/removeN.
        scalar tolerance() const;

        //- Print all four/five variable sets to Info for diagnostics.
        void reportVars() const;

        //- Union of listed_ ∪ nonSp_ ∪ usedAdaptSp_, returned as a
        //  hashedWordList. In non-EC mode just returns principalVars_.
        hashedWordList ecVars() const;

        //- Rebuild runningSpIndices_ from the current usedAdaptSp/listed
        //  partition. Must be called whenever membership changes
        //  (constructor, addVars, removeVars).
        void updateRunningSpIndices();

        //- Read-only access to the cached running-species index list.
        const labelList& runningSpIndices() const;

        //- Return runningNonSpFields_ as a list of volScalarField pointers.
        const List<volScalarField*>& runningNonSpFields() const;

        //- Return principalVars_ as a hashedWordList.
        const hashedWordList& principalVars() const { return principalVars_; }

        //- Return whether ecMode is enabled.
        bool ecMode() const { return ecMode_; }
};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline functions

#include "CCMVarsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include implementation only for templates
#ifdef NoRepository
    #include "CCMVars.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //