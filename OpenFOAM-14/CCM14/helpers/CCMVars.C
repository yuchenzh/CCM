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

#include "CCMVars.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class ThermoType>
Foam::CCMVars<ThermoType>::CCMVars
(
    const multicomponentMixture<ThermoType>& mixture,
    const PtrList<ThermoType>& specieThermos,
    const fvMesh& mesh,
    const dictionary& CCMdict
):
mixture_(mixture),
specieThermos_(specieThermos),
mesh_(mesh),
CCMdict_(CCMdict),
initMode_(CCMdict_.lookupOrDefault<word>("mode", "default")),
ecMode_(false),
nSlice_(CCMdict_.lookupOrDefault<label>("nSlice", 50)),
principalVars_(CCMdict_.lookupOrDefault<hashedWordList>("principalVars", hashedWordList())),
speciesList_(mixture_.species()),
ignoreMin_(CCMdict_.lookupOrDefault<scalar>("ignoreMin", 1e-12))
{

    // initialize EC mode, default to be false
    if (CCMdict_.isDict("ecMode"))
    {
        ecMode_ = CCMdict_.subDict("ecMode").lookupOrDefault<Switch>("enabled", false);
    }
    if (initMode_ ==  "default")
    {
        ecMode_ = true;
    }

    //initialization of principalVars first
    if (initMode_ == "default")
    {
        principalVars_.clear();
        forAll(speciesList_, i)
        {
            principalVars_.append(speciesList_[i]);
        }
        principalVars_.append("T");
    }

    // Safety floor: ecVars must always contain at least one variable so
    // that removeN can never empty it. T is a non-species, so it lands in
    // nonSpVars_ and is never touched by add/remove.
    if (!principalVars_.found("T"))
    {
        WarningInFunction
            << "principalVars_ does not contain T; appending automatically "
               "so that ecVars always retains at least one variable."
            << endl;
        principalVars_.append("T");
    }

    // output
    Info << "Principal Variables:" << endl;
    forAll(principalVars_,i)
    {
        Info << principalVars_[i] << " ";
        if (i%8 == 0 && i != 0)
        {
            Info << nl;
        }
    }
    Info << endl;
    Info << "nSlice     = " << nSlice_ << endl;
    Info << "ignoreMin  = " << ignoreMin_ << endl;


        // Initialize listed variables
    dictionary pvDict(CCMdict_.subOrEmptyDict("pvInfo"));
    List<Foam::keyType> keys = pvDict.keys();
    forAll(keys, i)
    {
        if (principalVars_.found(keys[i]))
        {
            listedVars_.insert(keys[i]);
        }
        else
        {
            FatalErrorIn("chemistryModel/helpers/CCMVars()")
                << "The variable " << keys[i] << " in pvInfo is not found in principalVars_"
                << exit(FatalError);
        }
    }
    
    // Partition the remainder of principalVars_:
    //   - species not in pvInfo  -> adaptSpVars_ (and unused initially)
    //   - non-species not in pvInfo -> nonSpVars_ (always active)
    // listed vars are already handled above and are skipped here.
    forAll(principalVars_,pi)
    {
        word var = principalVars_[pi];
        if (listedVars_.find(var) != listedVars_.end())
        {
            continue;
        }
        else if (speciesList_.found(var))
        {
            adaptSpVars_.insert(var);
            unusedAdaptSpVars_.insert(var);
        }
        else
        {
            nonSpVars_.insert(var);
        }
    }

    // initialize runningSpIndices_
    updateRunningSpIndices();

    // initialize runningNonSpFields_. It has a fixed size.
    forAll(principalVars_, pi)
    {
        const word& var = principalVars_[pi];
        if (!speciesList_.found(var))
        {
            runningNonSpFields_.append(&mesh_.lookupObjectRef<volScalarField>(var));
        }
    }

    // varIndex_[name] = position of `name` in principalVars_.
    // principalVars_ is the same on every rank (read from dict or built
    // from speciesList_), so this index map is bit-identical across ranks
    // and safe to use as the canonical slot for min/max/span.
    varIndex_.reserve(principalVars_.size() * 2);
    forAll(principalVars_, i)
    {
        varIndex_[principalVars_[i]] = static_cast<std::size_t>(i);
    }

    // Allocate min/max/span aligned with varIndex_; sentinel values so a
    // missed update is loudly wrong, not silently zero.
    minvals_.setSize(principalVars_.size(),  GREAT);
    maxvals_.setSize(principalVars_.size(), -GREAT);
    spans_.setSize(principalVars_.size(),    0.0);

    // Pre-fill listed slots from pvInfo. min/max here are just initial
    // estimates (updateMinMaxSpan will overwrite them every step), but
    // spans_ for listed vars stays pinned forever.

    if (listedVars_.size() > 0)
    {
       
        Info << "------------------------------------" << endl;
        Info << "Listed variables" << endl;
        Info << setw(15) << "species" 
        << setw(15) << "min" 
        << setw(15) << "max" 
        << setw(15) << "span" 
        << endl;
        for (const auto& var : listedVars_)
        {
            List<scalar> vals = pvDict.lookup(word(var));
            if (vals.size() != 3)
            {
                FatalErrorIn("CCMVars constructor")
                    << "pvInfo entry for " << var
                    << " must be (min max span)" << exit(FatalError);
            }
            std::size_t idx = varIndex_.at(var);
            minvals_[idx] = vals[0];
            maxvals_[idx] = vals[1];
            spans_[idx]   = vals[2];
            Info << setw(15) << var 
            << setw(15) << fixed << setprecision(2) << vals[0] 
            << setw(15) << fixed << setprecision(2) << vals[1] 
            << setw(15) << fixed << setprecision(2) << vals[2] 
            << endl;
        }
        Info << "------------------------------------" << endl;
    }
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
namespace Foam
{

    template<class ThermoType>
    labelList CCMVars<ThermoType>::sort(const scalarList& importance) const
    {
        if (importance.size() != speciesList_.size())
        {
            FatalErrorIn("CCMVars::sort")
                << "Importance size " << importance.size()
                << " does not match species list size " << speciesList_.size() << nl
                << "And this function is only intended to sort for the whole kinetics"
                << exit(FatalError);
        }
        labelList indices(importance.size());
        forAll(indices, i) { indices[i] = i; }
        std::sort(indices.begin(), indices.end(),
                    [&](label a, label b) { return importance[a] > importance[b]; });
        return indices;
    }


    // True iff no unused adaptSp species is above tolerance. Replaces the
    // VariableRanking::underControl flag used by the legacy EC loop.
    //
    // Precondition: sortSpeciesByStdNorm(Ystdnorm) was called with this same
    // Ystdnorm. Then sortedUnusedSpIndices_ is sorted descending by std, so
    // the first entry is the worst unused species -- if even it is within
    // tolerance, all the rest are too.
    template<class ThermoType>
    bool CCMVars<ThermoType>::ecOkay(const scalarField& Ystdnorm) const
    {
        if (sortedUnusedSpIndices_.empty()) return true;
        return Ystdnorm[sortedUnusedSpIndices_[0]] <= tolerance();
    }


    // Sort *all* species by their (normalized) std, then split the sorted
    // sequence into two streams: those currently in usedAdaptSpVars_ and
    // those in unusedAdaptSpVars_. Both result lists preserve the global
    // importance ordering, so addN() can pick from the top of `unused`
    // and removeN() from the bottom of `used`.
    template<class ThermoType>
    void  CCMVars<ThermoType>::sortSpeciesByStdNorm(const scalarList& stdvalues)
    {
        // Sort indices [0, nSpecie) by stdvalues descending. Both inputs
        // (the const initial array and stdvalues itself) are identical on
        // every rank, so the sorted output is too.
        labelList sortedIndices = sort(stdvalues);

        sortedUsedSpIndices_.setSize(usedAdaptSpVars_.size());
        sortedUnusedSpIndices_.setSize(unusedAdaptSpVars_.size());

        label usedCount = 0;
        label unusedCount = 0;

        forAll(sortedIndices, i)
        {
            label spIdx = sortedIndices[i];
            word speciesName = speciesList_[spIdx];
            if (hasVar(usedAdaptSpVars_, speciesName))
            {
                sortedUsedSpIndices_[usedCount] = spIdx;
                // usedStd[usedCount] = stdvalues[spIdx];
                usedCount++;
            }
            else if (hasVar(unusedAdaptSpVars_, speciesName))
            {
                sortedUnusedSpIndices_[unusedCount] = spIdx;
                // unusedStd[unusedCount] = stdvalues[spIdx];
                unusedCount++;
            }
        }

        Info << "\nCCM Species Sorting by Standard Deviation:" << endl;
        Info << "====================================" << endl;

        // --- Used: 2 columns ---
        Info << setw(15) << "Used Species"
             << setw(15) << "std" << endl;
        Info << "------------------------------------" << endl;
        forAll(sortedUsedSpIndices_, i)
        {
            label spIdx = sortedUsedSpIndices_[i];
            word speciesName = speciesList_[spIdx];
            Info << setw(15) << speciesName
                 << setw(15) << scientific << setprecision(4) << stdvalues[spIdx]
                 << endl;
        }
        Info << "------------------------------------" << endl;
        Info << nl;

        // --- Unused: 3 columns with Status ---
        Info << setw(15) << "Unused Species"
             << setw(15) << "std/max (%)"
             << setw(20) << "Status" << endl;
        Info << "------------------------------------" << endl;
        const scalar tol = tolerance();
        forAll(sortedUnusedSpIndices_, i)
        {
            label spIdx = sortedUnusedSpIndices_[i];
            word speciesName = speciesList_[spIdx];
            const word status =
                (stdvalues[spIdx] > tol) ? word("Candidate")
                                         : word("Within tolerance");
            Info << setw(15) << speciesName
                 << setw(15) << fixed << setprecision(2) << stdvalues[spIdx]*100 
                 << setw(20) << status
                 << endl;
        }
        Info << "------------------------------------" << endl;
    }



    // Promote at most n unused species into used. Walk
    // sortedUnusedSpIndices_ from the front (most important first) and
    // stop as soon as a species falls below tolerance() -- because the
    // list is sorted descending, everything after it is also below.
    template<class ThermoType>
    wordList CCMVars<ThermoType>::addN(label n, const scalarField& Ymaxstdnorm)
    {
        // sort Ymaxstd to update sortedUsedSpIndices_ and sortedUnusedSpIndices_
        sortSpeciesByStdNorm(Ymaxstdnorm);

        // initialize the list.
        wordList varsToAdd;

        // determine the numbers of vars to add
        label numToAdd = min(n, sortedUnusedSpIndices_.size());

        // check whether all the variabels are above tolerance
        for (label i = 0; i < numToAdd; i++)
        {
            label spIdx = sortedUnusedSpIndices_[i];
            if (Ymaxstdnorm[spIdx] > tolerance())
            {
                varsToAdd.append(speciesList_[spIdx]);
            }
            else
            {
                break;
            }
        }
        
        // output info
        Info << nl << "Adding " << varsToAdd.size() << " variables:";
        forAll(varsToAdd, i)
        {
            Info << varsToAdd[i] << " ";
        }
        Info << endl;

        // execute the addition
        addVars(varsToAdd);

        return varsToAdd;
    }

    // Demote at most n used species back to unused. Walk
    // sortedUsedSpIndices_ from the END (least important first) and stop
    // once a species rises to/above tolerance() -- everything after it
    // (more important) must stay.
    template<class ThermoType>
    wordList CCMVars<ThermoType>::removeN(label n, const scalarField& Ymaxstdnorm)
    {
        // sort Ymaxstd to update sortedUsedSpIndices_ and sortedUnusedSpIndices_
        sortSpeciesByStdNorm(Ymaxstdnorm);

        // initialize the list
        wordList varsToRemove;

        // determine the maximum vars to remove
        label numToRemove = min(n, sortedUsedSpIndices_.size());

        // check if there is any variable that is above tolerance
        // if so cannot remove
        for (label i = 0; i < numToRemove; i++)
        {
            label spIdx = sortedUsedSpIndices_[sortedUsedSpIndices_.size() - 1 - i];
            if (Ymaxstdnorm[spIdx] < tolerance())
            {
                varsToRemove.append(speciesList_[spIdx]);
            }
            else
            {
                break;
            }
        }

        // info output
        Info << nl << "Removing " << varsToRemove.size() << " variables:";
        forAll(varsToRemove, i)
        {
            Info << varsToRemove[i] << " ";
        }
        Info << endl;

        // execute the removal
        removeVars(varsToRemove);

        return varsToRemove;
    }


    // Report the real field min/max for listed variables, so the user can
    // sanity-check their pvInfo-pinned ranges against actual flow conditions.
    template<class ThermoType>
    void CCMVars<ThermoType>::reportMinMax() const
    {
        if (listedVars_.empty())
        {
            return;
        }

        Info << "\nCCM Listed Variable Analysis:" << endl;
        Info << "====================================" << endl;
        Info << setw(15) << "Variable"
                << setw(15) << "Min Value"
                << setw(15) << "Max Value"
                << endl;

        // Iterate principalVars_ for deterministic order across ranks.
        forAll(principalVars_, pi)
        {
            const word& varName = principalVars_[pi];
            if (listedVars_.find(varName) == listedVars_.end())
            {
                continue;
            }
            const volScalarField& field =
                mesh_.lookupObjectRef<volScalarField>(varName);
            Info << setw(15) << varName
                    << setw(15) << scientific << setprecision(4) << gMin(field)
                    << setw(15) << scientific << setprecision(4) << gMax(field)
                    << endl;
        }
        Info << "====================================" << endl;
    }

    template <class ThermoType>
    void CCMVars<ThermoType>::addVars(const wordList& varsToAdd)
    {
        forAll(varsToAdd, i)
        {
            word var = varsToAdd[i];
            if (hasVar(adaptSpVars_, var))
            {
                if (!hasVar(usedAdaptSpVars_, var))
                {
                    usedAdaptSpVars_.insert(var);
                    unusedAdaptSpVars_.erase(var);
                }
                else
                {
                    FatalErrorIn("CCMVars::addVars")
                        << "Trying to add variable " << var << " that is already in usedAdaptSpVars_"
                        << exit(FatalError);
                }
            }
        }
        updateRunningSpIndices();
    }

    template<class ThermoType>
    void CCMVars<ThermoType>::removeVars(const wordList& varsToRemove)
    {
        forAll(varsToRemove, i)
        {
            word var = varsToRemove[i];
            if (hasVar(adaptSpVars_, var))
            {
                if (hasVar(usedAdaptSpVars_, var))
                {
                    usedAdaptSpVars_.erase(var);
                    unusedAdaptSpVars_.insert(var);
                }
                else
                {
                    FatalErrorIn("CCMVars::removeVars")
                        << "Trying to remove variable " << var << " that is not in usedAdaptSpVars_"
                        << exit(FatalError);
                }
            }
        }
        updateRunningSpIndices();
    }


    // Refresh minvals_/maxvals_/spans_ for non-listed vars. Listed vars
    // keep their pvInfo-pinned values.
    //
    // Performance: a naive impl would do 2*N Allreduces (one for each
    // gMin/gMax). Instead we compute local min and -max for every var,
    // pack them into one scalarList, and do a single minOp reduction --
    // a single Allreduce regardless of how many variables we track.
    template<class ThermoType>
    void CCMVars<ThermoType>::updateMinMaxSpan()
    {
        const label N = static_cast<label>(minvals_.size());

        // localPack layout: [min_0, -max_0, min_1, -max_1, ...]
        // Negating max lets us reduce both with one minOp:
        //   global_min = minOp(local_min)
        //   global_max = -minOp(-local_max)
        scalarList localPack(2 * N, GREAT);

        forAll(principalVars_, pi)
        {
            const word var = principalVars_[pi];
            const std::size_t i = varIndex_.at(var);

            const volScalarField& f =
                mesh_.lookupObjectRef<volScalarField>(var);
            const Foam::scalarField& F = f.primitiveField();

            localPack[2 * i]     =  Foam::min(F);
            localPack[2 * i + 1] = -Foam::max(F);
        }

        // Single element-wise min reduction covers every variable.
        // (minOp<scalarList> does not exist: there is no element-wise
        // min(List,List) in OpenFOAM, so use the list combine helpers.)
        Pstream::listCombineGather(localPack, minEqOp<scalar>());
        Pstream::listCombineScatter(localPack);

        forAll(principalVars_, pi)
        {
            const word var = principalVars_[pi];
            const std::size_t i = varIndex_.at(var);


            // listedVars keep the pvInfo-pinned min/max/span (set at construct
            // time); adaptSpVars and nonSpVars get min/max/span refreshed
            // dynamically from the field every step.
            if (hasVar(adaptSpVars_, var) || hasVar(nonSpVars_, var))
            {
                minvals_[i] = localPack[2 * i];
                // Floor max at min + ignoreMin_ so consumers that do (max - min)
                // directly never hit zero
                maxvals_[i] = max(-localPack[2 * i + 1], minvals_[i] + ignoreMin_);
                spans_[i]   = max((maxvals_[i] - minvals_[i]) / nSlice_, ignoreMin_ / 10.0);
            }
        }
    }

}