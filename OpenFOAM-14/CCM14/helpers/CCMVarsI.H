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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
inline void Foam::CCMVars<ThermoType>::outputUnorderedSet
(
    const std::unordered_set<std::string>& vars
) const
{
    for (const auto& var : vars)
    {
        std::cout << var << " ";
    }
    std::cout << std::endl;    
}


template<class ThermoType>
inline bool Foam::CCMVars<ThermoType>::hasVar
(
    const std::unordered_set<std::string>& varSet,
    const Foam::word& var
) const
{
    return varSet.find(var) != varSet.end();
}


template<class ThermoType>
inline void Foam::CCMVars<ThermoType>::outputHashedWordList
(
    const Foam::hashedWordList& list
) const
{
    forAll(list, i)
    {
        Info << list[i] << " ";
    }
    Info << endl;
}

template<class ThermoType>
inline Foam::scalar Foam::CCMVars<ThermoType>::tolerance() const
{
    return 1.0/nSlice_;
}


template<class ThermoType>
inline void Foam::CCMVars<ThermoType>::reportVars() const
{
    Info << "====================================" << nl;
    Info << "non species variables:";
    outputUnorderedSet(nonSpVars_);

    Info << "listed variables:" ;
    outputUnorderedSet(listedVars_);

    Info << "adaptive species variables";
    outputUnorderedSet(adaptSpVars_);

    Info << "used adaptive species variables:";
    outputUnorderedSet(usedAdaptSpVars_);

    Info << "unused adaptive species varaibles:";
    outputUnorderedSet(unusedAdaptSpVars_);
    Info << "====================================" << nl;
    Info << endl;
}

template<class ThermoType>
inline Foam::hashedWordList Foam::CCMVars<ThermoType>::ecVars() const
{
    if (!ecMode_)
    {
        return principalVars_;
    }

    // ecVars = listedVars_ ∪ nonSpVars_ ∪ usedAdaptSpVars_
    //        = principalVars_ \ unusedAdaptSpVars_
    // Iterating principalVars_ (deterministic across ranks) instead of
    // the std::unordered_set members keeps ecVars order bit-identical
    // on every rank, which distributeReactionEntry's cell-grouping
    // encoding relies on.
    wordList ecVars(principalVars_.size() - unusedAdaptSpVars_.size());
    label count = 0;
    forAll(principalVars_, i)
    {
        const word& name = principalVars_[i];
        if (!hasVar(unusedAdaptSpVars_, name))
        {
            ecVars[count++] = name;
        }
    }
    return hashedWordList(ecVars);
}


template<class ThermoType>
inline void Foam::CCMVars<ThermoType>::updateRunningSpIndices()
{
    runningSpIndices_.clear();

    // When ecMode_ is off every adapt species sits in unusedAdaptSpVars_
    // by construction, so the "exclude unused" path would drop them all;
    // keep the two branches separate to preserve that intent.
    if (ecMode_ == false)
    {
        forAll(principalVars_, pi)
        {
            const word& var = principalVars_[pi];
            if (speciesList_.found(var))
            {
                runningSpIndices_.append(speciesList_[var]);
            }
        }
    }
    else
    {
        forAll(principalVars_, pi)
        {
            const word& var = principalVars_[pi];
            if (speciesList_.found(var) && !hasVar(unusedAdaptSpVars_, var))
            {
                runningSpIndices_.append(speciesList_[var]);
            }
        }
    }
}


template<class ThermoType>
inline const Foam::labelList& Foam::CCMVars<ThermoType>::runningSpIndices() const
{
    return runningSpIndices_;
}


template<class ThermoType>
inline const Foam::List<Foam::volScalarField*>& Foam::CCMVars<ThermoType>::runningNonSpFields() const
{
    return runningNonSpFields_;
}