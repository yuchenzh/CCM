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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
inline const Foam::multicomponentMixture<ThermoType>&
Foam::chemistryModels::CCM<ThermoType>::mixture() const
{
    return mixture_;
}


template<class ThermoType>
inline const Foam::PtrList<Foam::Reaction<ThermoType>>&
Foam::chemistryModels::CCM<ThermoType>::reactions() const
{
    return reactions_;
}


template<class ThermoType>
inline const Foam::PtrList<ThermoType>&
Foam::chemistryModels::CCM<ThermoType>::specieThermos() const
{
    return specieThermos_;
}


template<class ThermoType>
inline Foam::label
Foam::chemistryModels::CCM<ThermoType>::nReaction() const
{
    return reactions_.size();
}


template<class ThermoType>
inline const Foam::word&
Foam::chemistryModels::CCM<ThermoType>::reactionName
(
    const label reactioni
) const
{
    return reactions_[reactioni].name();
}


template<class ThermoType>
inline const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::fvMesh>>&
Foam::chemistryModels::CCM<ThermoType>::RR() const
{
    return RR_;
}


// FNV-1a hashing now lives in MPHashTable::encode; these forward to it so
// the table and every CCM call site hash identical bytes the same way.
template<class ThermoType>
inline uint Foam::chemistryModels::CCM<ThermoType>::encode(const word& s)
{
    return MPHashTable::encode(s.data(), label(s.size()));
}

template<class ThermoType>
inline uint Foam::chemistryModels::CCM<ThermoType>::encode(const char* data, label len)
{
    return MPHashTable::encode(data, len);
}

template<class ThermoType>
inline Foam::string Foam::chemistryModels::CCM<ThermoType>::toBase256Word(uint value)
{
    if (value > 49729) 
    {
        FatalErrorIn("CCMchemistryModel::toBase256Word")
            << "Value exceeds the maximum limit of 49728 (223*223-1)" << exit(FatalError);
    }
    
    char high = static_cast<char>(value / 223 + 32);  // high decimal 
    char low  = static_cast<char>(value % 223 + 32);  // low decimal

    string returnValue = (value > 222) ? (string() + high + low) : (string() + low);

    return returnValue;
}


template <class ThermoType>
inline void Foam::chemistryModels::CCM<ThermoType>::listToBase256Word
(
    const labelList& valueList,
    string& encodedString
)
{
    forAll(valueList, i)
    {
        uint value = valueList[i];
        if (value > 49729)
        {
            FatalErrorIn("CCMchemistryModel::listToBase256Word")
                << "Value at index " << i << " exceeds the maximum limit of 49728 (223*223-1)" << exit(FatalError);
        }

        char high = static_cast<char>(value / 223 + 32);  // high decimal
        char low  = static_cast<char>(value % 223 + 32);  // low decimal

        encodedString[2*i+0] = high;
        encodedString[2*i+1] = low;
    }
}

// ************************************************************************* //
