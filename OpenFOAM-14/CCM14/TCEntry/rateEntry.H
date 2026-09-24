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

#ifndef rateEntry_H
#define rateEntry_H

#include "fvCFD.H"
#include "MPHashTable.H"

namespace Foam
{

class rateEntry
{
    char* entryPtr_;
    label nSpecies_;

    // Offsets (bytes from entryPtr_) of every field.
    label OFF_RR_, OFF_dt_;


    // Private helpers

        template<typename T>
        T& slot_(label off) const
        {
            return *reinterpret_cast<T*>(entryPtr_ + off);
        }

        scalar* arr_(label off) const
        {
            return reinterpret_cast<scalar*>(entryPtr_ + off);
        }

        void computeOffsets_()
        {
            auto* e = reinterpret_cast<MPHashTable::Entry*>(entryPtr_);
            const label slot = label(MemoryPool::alignUnit_);
            label off = MPHashTable::headerBytes()
                      + MemoryPool::alignUp(e->keyLength);
            OFF_RR_ = off;  off += nSpecies_ * slot;
            OFF_dt_ = off;
        }


public:

    // Static helpers

        //- User-content bytes after the Entry header. This is the value
        //  passed to MPHashTable::insertManual(key, size). The aligned
        //  key region + RR + dt; the tail pad is added by insertManual
        //  via alignUp.
        static label userSize(label keyLen, label nSpecies)
        {
            const label slot = label(MemoryPool::alignUnit_);
            return MemoryPool::alignUp(keyLen) + (nSpecies + 1) * slot;
        }

        //- Backwards-compatible alias for callers that already say
        //  insertSize(...). Same value as userSize(...).
        static label insertSize(label keyLen, label nSpecies)
        {
            return userSize(keyLen, nSpecies);
        }


    // Constructors

        rateEntry(void* entry, label nSpecies)
        :
            entryPtr_(static_cast<char*>(entry)),
            nSpecies_(nSpecies)
        {
            computeOffsets_();
        }


    // Initialise a fresh block.

        void initialise
        (
            const scalarField& RR,
            scalar deltaTChem
        )
        {
            std::memcpy
            (
                arr_(OFF_RR_), RR.cdata(),
                size_t(nSpecies_) * sizeof(scalar)
            );
            slot_<scalar>(OFF_dt_) = deltaTChem;
        }


    // Safe getters (deep copy)

        string getKey() const
        {
            auto* e = reinterpret_cast<MPHashTable::Entry*>(entryPtr_);
            return string
            (
                static_cast<const char*>(e->data()),
                size_t(e->keyLength)
            );
        }
        scalarField getRR() const
        {
            scalarField f(nSpecies_);
            const scalar* r = arr_(OFF_RR_);
            for (label i = 0; i < nSpecies_; ++i) f[i] = r[i];
            return f;
        }
        scalar getDeltaTChem() const { return slot_<scalar>(OFF_dt_); }

        //- Raw (non-owning) pointer to the nSpecies RR values inside the
        //  block. Lets the tiled RR assignment copy without allocating a
        //  scalarField per cell.
        const scalar* rawRR() const { return arr_(OFF_RR_); }

        //- Total bytes this block occupies, read from the Entry header
        //  (MPHashTable writes it on insertManual).
        label blockSize() const
        {
            return reinterpret_cast<MPHashTable::Entry*>(entryPtr_)->blockSize;
        }


    // Block ingestion

        //- Clear \c dst and repopulate it from a packed buffer of
        //  rateEntry blocks. Thin wrapper over MPHashTable::decodeFrom
        //  (the decode is byte-level and needs no rateEntry knowledge).
        static label decodeInto
        (
            MPHashTable& dst,
            const void* blocksStart,
            label nEntries
        )
        {
            return dst.decodeFrom(blocksStart, nEntries);
        }
};

} // End namespace Foam

#endif

// ************************************************************************* //
