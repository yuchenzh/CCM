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

#ifndef TCEntry_H
#define TCEntry_H

#include "fvCFD.H"
#include "MPHashTable.H"

namespace Foam
{

class TCEntry
{
    char* entryPtr_;
    label nSpecies_;

    // Offsets (bytes from entryPtr_) of every debug-independent field.
    // Ystd lives at OFF_debug_ + slot (only valid when getDebug() is true).
    label OFF_Y_, OFF_T_, OFF_p_, OFF_dt_, OFF_rho0_, OFF_rho_,
          OFF_count_, OFF_debug_;


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

        //- Offset of Ystd region (only meaningful if debug=true).
        label offYstd_() const
        {
            return OFF_debug_ + label(MemoryPool::alignUnit_);
        }

        void computeOffsets_()
        {
            auto* e = reinterpret_cast<MPHashTable::Entry*>(entryPtr_);
            const label slot = label(MemoryPool::alignUnit_);
            label off = MPHashTable::headerBytes()
                      + MemoryPool::alignUp(e->keyLength);
            OFF_Y_     = off;  off += nSpecies_ * slot;
            OFF_T_     = off;  off += slot;
            OFF_p_     = off;  off += slot;
            OFF_dt_    = off;  off += slot;
            OFF_rho0_  = off;  off += slot;
            OFF_rho_   = off;  off += slot;
            OFF_count_ = off;  off += slot;
            OFF_debug_ = off;
        }


public:

    // Static helpers

        //- User-content bytes after the Entry header. This is the value
        //  passed to MPHashTable::insertManual(key, size). The aligned
        //  key region + payload; the tail pad is added by insertManual
        //  via alignUp.
        static label userSize(label keyLen, label nSpecies, bool debug)
        {
            const label slot = label(MemoryPool::alignUnit_);
            const label nSlots = (debug ? 2*nSpecies : nSpecies) + 7;
            return MemoryPool::alignUp(keyLen) + nSlots * slot;
        }

        //- Backwards-compatible alias.
        static label insertSize(label keyLen, label nSpecies, bool debug)
        {
            return userSize(keyLen, nSpecies, debug);
        }


    // Constructors

        //- Attach to an existing MPHashTable block.
        TCEntry(void* entry, label nSpecies)
        :
            entryPtr_(static_cast<char*>(entry)),
            nSpecies_(nSpecies)
        {
            computeOffsets_();
        }


    // Initialise a fresh block. Ystd is written only when debug=true; the
    // caller must have requested a debug-sized block from insertManual.

        void initialise
        (
            const scalarField& Y,
            scalar T, scalar p, scalar dtChem,
            scalar rho0, scalar rho,
            label count,
            bool debug,
            const scalarField& Ystd = scalarField(0)
        )
        {
            const size_t bytes = size_t(nSpecies_) * sizeof(scalar);
            std::memcpy(arr_(OFF_Y_), Y.cdata(), bytes);
            slot_<scalar> (OFF_T_)     = T;
            slot_<scalar> (OFF_p_)     = p;
            slot_<scalar> (OFF_dt_)    = dtChem;
            slot_<scalar> (OFF_rho0_)  = rho0;
            slot_<scalar> (OFF_rho_)   = rho;
            slot_<int64_t>(OFF_count_) = count;
            slot_<int64_t>(OFF_debug_) = debug ? 1 : 0;

            if (debug)
            {
                scalar* s = arr_(offYstd_());
                if (Ystd.size() == nSpecies_)
                    std::memcpy(s, Ystd.cdata(), bytes);
                else
                    std::memset(s, 0, bytes);
            }
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
        scalarField getY() const
        {
            scalarField f(nSpecies_);
            const scalar* y = arr_(OFF_Y_);
            for (label i = 0; i < nSpecies_; ++i) f[i] = y[i];
            return f;
        }
        //- Returns empty scalarField when debug=false (matching the result
        //  of ReactionEntry's non-debug operator+).
        scalarField getYstd() const
        {
            if (!getDebug()) return scalarField(0);
            scalarField f(nSpecies_);
            const scalar* s = arr_(offYstd_());
            for (label i = 0; i < nSpecies_; ++i) f[i] = s[i];
            return f;
        }
        scalar getT()      const { return slot_<scalar>(OFF_T_); }
        scalar getP()      const { return slot_<scalar>(OFF_p_); }
        scalar getDtChem() const { return slot_<scalar>(OFF_dt_); }
        scalar getRho0()   const { return slot_<scalar>(OFF_rho0_); }
        scalar getRho()    const { return slot_<scalar>(OFF_rho_); }
        label  getN()      const { return label(slot_<int64_t>(OFF_count_)); }
        bool   getDebug()  const { return slot_<int64_t>(OFF_debug_) != 0; }

        //- Total bytes this block occupies, read from the Entry header
        //  (MPHashTable writes it on insertManual).
        label blockSize() const
        {
            return reinterpret_cast<MPHashTable::Entry*>(entryPtr_)->blockSize;
        }


    // Hot path: in-place merge — *this <- weighted_avg(*this, other).
    // Ystd is only touched when both sides are debug=true.

        void mergeInPlace(const TCEntry& other)
        {
            const scalar c1 = scalar(getN());
            const scalar c2 = scalar(other.getN());
            const scalar tot = c1 + c2;
            const scalar r1 = c1/tot;
            const scalar r2 = c2/tot;
            const bool dbg = getDebug() && other.getDebug();

            scalar* y1 = arr_(OFF_Y_);
            const scalar* y2 = other.arr_(other.OFF_Y_);

            if (dbg)
            {
                // Rare path — readability over allocation cost. Mirrors
                // ReactionEntry::operator+ Welford merge term-for-term.
                const scalarField Y1 = getY();
                const scalarField Y2 = other.getY();
                const scalarField S1 = (c1 == 1)
                    ? scalarField(nSpecies_, 0.0) : getYstd();
                const scalarField S2 = (c2 == 1)
                    ? scalarField(nSpecies_, 0.0) : other.getYstd();
                const scalar invSum = 1.0/(1.0/c1 + 1.0/c2);
                const scalarField M
                    = sqr(S1)*(c1 - 1) + sqr(S2)*(c2 - 1)
                    + sqr(Y1 - Y2)*invSum;
                const scalarField Yavg = Y1*r1 + Y2*r2;
                const scalarField Ystd = sqrt(M/(tot - 1));

                scalar* yDst = arr_(OFF_Y_);
                scalar* sDst = arr_(offYstd_());
                for (label i = 0; i < nSpecies_; ++i)
                {
                    yDst[i] = Yavg[i];
                    sDst[i] = Ystd[i];
                }
            }
            else
            {
                // Hot path — pointer loop, zero allocation, vectorisable.
                for (label i = 0; i < nSpecies_; ++i)
                {
                    y1[i] = y1[i]*r1 + y2[i]*r2;
                }
            }

            slot_<scalar>(OFF_T_)    = getT()      * r1 + other.getT()      * r2;
            slot_<scalar>(OFF_p_)    = getP()      * r1 + other.getP()      * r2;
            slot_<scalar>(OFF_dt_)   = getDtChem() * r1 + other.getDtChem() * r2;
            slot_<scalar>(OFF_rho0_) = getRho0()   * r1 + other.getRho0()   * r2;
            slot_<scalar>(OFF_rho_)  = getRho()    * r1 + other.getRho()    * r2;
            slot_<int64_t>(OFF_count_) = int64_t(tot);
            slot_<int64_t>(OFF_debug_) = dbg ? 1 : 0;
        }


    // Block ingestion

        //- Clear \c dst and repopulate it from a packed buffer of
        //  TCEntry blocks. Thin wrapper over MPHashTable::decodeFrom
        //  (the decode is byte-level and needs no TCEntry knowledge).
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
