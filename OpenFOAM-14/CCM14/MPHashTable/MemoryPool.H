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

#ifndef MemoryPool_H
#define MemoryPool_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class MemoryPool Declaration
\*---------------------------------------------------------------------------*/

class MemoryPool
{
    // Private Data

        //- Pointer to the contiguous backing buffer
        char* pool_;

        //- Current capacity of the backing buffer, in bytes
        label size_;

        //- High-water mark: bytes handed out since the last clear()
        label used_;

        //- Growth increment used whenever the pool overflows.
        //  Captured from the constructor argument and therefore equal to
        //  the initial pool size.
        label enlargeUnit_;


public:

        //- Alignment unit, in bytes. All allocations are rounded up to a
        //  multiple of this value so the buffer stays aligned for 8-byte
        //  types (scalar, pointer, label on Int64). Public so callers
        //  (e.g. MPHashTable, TCEntry) can pad to the same unit.
        constexpr static size_t alignUnit_ = 8;


    // static helper functions
    template<typename T>
    static T* align(T* ptr, size_t alignment = alignUnit_)
    {
        std::uintptr_t p = reinterpret_cast<std::uintptr_t>(ptr);
        p = (p + alignment - 1) & ~(alignment - 1);
        return reinterpret_cast<T*>(p);
    }

    //- Round a byte count up to the next multiple of \c alignUnit_.
    //  Used by callers computing block sizes that must stay aligned.
    static label alignUp(label bytes)
    {
        constexpr label A = label(alignUnit_);
        return (bytes + A - 1) & ~(A - 1);
    }

    // Constructors

        //- Construct from initial pool size (also used as growth increment).
        //  \c poolSize must be at least \c alignUnit_ bytes.
        MemoryPool(label poolSize)
        :
            size_(poolSize),
            used_(0),
            enlargeUnit_(poolSize)
        {
            if (poolSize < label(alignUnit_))
            {
                FatalErrorInFunction
                    << "MemoryPool: initial size " << poolSize
                    << " smaller than alignment unit " << alignUnit_
                    << " bytes" << exit(FatalError);
            }

            // new char[N] is guaranteed to return memory aligned to
            // alignof(std::max_align_t) (16 bytes on x86-64), so the
            // buffer start already satisfies alignUnit_ - no extra work.
            pool_ = new char[size_];
        }


    //- Destructor
    ~MemoryPool()
    {
        delete[] pool_;
    }


    // Member Functions

        // Access

            //- Current capacity (bytes)
            const label& size() const
            {
                return size_;
            }

            //- Bytes handed out since the last clear()
            const label& used() const
            {
                return used_;
            }

            //- get pointer based on offset
            char* get(label offset) const
            {
                return pool_ + offset;
            }

            //- start
            char* start() const
            {
                return pool_;
            }


        // Allocation

            //- Return a pointer to \c bytes of storage inside the pool.
            //  Grows the pool by \c enlargeUnit_ on overflow.
            //  The request is rounded up to \c alignUnit_ so the next
            //  allocation also lands on an aligned boundary.
            //  \warning Previously returned pointers become dangling on
            //  growth - hold offsets instead if you need stable references.
            char* allocate(label bytes)
            {
                // Round request up so used_ stays a multiple of alignUnit_.
                const label A = label(alignUnit_);
                bytes = (bytes + A - 1) & ~(A - 1);

                if (used_ + bytes > size_)
                {
                    Info<< "Enlarging memory pool from " << size_
                        << " to " << size_ + enlargeUnit_ << endl;

                    // Grow-and-copy. The new buffer is fresh from new[],
                    // so its start is already aligned.
                    char* newPool = new char[size_ + enlargeUnit_];
                    std::memcpy(newPool, pool_, used_);
                    delete[] pool_;

                    pool_  = newPool;
                    size_ += enlargeUnit_;
                }

                char* ptr = pool_ + used_;
                used_    += bytes;
                return ptr;
            }

            //- Rewind the high-water mark without releasing memory.
            //  Bytes are left untouched and will be overwritten by the next
            //  allocate(). Cheap; intended for per-step reset.
            void clear()
            {
                used_ = 0;
            }

            //- Release the current buffer and re-allocate one of the initial
            //  size. Use to recover memory after a spike.
            void clearAndShrink()
            {
                delete[] pool_;
                pool_ = nullptr;

                size_ = enlargeUnit_;
                pool_ = new char[size_];
                used_ = 0;
            }

            //- Per-allocation free is not supported by a bump allocator.
            void deallocate()
            {
                NotImplemented;
            }


        // Debug

            //- Write current size and usage to Info
            void printInfo() const
            {
                Info<< "MemoryPool: size = " << size_
                    << ", used = " << used_
                    << ", free = " << size_ - used_
                    << endl;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
