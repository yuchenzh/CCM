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

#ifndef FASTHASHTABLE_H
#define FASTHASHTABLE_H
#include "MemoryPool.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
/*---------------------------------------------------------------------------*\
                         Class MPHashTable Declaration
\*---------------------------------------------------------------------------*/

class MPHashTable
{

public:
struct Entry
{
    //- Number of bytes in the key (key data follows the header).
    label keyLength;

    //- Total footprint of this block in bytes — alignUp(sizeof(Entry))
    //  + alignUp(userSize). Written by insertManual; lets external
    //  consumers walk a packed buffer of blocks without re-deriving
    //  sizes from the payload type.
    label blockSize;

    //- Bucket-chain link, byte offset into MemoryPool. -1 = end of chain.
    label next = -1;

    //- Pads sizeof(Entry) up to a multiple of alignUnit_ (8) so the key
    //  immediately following the header lands on an aligned offset.
    //  With WM_LABEL_SIZE=32 a 3-label Entry would be 12 B (not a
    //  multiple of 8); the 4th field brings it to 16 B.
    label pad_ = 0;


    // Helpers — convert between Entry* and the user-data pointer that
    // sits immediately after the (alignUp'd) header.

        //- Pointer to the user content (key + payload), i.e. \c this
        //  advanced by alignUp(sizeof(Entry)).
        void* data()
        {
            return reinterpret_cast<char*>(this)
                 + MemoryPool::alignUp(label(sizeof(Entry)));
        }

        const void* data() const
        {
            return reinterpret_cast<const char*>(this)
                 + MemoryPool::alignUp(label(sizeof(Entry)));
        }

        //- Inverse of data(): recover the Entry header given a pointer
        //  to the user content of the same block.
        static Entry* fromData(void* d)
        {
            return reinterpret_cast<Entry*>
            (
                static_cast<char*>(d)
              - MemoryPool::alignUp(label(sizeof(Entry)))
            );
        }
};

private:

    // Private Data

        //- Pointer to the contiguous backing buffer
        MemoryPool memory_;

        //- buckets to store the first element of each hash chain
        labelList buckets_;

        //- number of buckets
        label nBuckets_;

        //- number of Entries
        mutable label nEntries_ = 0;



    //- Raw chain walk: caller supplies the key bytes, its length, and a
    //  precomputed FNV hash. This is the single core used by every other
    //  lookup/find/insert overload — it lets hot paths (grouping / RR
    //  assignment) point straight at a contiguous key arena, with no
    //  std::string built per cell and no re-hashing.
    Entry* lookup
    (
        const char* keyData,
        label keyLen,
        label& bucketIndex,
        uint hash
    )
    {
        bucketIndex = label(hash % uint(nBuckets_));
        label off = buckets_[bucketIndex];
        while (off != -1)
        {
            Entry* entry  = reinterpret_cast<Entry*>(memory_.get(off));
            char*  keyPtr = static_cast<char*>(entry->data());
            if (entry->keyLength == keyLen
             && std::memcmp(keyPtr, keyData, size_t(keyLen)) == 0)
            {
                return entry;
            }
            off = entry->next;
        }
        return nullptr;
    }

    Entry* lookup(const string& key, label& bucketIndex, uint hash)
    {
        return lookup(key.data(), label(key.size()), bucketIndex, hash);
    }

    Entry* lookup(const string& key, label& bucketIndex)
    {
        return lookup(key, bucketIndex, encode(key));
    }



public:

    //- Footprint of the Entry header inside a block (always
    //  alignUp(sizeof(Entry))). The key starts at this offset.
    static label headerBytes()
    {
        return MemoryPool::alignUp(label(sizeof(Entry)));
    }

    // static hash function (FNV-1a). Canonical implementation lives here;
    // CCMchemistryModel::encode just forwards to these so the table and
    // every caller hash identical bytes the same way.
    static inline uint encode(const char* data, label len)
    {
        const uint FNV_prime     = 16777619u;     // 0x01000193
        const uint offset_basis  = 2166136261u;   // 0x811C9DC5

        uint hash = offset_basis;

        for (label i = 0; i < len; ++i)
        {
            hash ^= static_cast<unsigned char>(data[i]);  // XOR with byte
            hash *= FNV_prime;                            // mul by FNV prime
        }

        return hash;
    }

    static inline uint encode(const string& s)
    {
        return encode(s.data(), label(s.size()));
    }



    // Constructors

        MPHashTable(label nBuckets, label initialPoolSize)
        :
            memory_(initialPoolSize),
            buckets_(nBuckets, -1),
            nBuckets_(nBuckets)
        {}

//- Destructor
    ~MPHashTable(){}

    // Access

        //- Pointer to the first block in the pool. Use together with
        //  nEntries() and blockSize headers to walk the table linearly.
        char* poolStart() const { return memory_.start(); }

        //- Total bytes handed out since the last clear().
        label poolUsed()  const { return memory_.used(); }

    //- Member functions
    void* insertManual(const string& key, label size)
    {
        return insertManual(key.data(), label(key.size()), size, encode(key));
    }

    //- insertManual variant taking a precomputed FNV hash for the key,
    //  so the caller can avoid re-hashing.
    void* insertManual(const string& key, label size, uint hash)
    {
        return insertManual(key.data(), label(key.size()), size, hash);
    }

    //- Raw insertManual core: key bytes + length straight from a caller
    //  arena, with a precomputed FNV hash. No std::string involved.
    void* insertManual
    (
        const char* keyData,
        label keyLen,
        label size,
        uint hash
    )
    {
        label bucketIndex;
        if (Entry* existing = lookup(keyData, keyLen, bucketIndex, hash))
        {
            Info << "Key already exists, skipping insert" << endl;
            return existing;
        }

        // alignUp(header) + alignUp(user) — independently rounded so the
        // key (which sits right after the header) lands on an aligned
        // offset, and the next block (which sits right after the user
        // content) does too.
        const label total = headerBytes() + MemoryPool::alignUp(size);

        char* ptr = memory_.allocate(total);
        Entry* entry = reinterpret_cast<Entry*>(ptr);
        entry->keyLength = keyLen;
        entry->blockSize = total;
        entry->next      = buckets_[bucketIndex];
        // initialize key
        std::memcpy(entry->data(), keyData, size_t(keyLen));

        buckets_[bucketIndex] = label(ptr - memory_.start());

        nEntries_++;
        return entry;
    }

    void clear()
    {
        buckets_ = labelList(nBuckets_, -1);
        memory_.clear();
        nEntries_ = 0;
    }


    //- Member functions
    void* find(const string& key)
    {
        label bucketIndex;
        return lookup(key, bucketIndex);
    }

    //- find variant taking a precomputed FNV hash for the key.
    void* find(const string& key, uint hash)
    {
        label bucketIndex;
        return lookup(key, bucketIndex, hash);
    }

    //- Raw find: key bytes + length from a caller arena + precomputed hash.
    void* find(const char* keyData, label keyLen, uint hash)
    {
        label bucketIndex;
        return lookup(keyData, keyLen, bucketIndex, hash);
    }

    void* operator[](const string& key)
    {
        return find(key);
    }

    label nEntries() const
    {
        return nEntries_;
    }
    void bucketValues() const
    {
        forAll(buckets_, i)
        {
            Info << buckets_[i] << " ";
        }
        Info << endl;
    }


    // Block ingestion

        //- Clear this table and repopulate it from a packed buffer of
        //  blocks (e.g. an MPI All-to-all receive buffer). Memory is
        //  reused; nothing is released. Returns the number of blocks
        //  ingested (== \c nEntries on success).
        //
        //  Each source block already carries its own blockSize in its
        //  Entry header, so the walk needs no payload-type knowledge —
        //  rateEntry::decodeInto / TCEntry::decodeInto are thin wrappers
        //  that just forward here.
        label decodeFrom(const void* blocksStart, label nEntries)
        {
            clear();
            const char* src = static_cast<const char*>(blocksStart);
            for (label i = 0; i < nEntries; ++i)
            {
                auto* hdr = reinterpret_cast<const Entry*>(src);
                const label bs       = hdr->blockSize;
                const label userSize = bs - headerBytes();

                // Reconstruct the key as a string and let insertManual
                // allocate a fresh block of identical footprint in our
                // pool. Then copy the user content verbatim — the
                // contents are layout-compatible because the source
                // block was produced by the same Entry layout (same
                // headerBytes, same alignUp rules).
                const string key
                (
                    static_cast<const char*>(hdr->data()),
                    size_t(hdr->keyLength)
                );

                Entry* dst = static_cast<Entry*>(insertManual(key, userSize));
                std::memcpy(dst->data(), hdr->data(), size_t(userSize));

                src += bs;
            }
            return nEntries_;
        }


};



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
