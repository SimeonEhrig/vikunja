//
// Created by mewes30 on 19.12.18.
//

#pragma once

#include <vikunja/mem/iterator/PolicyBasedBlockIterator.hpp>
#include <alpaka/alpaka.hpp>

namespace vikunja {
namespace reduce {
namespace detail {


    struct LinearMemAccessPolicy {
        template<typename TAcc, typename TIdx>
        ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE
        static auto getStartIndex(TAcc const &acc, TIdx const &problemSize, TIdx const &blockSize) -> TIdx const {
            auto gridDimension{alpaka::workdiv::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0]};
            auto indexInBlock{alpaka::idx::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0]};
            auto gridSize = gridDimension * blockSize;
            // TODO: catch overflow
            return (problemSize * indexInBlock) / gridSize;
        }

        template<typename TAcc, typename TIdx>
        ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE
        static auto getEndIndex(TAcc const &acc, TIdx const &problemSize, TIdx const &blockSize) -> TIdx const {
            auto gridDimension{alpaka::workdiv::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0]};
            auto indexInBlock{alpaka::idx::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0]};
            auto gridSize = gridDimension * blockSize;
            // TODO: catch overflow
            return (problemSize * indexInBlock + problemSize) / gridSize;
        }

        template<typename TAcc, typename TIdx>
        ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE
        static auto getStepSize(TAcc const &acc, TIdx const &problemSize, TIdx const &blockSize) {
            return 1;
        }
    };

    struct GridStridingMemAccessPolicy {
        template<typename TAcc, typename TIdx>
        ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE
        static auto getStartIndex(TAcc const &acc, TIdx const &problemSize, TIdx const &blockSize) -> TIdx const {
            return alpaka::idx::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0];
        }

        template<typename TAcc, typename TIdx>
        ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE
        static auto getEndIndex(TAcc const &acc, TIdx const &problemSize, TIdx const &blockSize) -> TIdx const {
            return problemSize;
        }

        template<typename TAcc, typename TIdx>
        ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE
        static auto getStepSize(TAcc const &acc, TIdx const &problemSize, TIdx const &blockSize) {
            auto gridDimension{alpaka::workdiv::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0]};
            return gridDimension * blockSize;
        }
    };


    template<typename TRed, uint64_t size>
    struct sharedStaticArray
    {
        TRed data[size];

        ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE TRed &operator[](uint64_t index) {
            return data[index];
        }
        ALPAKA_FN_HOST_ACC ALPAKA_FN_INLINE const TRed &operator[](uint64_t index) const {
            return data[index];
        }
    };

    template<uint64_t TBlockSize, typename TMemAccessPolicy, typename TRed, typename TFunc>
    struct BlockThreadReduceKernel {
        template<typename TAcc, typename TIdx,
                typename TInputIterator, typename TOutputIterator>
        ALPAKA_FN_ACC void operator()(TAcc const &acc,
                TInputIterator const * const source,
                TOutputIterator destination,
                TIdx const &n,
                TFunc func) const  {
            // Shared Mem
            auto &sdata(
                    alpaka::block::shared::st::allocVar<sharedStaticArray<TRed, TBlockSize>,
                    __COUNTER__>(acc));


            // blockIdx.x
            auto blockIndex{alpaka::idx::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0]};
            // threadIdx.x
            auto threadIndex{alpaka::idx::getIdx<alpaka::Block, alpaka::Threads>(acc)[0]};
            // gridDim.x
            auto gridDimension{alpaka::workdiv::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0]};

            // blockIdx.x * TBlocksize + threadIdx.x
            auto indexInBlock{alpaka::idx::getIdx<alpaka::Grid, alpaka::Threads>(acc)[0]};

            using MemPolicy = TMemAccessPolicy;
            vikunja::mem::iterator::PolicyBasedBlockIterator<MemPolicy, TAcc, TInputIterator> iter{source, acc, n, TBlockSize};
            // grid striding
            /*auto startIndex = indexInBlock;
            auto endIndex = n;
            auto stepSize = gridSize;*/

            // sequential
            auto startIndex = MemPolicy::getStartIndex(acc, n, TBlockSize);
            auto endIndex = MemPolicy::getEndIndex(acc, n, TBlockSize);
            auto stepSize = MemPolicy::getStepSize(acc, n, TBlockSize);
            //std::cout << "startIndex: " << startIndex << ", endIndex: " << endIndex << "\n";
           /* if(iter >= iter.end()) {
                return;
            }
            auto tSum = *iter;
            ++iter;
            while(iter < iter.end()) {
                tSum = func(tSum, *iter);
                ++iter;
            }*/

            auto i = startIndex;
            if(i >= n) {
                return;
            }

            auto tSum = *(source + i);
            i += stepSize;
            // Level 1: Grid reduce, reading from global memory
            while(i < endIndex) {
                tSum = func(tSum, *(source + i));
                i += stepSize;
            }
            if(threadIndex < n) {
                sdata[threadIndex] = tSum;
            }

            alpaka::block::sync::syncBlockThreads(acc);

            // blockReduce
            // unroll for better performance
            for(TIdx bs = TBlockSize, bSup = (TBlockSize + 1) / 2;
            bs > 1; bs = bs / 2, bSup = (bs + 1) / 2) {
                bool condition = threadIndex < bSup && // only first half of block is working
                         (threadIndex + bSup) < TBlockSize && // index for second half must be in bounds
                         (indexInBlock + bSup) < n; // if element in second half has ben initialized before
                if(condition) {
                    sdata[threadIndex] = func(sdata[threadIndex], sdata[threadIndex + bSup]);
                }
                alpaka::block::sync::syncBlockThreads(acc);
            }
            if(threadIndex == 0) {
                *(destination + blockIndex) = sdata[0];
            }
        }
    };
}
}
}