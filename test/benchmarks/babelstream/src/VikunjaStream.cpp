// Copyright (c) 2015-16 Tom Deakin, Simon McIntosh-Smith,
// University of Bristol HPC
//
// For full license terms please see the LICENSE file distributed with this
// source code
//
// Cupla version created by Jeff Young in 2021
// Ported from cupla to alpaka by Bernhard Manfred Gruber in 2022

#include "VikunjaStream.h"

#include <vikunja/access/BlockStrategy.hpp>
#include <vikunja/operators/operators.hpp>
#include <vikunja/reduce/reduce.hpp>
#include <vikunja/transform/transform.hpp>
#include <vikunja/workdiv/BlockBasedWorkDiv.hpp>

#include <numeric>
#include <tuple>

constexpr auto TBSIZE = 1024;
constexpr auto DOT_NUM_BLOCKS = 256;

template<typename T>
VikunjaStream<T>::VikunjaStream(Idx arraySize, Idx deviceIndex)
    : arraySize(arraySize)
    , devHost(alpaka::getDevByIdx<DevHost>(0u))
    , devAcc(alpaka::getDevByIdx<Acc>(deviceIndex))
    , sums(alpaka::allocBuf<T, Idx>(devHost, DOT_NUM_BLOCKS))
    , d_a(alpaka::allocBuf<T, Idx>(devAcc, arraySize))
    , d_b(alpaka::allocBuf<T, Idx>(devAcc, arraySize))
    , d_c(alpaka::allocBuf<T, Idx>(devAcc, arraySize))
    , d_sum(alpaka::allocBuf<T, Idx>(devAcc, DOT_NUM_BLOCKS))
    , queue(devAcc)
{
    if(arraySize % TBSIZE != 0)
        throw std::runtime_error("Array size must be a multiple of " + std::to_string(TBSIZE));
    std::cout << "Using alpaka device " << alpaka::getName(devAcc) << std::endl;
}

template<typename T>
struct InitFunctor
{
private:
    T const init_value;

public:
    InitFunctor(T const init_value) : init_value(init_value)
    {
    }

    ALPAKA_FN_HOST_ACC T operator()(T const&) const
    {
        return init_value;
    }
};

template<typename T>
void VikunjaStream<T>::init_arrays(T initA, T initB, T initC)
{
    vikunja::transform::deviceTransform<Acc>(
        devAcc,
        queue,
        static_cast<Idx>(arraySize),
        alpaka::getPtrNative(d_a),
        alpaka::getPtrNative(d_a),
        InitFunctor{initA});
    vikunja::transform::deviceTransform<Acc>(
        devAcc,
        queue,
        static_cast<Idx>(arraySize),
        alpaka::getPtrNative(d_b),
        alpaka::getPtrNative(d_b),
        InitFunctor{initB});
    vikunja::transform::deviceTransform<Acc>(
        devAcc,
        queue,
        static_cast<Idx>(arraySize),
        alpaka::getPtrNative(d_c),
        alpaka::getPtrNative(d_c),
        InitFunctor{initC});

    alpaka::wait(queue);
}

template<typename T>
void VikunjaStream<T>::read_arrays(std::vector<T>& a, std::vector<T>& b, std::vector<T>& c)
{
    alpaka::memcpy(queue, alpaka::createView(devHost, a), d_a);
    alpaka::memcpy(queue, alpaka::createView(devHost, b), d_b);
    alpaka::memcpy(queue, alpaka::createView(devHost, c), d_c);
}

template<typename T>
void VikunjaStream<T>::copy()
{
    vikunja::transform::deviceTransform<Acc>(
        devAcc,
        queue,
        static_cast<Idx>(arraySize),
        alpaka::getPtrNative(d_a),
        alpaka::getPtrNative(d_c),
        [] ALPAKA_FN_HOST_ACC(T const data) { return data; });

    alpaka::wait(queue);
}


template<typename T>
void VikunjaStream<T>::mul()
{
    vikunja::transform::deviceTransform<Acc>(
        devAcc,
        queue,
        static_cast<Idx>(arraySize),
        alpaka::getPtrNative(d_c),
        alpaka::getPtrNative(d_b),
        [] ALPAKA_FN_HOST_ACC(T const data)
        {
            const T scalar = startScalar;
            return scalar * data;
        });

    alpaka::wait(queue);
}

template<typename T>
void VikunjaStream<T>::add()
{
    vikunja::transform::deviceTransform<Acc>(
        devAcc,
        queue,
        static_cast<Idx>(arraySize),
        alpaka::getPtrNative(d_a),
        alpaka::getPtrNative(d_b),
        alpaka::getPtrNative(d_c),
        [] ALPAKA_FN_HOST_ACC(T const a, T const b) { return a + b; });

    alpaka::wait(queue);
}

struct TriadKernel
{
    template<typename TAcc, typename T>
    ALPAKA_FN_ACC void operator()(TAcc const& acc, T* a, T const* b, T const* c) const
    {
        const T scalar = startScalar;
        auto const [i] = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc);
        a[i] = b[i] + scalar * c[i];
    }
};

template<typename T>
void VikunjaStream<T>::triad()
{
    vikunja::transform::deviceTransform<Acc>(
        devAcc,
        queue,
        static_cast<Idx>(arraySize),
        alpaka::getPtrNative(d_b),
        alpaka::getPtrNative(d_c),
        alpaka::getPtrNative(d_a),
        [] ALPAKA_FN_HOST_ACC(T const b, T const c)
        {
            const T scalar = startScalar;
            return b + scalar * c;
        });

    alpaka::wait(queue);
}

template<typename TType1, typename TType2>
class Zip_Iterator
{
public:
    using Type1 = std::decay_t<TType1>;
    using Type2 = std::decay_t<TType2>;

    using difference_type = Idx;
    using value_type = std::tuple<Type1, Type2>;
    using pointer = std::tuple<Type1, Type2>*;
    using reference = std::tuple<Type1, Type2>&;
    using iterator_category = std::random_access_iterator_tag;

private:
    Type1* m_type_ptr1;
    Type2* m_type_ptr2;

public:
    Zip_Iterator(TType1* ptr1, TType2* ptr2) : m_type_ptr1(ptr1), m_type_ptr2(ptr2)
    {
    }
    ALPAKA_FN_HOST_ACC auto operator[](const difference_type i) const
    {
        return std::make_tuple(m_type_ptr1[i], m_type_ptr2[i]);
    }
};


template<typename T>
void VikunjaStream<T>::nstream()
{
    Zip_Iterator zipped_input(alpaka::getPtrNative(d_b), alpaka::getPtrNative(d_c));

    vikunja::transform::deviceTransform<Acc>(
        devAcc,
        queue,
        static_cast<Idx>(arraySize),
        alpaka::getPtrNative(d_a),
        zipped_input,
        alpaka::getPtrNative(d_a),
        [] ALPAKA_FN_HOST_ACC(T const a, std::tuple<T, T> const bc)
        {
            auto b = std::get<0>(bc);
            auto c = std::get<1>(bc);
            const T scalar = startScalar;
            return a + b + scalar * c;
        });

    alpaka::wait(queue);
}

struct DotKernel
{
    template<typename TAcc, typename T>
    ALPAKA_FN_ACC void operator()(TAcc const& acc, T const* a, T const* b, T* sum, int arraySize) const
    {
        // TODO - test if sharedMem bug is affecting performance here
        auto& tb_sum = alpaka::declareSharedVar<T[TBSIZE], __COUNTER__>(acc);

        auto [i] = alpaka::getIdx<alpaka::Grid, alpaka::Threads>(acc);
        auto const [local_i] = alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc);
        auto const [totalThreads] = alpaka::getWorkDiv<alpaka::Grid, alpaka::Threads>(acc);

        tb_sum[local_i] = 0.0;
        for(; i < arraySize; i += totalThreads)
            tb_sum[local_i] += a[i] * b[i];

        auto const [blockDim] = alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc);
        for(int offset = blockDim / 2; offset > 0; offset /= 2)
        {
            alpaka::syncBlockThreads(acc);
            if(local_i < offset)
                tb_sum[local_i] += tb_sum[local_i + offset];
        }

        auto const [blockIdx] = alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc);
        if(local_i == 0)
            sum[blockIdx] = tb_sum[local_i];
    }
};

template<typename T>
T VikunjaStream<T>::dot()
{
#if 0
    Zip_Iterator zipped_input(alpaka::getPtrNative(d_a), alpaka::getPtrNative(d_b));

    auto transform_functor = [] ALPAKA_FN_HOST_ACC(T const i) -> T { return i; };
    auto reduce_functor = [] ALPAKA_FN_HOST_ACC(T const sum, std::tuple<T, T> const ab) -> T
    {
        auto a = std::get<0>(ab);
        auto b = std::get<1>(ab);
        return sum + a * b;
    };

    return vikunja::reduce::deviceTransformReduce<
        Acc,
        vikunja::workdiv::BlockBasedPolicy<Acc>,
        vikunja::MemAccess::MemAccessPolicy<Acc>,
        decltype(transform_functor),
        decltype(reduce_functor),
        decltype(zipped_input),
        decltype(devAcc),
        decltype(devHost),
        decltype(queue),
        Idx,
        vikunja::operators::UnaryOp<Acc, decltype(transform_functor), T>,
        vikunja::operators::BinaryOp<
            Acc,
            decltype(reduce_functor),
            T,
            typename std::iterator_traits<decltype(zipped_input)>::value_type>,
        T>(devAcc, devHost, queue, static_cast<Idx>(arraySize), zipped_input, transform_functor, reduce_functor);

#endif
    auto const workdiv = WorkDiv{DOT_NUM_BLOCKS, TBSIZE, 1};
    // auto const workdiv = alpaka::getValidWorkDiv(devAcc, DOT_NUM_BLOCKS * TBSIZE);
    alpaka::exec<Acc>(
        queue,
        workdiv,
        DotKernel{},
        alpaka::getPtrNative(d_a),
        alpaka::getPtrNative(d_b),
        alpaka::getPtrNative(d_sum),
        arraySize);
    alpaka::wait(queue);

    alpaka::memcpy(queue, sums, d_sum);
    T const* sumPtr = alpaka::getPtrNative(sums);
    // TODO(bgruber): replace by std::reduce, when gcc 9.3 is the baseline
    return std::accumulate(sumPtr, sumPtr + DOT_NUM_BLOCKS, static_cast<T>(0));
}

void listDevices()
{
    auto const count = alpaka::getDevCount<Acc>();
    std::cout << "Devices:" << std::endl;
    for(int i = 0; i < count; i++)
        std::cout << i << ": " << getDeviceName(i) << std::endl;
}

std::string getDeviceName(int deviceIndex)
{
    return alpaka::getName(alpaka::getDevByIdx<Acc>(deviceIndex));
}

std::string getDeviceDriver(int device)
{
    return "Not supported";
}

template class VikunjaStream<float>;
template class VikunjaStream<double>;
