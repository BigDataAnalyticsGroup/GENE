#pragma once

#include <utility>
#include <functional>
#include <cstdint>

// Define hash function for pairs based on boost::hash
template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

struct pair_hash
{
    template<typename S, typename T>
    inline std::size_t operator()(const std::pair<S, T> & v) const
    {
         std::size_t seed = 0;
         hash_combine<S>(seed, v.first);
         hash_combine<T>(seed, v.second);
         return seed;
    }
};

struct Murmur3
{
    uint32_t operator()(uint32_t v) const {
        v ^= v >> 16;
        v *= 0x85ebca6b;
        v ^= v >> 13;
        v *= 0xc2b2ae35;
        v ^= v >> 16;
        return v;
    }

    uint64_t operator()(uint64_t v) const {
        v ^= v >> 33;
        v *= 0xff51afd7ed558ccd;
        v ^= v >> 33;
        v *= 0xc4ceb9fe1a85ec53;
        v ^= v >> 33;
        return v;
    }
};

// Knuth multiplicative hash
// https://stackoverflow.com/questions/11871245/knuth-multiplicative-hash
struct Multiplicative_hash
{
    uint32_t operator()(uint32_t v) const {
        return v * UINT32_C(2654435761);
    }
};
