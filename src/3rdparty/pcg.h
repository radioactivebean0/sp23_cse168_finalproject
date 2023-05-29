#pragma once
#ifndef pcg_h
#define pcg_h
#include <cstdint>

// Lajolla uses a random number generator called PCG https://www.pcg-random.org/
// which is a very lightweight random number generator based on a simple
// postprocessing of a standard linear congruent generators.
// PCG has generally good statistical properties and is much cheaper to compute compared
// to more advanced RNG e.g., Merserner Twister.
// Highly recommend Melissa O'Neill's talk about PCG https://www.youtube.com/watch?v=45Oet5qjlms

// A crucial feature of PCG is that it allows multiple "streams": 
// given a seed, we can initialize many different streams of RNGs 
// that have different independent random numbers.

struct pcg32_state {
    uint64_t state;
    uint64_t inc;
};

// http://www.pcg-random.org/download.html
inline uint32_t next_pcg32(pcg32_state &rng) {
    uint64_t oldstate = rng.state;
    // Advance internal state
    rng.state = oldstate * 6364136223846793005ULL + (rng.inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = uint32_t(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = uint32_t(oldstate >> 59u);
    return uint32_t((xorshifted >> rot) | (xorshifted << ((-rot) & 31)));
}

// https://github.com/wjakob/pcg32/blob/master/pcg32.h#L47
inline pcg32_state init_pcg32(uint64_t stream_id = 1, uint64_t seed = 0x853c49e6748fea9bULL) {
    pcg32_state s;
    s.state = 0U;
    s.inc = (stream_id << 1u) | 1u;
    next_pcg32(s);
    s.state += seed;
    next_pcg32(s);
    return s;
}

template <typename T>
T next_pcg32_real(pcg32_state &rng) {
    return T(0);
}

// https://github.com/wjakob/pcg32/blob/master/pcg32.h
template <>
inline float next_pcg32_real(pcg32_state &rng) {
    union {
        uint32_t u;
        float f;
    } x;
    x.u = (next_pcg32(rng) >> 9) | 0x3f800000u;
    return x.f - 1.0f;
}

// https://github.com/wjakob/pcg32/blob/master/pcg32.h
template <>
inline double next_pcg32_real(pcg32_state &rng) {
    union {
        uint64_t u;
        double d;
    } x;
    x.u = ((uint64_t) next_pcg32(rng) << 20) | 0x3ff0000000000000ULL;
    return x.d - 1.0;
}

/// Generate a uniformly distributed unsigned 32-bit random number
inline uint32_t nextUInt(pcg32_state &rng) {
    uint64_t oldstate = rng.state;
    rng.state = oldstate * 6364136223846793005ULL + (rng.inc|1);
    uint32_t xorshifted = (uint32_t) (((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = (uint32_t) (oldstate >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
}

/// Generate a uniformly distributed number, r, where 0 <= r < bound
inline uint32_t nextUInt(pcg32_state &rng, uint32_t bound) {
    // To avoid bias, we need to make the range of the RNG a multiple of
    // bound, which we do by dropping output less than a threshold.
    // A naive scheme to calculate the threshold would be to do
    //
    //     uint32_t threshold = 0x100000000ull % bound;
    //
    // but 64-bit div/mod is slower than 32-bit div/mod (especially on
    // 32-bit platforms).  In essence, we do
    //
    //     uint32_t threshold = (0x100000000ull-bound) % bound;
    //
    // because this version will calculate the same modulus, but the LHS
    // value is less than 2^32.

    uint32_t threshold = (~bound+1u) % bound;

    // Uniformity guarantees that this loop will terminate.  In practice, it
    // should usually terminate quickly; on average (assuming all bounds are
    // equally likely), 82.25% of the time, we can expect it to require just
    // one iteration.  In the worst case, someone passes a bound of 2^31 + 1
    // (i.e., 2147483649), which invalidates almost 50% of the range.  In
    // practice, bounds are typically small and only a tiny amount of the range
    // is eliminated.
    for (;;) {
        uint32_t r = nextUInt(rng);
        if (r >= threshold)
            return r % bound;
    }
}
#endif