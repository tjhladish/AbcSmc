
#ifndef ABCSMC_PARRNG_H
#define ABCSMC_PARRNG_H

#include <vector>
#include <map>
#include <utility>
#include <cstddef>
#include <cassert>
#include <memory>

using std::vector;
using std::map;
using std::pair;

namespace ABC {

// ParRNG: a state-machine, `prng`, for uniform "sampling" of priors, posteriors, and pseudo parameters.
// Must be passed around by reference.
// When unlocked, the indexed parameters will increment (or wrapping back to 0, if at max index)
// their state. Incrementing a pseudo parameter will cause the prng to lock, meaning that future
// "draws" will not increment indexes, until unlocked.
//
// When using a mixture of priors, posteriors, and pseudo parameters, the repeated use of the
// prng is intended to
//   1. randomly sample from the priors
//   2. sequentially sample from the pseudo parameters
//   3. sequentially sample from the posteriors
//   4. ...however, to sample all combinations of pseudo parameters, before moving on to the next posterior sample
//
// So:
//
// Priors access and advance the RNG
// Posterior parameters use the posterior index. If prng is unlocked, either increment or reset the posterior index
// Pseudo parameters use the pseudo map to look up their index. If prng is unlocked:
//   1. unless index == maxindex, increment index in the map and lock the prng.
//   2. if index == maxindex, reset index to 0; do not lock the prng.
template<typename Par, typename RNG>
struct ParRNG {
    typedef std::shared_ptr<Par> ParPtr;
    ParRNG(
        RNG* rng, const vector<ParPtr> &mpars, const size_t posterior_size
    ) : _rng(rng), posteriorMaxIdx(posterior_size-1) {
        // this builds up the map of pseudo parameters => counters
        for (auto p : mpars) {                                 // for each parameter ...
            if ((not (p->isPosterior())) and (p->state_size() != 0)) { // if this isn't a posterior, and it has state size > 0 ...
                _pseudo[p] = { 0, p->state_size() - 1 };              // ... then register it 
            }               
        }
    }

    RNG* rng() const { return _rng; }
    void unlock() { lock = false; }

    size_t pseudo(const ParPtr p);
    size_t posterior();

    private:
        RNG* _rng;
        map<const ParPtr, pair<size_t, size_t>> _pseudo;
        bool lock = false;
        size_t posteriorIdx = 0;
        size_t posteriorMaxIdx;
};

template<typename Par, typename RNG>
size_t ParRNG<Par, RNG>::pseudo(const ParPtr p) {
    assert(_pseudo.count(p) == 1);
    auto ret = _pseudo[p];
    if (not lock) {
        if (ret.first < ret.second) { _pseudo[p].first++; lock = true; } else { _pseudo[p].first = 0; }
    }
    return ret.first;
}

template<typename Par, typename RNG>
size_t ParRNG<Par, RNG>::posterior() {
    auto ret = posteriorIdx;
    if (not lock) {
        if (posteriorIdx < posteriorMaxIdx) { posteriorIdx++; } else { posteriorIdx = 0; }
    }
    return ret;
}

} // namespace ABC

#endif // ABCSMC_PARRNG_H