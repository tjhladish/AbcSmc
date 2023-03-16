#include <AbcSmc/ParRNG.h>
#include <gsl/gsl_rng.h>
#include <vector>
#include <iostream> // cout, endl

using namespace ABC;
using namespace std;

struct MockPar;

typedef ParRNG<MockPar, gsl_rng> PRNG;

struct MockPar {
    MockPar(const size_t max_index) : mi(max_index) {}

    const size_t mi;
    virtual bool isPosterior() const = 0;
    virtual double sample(PRNG& rng) const = 0;
    size_t max_index() const { return mi; }
};

struct MockPrior : MockPar {
    MockPrior() : MockPar(0) {}
    double sample(PRNG& rng) const override { return gsl_rng_uniform(rng.rng()); }
    bool isPosterior() const override { return false; }
};

struct MockPseudo : MockPar {
    MockPseudo(const std::vector<double> & vals) : MockPar(vals.size() - 1), vs(vals) {}
    bool isPosterior() const override { return false; }
    double sample(PRNG& rng) const override { return vs[rng.pseudo(this)]; }
    private:
        const std::vector<double> vs;
};

struct MockPost : MockPar {
    MockPost(const size_t max_rank) : MockPar(max_rank) {}
    bool isPosterior() const override { return true; }
    double sample(PRNG& rng) const override { return static_cast<double>(rng.posterior()); }
};

int main() {
    gsl_rng * RNG = gsl_rng_alloc (gsl_rng_taus2);
    const size_t postsize = 10;
    vector<MockPar*> mps = { new MockPrior(), new MockPseudo({1, 2, 3}), new MockPseudo({5, 4, 3, 2, 1}), new MockPost(postsize) };
    auto prng = ParRNG(RNG, mps, postsize);
    for (size_t i = 0; i < postsize*3*5; ++i) {
        prng.unlock();
        for (size_t parIdx = 0; parIdx < mps.size(); ++parIdx) {
            if (not mps[parIdx]->isPosterior()) {
                if (mps[parIdx]->max_index() == 0) {
                    cout << parIdx << " (prior): " << mps[parIdx]->sample(prng) << endl;
                } else {
                    cout << parIdx << " (pseudo): " << mps[parIdx]->sample(prng) << endl;
                }
            }
        }
        for (size_t parIdx = 0; parIdx < mps.size(); ++parIdx) {
            if (mps[parIdx]->isPosterior()) {
                cout << parIdx << " (posterior): " << mps[parIdx]->sample(prng) << endl;
            }
        }
    }

    gsl_rng_free(RNG);
}