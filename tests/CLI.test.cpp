#include <AbcSmc/CLI.h>
#include <string>

using namespace ABC;
using namespace std;

// Runs the ABC-SMC process, for some object that implements the ABCSMC *verbs*:
// - parse(a configuration file)
// - build(a verbosity level)
// - process(a random number generator, a verbosity level)
// - evaluate(a random number generator, a buffer size, a verbosity level)
struct MockABCSMC {
    MockABCSMC() {}
    void parse(const string &config_file, const size_t verbosity = 0) { 
        if (verbosity > 0) {
            cout << "parse(" << config_file << ")" << endl;
        }
    }
    void build(const size_t verbosity = 0) {
        if (verbosity > 0) {
            cout << "building, verbosity = " << verbosity << endl;
        }
    }
    void process(const gsl_rng* /* RNG */, const size_t verbosity = 0) { 
        if (verbosity > 0) {
            cout << "processing, verbosity = " << verbosity <<  endl;
        }
    }
    void evaluate(const optional<size_t> buffer_size, const size_t verbosity = 0) {
        if (verbosity > 0) {
            if (buffer_size.has_value()) {
                cout << "evaluating, buffer_size = " << buffer_size.value() << endl;
            } else {
                cout << "evaluating, buffer_size = ALL" << endl;
            }
        }
    }
};

int main() {
    gsl_rng * RNG = gsl_rng_alloc (gsl_rng_taus2);
    MockABCSMC* abcsmc = new MockABCSMC();

    const char* useargs[] = {"./CLI.test", "-h"};
    cerr << "Should print usage message:" << endl;
    auto args = ABC::parse_args(2, useargs);
    cerr << endl;

    cerr << "Should be building:" << endl;
    const char* bargs[] = {"./CLI.test", "config.json", "-v", "-b"};
    args = ABC::parse_args(4, bargs);
    ABC::run(abcsmc, args, RNG);
    cerr << endl;

    cerr << "Should be building + processing:" << endl;
    const char* sargs[] = {"./CLI.test", "config.json", "-v", "-s"};
    args = ABC::parse_args(4, sargs);
    ABC::run(abcsmc, args, RNG);
    cerr << endl;

    cerr << "Should be building + processing + evaluating:" << endl;
    const char* aargs[] = {"./CLI.test", "config.json", "-v", "-a"};
    args = ABC::parse_args(4, aargs);
    ABC::run(abcsmc, args, RNG);
    cerr << endl;

    cerr << "Should be processing:" << endl;
    const char* pargs[] = {"./CLI.test", "config.json", "-v", "-p"};
    args = ABC::parse_args(4, pargs);
    ABC::run(abcsmc, args, RNG);
    cerr << endl;

    cerr << "Should be evaluating:" << endl;
    const char* eargs[] = {"./CLI.test", "config.json", "-v", "-e"};
    args = ABC::parse_args(4, eargs);
    ABC::run(abcsmc, args, RNG);
    cerr << endl;

    cerr << "Should be evaluating + processing:" << endl;
    const char* cargs[] = {"./CLI.test", "config.json", "-v", "-c"};
    args = ABC::parse_args(4, cargs);
    ABC::run(abcsmc, args, RNG);
    cerr << endl;

    gsl_rng_free(RNG);
}