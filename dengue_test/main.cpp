#include "AbcSmc.h"

using namespace std;

int main(int argc, char* argv[]) {

    if (argc != 2) {
        cerr << "\n\tUsage: ./abc abc_config_file.json\n\n";
        return 100;
    }
    
    const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);
    gsl_rng_set(RNG, time (NULL) * getpid()); // seed the rng using sys time and the process id

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));
    abc->run(RNG);

    return 0;
}
