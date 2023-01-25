#include <dlfcn.h>
#include "AbcSmc.h"

void usage() {
    std::cerr << std::endl << 
    "\tUsage: ./abc_dynamic config.json simulator.so" <<
    std::endl << std::endl;
}


int main(int argc, char* argv[]) {

    if (not (argc == 3 or argc == 5) ) {
        usage();
        exit(100);
    }

    int buffer_size = 1;

    if ( argc == 5 ) {
        buffer_size = atoi(argv[4]);
    }

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));
    abc->set_simulation(new AbcFPtr(argv[2]));
    abc->simulate_next_particles(buffer_size);

    return 0;
}
