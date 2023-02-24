
#include "AbcSmc.h"
#include "examples.h"
#include "dice.h"

int main(const int argc, const char* argv[]) {

    // convenience method for parsing arguments; from examples.h
    CLIArgs args = parse_args("abc_sql", argc, argv);

    AbcSmc* abc = new AbcSmc();
    abc->set_simulator(simulator);
    
    do_abc(abc, args, RNG);

    return 0;
}
