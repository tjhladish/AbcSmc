#ifndef ABCSIM_H
#define ABCSIM_H

#include "AbcUtil.h"
#include <vector>

using std::vector;

typedef vector<ABC::float_type> (*SimFunc)(vector<ABC::float_type>, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par*);

#endif