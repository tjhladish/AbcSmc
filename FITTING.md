# Notes on FITTING mode flow

Setup => Simulate => Process [=> Simulation => Process [...]] - ends with 

## Setup

### Current

### Desired

## Simulate

### Current

### Desired

## Process

### Current

read_SMC_sets_from_database - receives DB + reference to vector-of-vector of serials (will become: smcsets of serials)

side effects: fills _particle_parameters, _particle_metrics - vectors of Mat2Ds, which get fed into PLS step

if any *done* smcset is complete, but nothing is "in" the posterior (i.e. all -1) => perform PLS => write the rankings (_filter_particles) to that SMC set (note: presumably this is only the "last" smcset, but current implementation will do whatever it gets)

_filter_particles - receives metrics, parameters, next_set_size (and smcset id).
side effects: fills _predictive_priors (by smcset id, via _fp_helper), writes to database (via _fp_helper)

_fp_helper - receives smcset id, distances, size
side effects: 

### Desired

desired behavior: split into stages:

read_SMC_complete (bool for last vs all)
 - read all (last only) SMC sets with any members assigned to posterior - side effect: fill particle_parameters, particle_metrics
 - confirm all sets => if not, warn => exit; if yes, read into _parameters, _metrics

rank_SMC_last
 - read last SMC set with any members done
 - confirm that set is complete => if not, warn => exit
 - check if already ranked => if yes, warn (overwrite option?) => exit
 - write to storage (maybe bool controlled?)

summarize_SMC
 - output summary stats for whatever is in _parameters, _metrics?

### Backwards Compatibility

 - would need to write AbcSQLite methods in terms of an arbitrary Db reference (might be fine? problem for abcstorage wanting to be aware of params metrics / etc - but good for testing those methods without having to construct a full object, and then the "convenience" versions on generic AbcStorage can have their overrides be written in terms of private members + those static methods) 

then read_SMC_sets_from_database is:
 - read_SMC_complete => rank_SMC_last => summarize_SMC
