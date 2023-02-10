# AbcSmc
*Sequential Monte Carlo Approximate Bayesian Computation with Partial Least Squares*

`AbcSmc` is a combined parameter estimation and scenario projection library implemented in C++. `AbcSmc` specifically focuses on:

 - complex, computationally demanding, stochastic models,
 - fitted to disparate types of empirical data, which defy convenient likelihoods,
 - and ambiguity as to what parameters are key and what metrics are informative.
 
 `AbcSmc` use partial least squares to address problems arising from parameters and/or empirical metrics that co-vary or are unidentifiable (parameters) or uninformative (metrics). `AbcSmc` is intended to support use in high performance computation (e.g. cluster or supercomputer) environments, consistent with its focus on models with longer run times and large memory footprints. `AbcSmc` includes a convenient means of distributing and gathering work in HPC environments: the program pulls jobs from and writes output to a standardized SQL database, and implements a dynamic load balancing scheme to compensate for variable simulation run times and hardware failures. `AbcSmc` uses SQLite for the database, which has a host of library wrappers in many programming languages, making post processing analysis smoother, and simplifying data management as SQLite databases are single, self-contained, portable files.

## Peer-reviewed publications that have used AbcSmc:

Hladish, T.J., C.A.B. Pearson, D.P. Rojas, K.B. Toh, P. Manrique-Saide, G.M. Vazquez-Prokopec, M.E. Halloran, I.M. Longini (2020) Designing effective control of dengue with combined interventions. Proc Natl Acad Sci U S A.  2020 Jan 23. pii: 201903496. [doi:10.1073/pnas.1903496117](https://doi.org/10.1073/pnas.1903496117)

Hladish, T.J., C.A.B. Pearson, D.P. Rojas, H. Gómez-Dantés, M.E. Halloran, G.M. Vazquez-Prokopec, I.M. Longini (2018) Forecasting the effectiveness of indoor residual spraying for reducing dengue burden.  PLoS Negl Trop Dis 12(6): e0006570. [doi:10.1371/journal.pntd.0006570](https://doi.org/10.1371/journal.pntd.0006570)

Flasche, S., M. Jit, I. Rodríguez-Barraquer, L. Coudeville, M. Recker, K. Koelle, G. Milne, T.J. Hladish, A. Perkins, D.A.T. Cummings, I. Dorigatti, D.J. Laydon, G. España, J. Kelso, I. Longini, J. Lourenco, C.A.B. Pearson, R.C. Reiner, L. Mier-y-Terán-Romero, K. Vannice, N. Ferguson (2016) The long term safety, public health impact, and cost effectiveness of routine vaccination with a recombinant, live-attenuated dengue vaccine (Dengvaxia): a model comparison study. PLoS Medicine 13(11): e1002181. [doi:10.1371/journal.pmed.1002181](https://doi.org/10.1371/journal.pmed.1002181)

Hladish, T.J., C.A.B. Pearson, D.L. Chao, D.P. Rojas, G.L. Recchia, H. Gómez-Dantés, M.E. Halloran, J.R.C. Pulliam, I.M. Longini (2016) Projected impact of dengue vaccination in Yucatán, Mexico. PLoS Negl Trop Dis 10(5): e0004661. [doi:10.1371/journal.pntd.0004661](https://doi.org/10.1371/journal.pntd.0004661)

## License
AbcSmc is licensed under GPLv3 (or later, at your discretion); see license [here](LICENSE).

## Database Schema

`AbcSmc` divides its database according to its different modes and tasks. The primary identifier representing a calculation is the `serial`, which broadly links across all of the main tables. The database is also prepopulated with [views](https://www.sqlite.org/lang_createview.html): these views support backwards compatibility and direct human-readability of the database. These views manage materialization of underlying representations, and provide a convenient way to access "wide" perspective data while maintaining the benefits of "long" format data. In general, the views have human-meaningful names, and the underlying tables have algorithmically generated, shorthand names.

### Job Management

For all uses, `AbcSmc` manages tracking of work with a `rjob` table and a corresponding `job` view. The `job` view concerns:

 - `serial` associated with a calculation.
 - `smcSet`: indicating which round of sequential monte carlo the serial is associated with. This will be `NULL` for projection mode.
 - `particleIdx`: which underlying "particle" the job concerns. This is `NULL` for fitting mode.
 - `startTime` and `duration`: time-keeping records for tracking performance characteristics of the simulator (integer seconds; since 1970-01-01 and since `startTime`, respectively)
 - `status`: the current job status, either `Q`ueued, `P`aused, or `R`unning. 

The `rjob` table does not include the additional particle meta data (`smcSet`, `particleIdx`, `posterior`) and encodes run status with a number based on `rjob_status` table.

TODO: suppress view columns based on mode?

### Settings

In general, simulations are expressed in terms of variable values. Some of these are fit by the `AbcSmc` process, others may be single values that need to be set (but might vary more broadly in the use of the model), and others might vary systematically as part of e.g. a scenario analysis.

The `rset_*` tables manage the latter inputs: those which are *not* fitted, but rather have single or systematically varying values. The `setting` view materializes this table into a view of complete "setting" combination for use a simulation. The `setting` view consists of:

 - `id`: the unique id defining this combination of settings
 - `...`: assorted columns corresponding to the settings used in the simulation. Their values correspond to the values provided to simulator for runs corresponding to this setting id. Any settings that aren't pertinent to a particular id are `NULL`.

TODO: future normalization approach:

### Priors

### Posteriors

The `rset_*` tables consist of `rset_ref` (reference values, generally for non-varying inputs), `rset_...` (a table each for the varying parameters), and `rset_combos` (having `id`, `ref_id`, and then `..._id` columns, which can then be joined on `rset_ref` and all of the `rset_...` tables to create `setting` view).

The SMC-ABC-PLS fitting process


https://stackoverflow.com/questions/11563869/update-multiple-rows-with-different-values-in-a-single-sql-query

https://stackoverflow.com/questions/2444708/sqlite-long-to-wide-formats

https://www.sqlite.org/lang_createview.html