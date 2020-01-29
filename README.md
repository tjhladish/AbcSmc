# AbcSmc
Sequential Monte Carlo Approximate Bayesian Computation with Partial Least Squares

AbcSmc is a parameter estimation library implemented in C++ that has been developed to enable fitting complex stochastic models to disparate types of empirical data.  We use partial least squares to address problems arising from parameters and/or empirical metrics that co-vary or are unidentifiable (parameters) or uninformative (metrics).  Because of the long running times, often requiring many processor-core years of computation, AbcSmc is particularly well-suited to being used in high performance (e.g. cluster or supercomputer) environments. AbcSmc includes a convenient means of distributing and gathering work in HPC environments: the program pulls jobs from and writes output to a standardized SQL database, and implements a dynamic load balancing scheme to compensate for variable simulation run times and hardware failures. AbcSmc uses SQLite for the database, for portability of data.

## Peer-reviewed publications that have used AbcSmc:

Hladish, T.J., C.A.B. Pearson, D.P. Rojas, K.B. Toh, P. Manrique-Saide, G.M. Vazquez-Prokopec, M.E. Halloran, I.M. Longini (2020) Designing effective control of dengue with combined interventions. Proc Natl Acad Sci U S A.  2020 Jan 23. pii: 201903496. [doi:10.1073/pnas.1903496117](https://doi.org/10.1073/pnas.1903496117)

Hladish, T.J., C.A.B. Pearson, D.P. Rojas, H. Gómez-Dantés, M.E. Halloran, G.M. Vazquez-Prokopec, I.M. Longini (2018) Forecasting the effectiveness of indoor residual spraying for reducing dengue burden.  PLoS Negl Trop Dis 12(6): e0006570. [doi:10.1371/journal.pntd.0006570](https://doi.org/10.1371/journal.pntd.0006570)

Flasche, S., M. Jit, I. Rodríguez-Barraquer, L. Coudeville, M. Recker, K. Koelle, G. Milne, T.J. Hladish, A. Perkins, D.A.T. Cummings, I. Dorigatti, D.J. Laydon, G. España, J. Kelso, I. Longini, J. Lourenco, C.A.B. Pearson, R.C. Reiner, L. Mier-y-Terán-Romero, K. Vannice, N. Ferguson (2016) The long term safety, public health impact, and cost effectiveness of routine vaccination with a recombinant, live-attenuated dengue vaccine (Dengvaxia): a model comparison study. PLoS Medicine 13(11): e1002181. [doi:10.1371/journal.pmed.1002181](https://doi.org/10.1371/journal.pmed.1002181)

Hladish, T.J., C.A.B. Pearson, D.L. Chao, D.P. Rojas, G.L. Recchia, H. Gómez-Dantés, M.E. Halloran, J.R.C. Pulliam, I.M. Longini (2016) Projected impact of dengue vaccination in Yucatán, Mexico. PLoS Negl Trop Dis 10(5): e0004661. [doi:10.1371/journal.pntd.0004661](https://doi.org/10.1371/journal.pntd.0004661)

## License
AbcSmc is licensed under GPLv3 (or later, at your discretion); see license [here](LICENSE).
