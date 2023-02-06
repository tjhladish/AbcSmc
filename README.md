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

## Usage

### JSON Configuration File

`AbcSmc` uses a JSON-formatted file to specify run configuration. `AbcSmc` is essentially used in two-modes: fitting and projecting. The typical research use case uses both modes with the same simulation: first fitting parameter distributions, then using the fitted values for projections (often under varying scenarios).

#### Top Level Keys

TODO list required / optional / computed top level keys

#### Parameters / Metrics

The simulation also fundamental revolves around "parameters" and "metrics". Relative to the model and usage of `AbcSmc` in projection mode, these are essentially the inputs and outputs (respectively) - "parameters" are provided to the simulator, it runs, and then returns "metrics". However, in fitting mode, the "parameters" can also be thought of the outputs - they are the posterior distributions that `AbcSmc` is attempting to find - while the "metrics" are inputs, the targets for fitting to reproduce.

In addition to structural schema, differences in usage mode put different constraints on "parameters". Both "parameters" and "metrics" should be a JSON array of objects.

The parameter objects have the following JSON schema structure:

```JSON
{
    "type" : "object",
    "properties" : {
        "name" : {
            "type" : "string"
        },
        "short_name" : {
            "type" : "string"
        },
        "num_type" : {
            "type" : "string",
            "enum" : ["UNIFORM", "NORMAL", "PSEUDO", "POSTERIOR"]
        },
        "num_type" : {
            "type" : "string",
            "enum" : ["INT", "FLOAT"],
            "default" : "FLOAT"
        },
        "par1" : {
            "type" : ["integer", "number"]
        },
        "par2": {
            "type" : ["integer", "number"]
        },
        "step": {
            "type" : ["integer", "number"],
            "default" : 1
        },
        "values": {
            "type" : "array",
            "items" : ["integer", "number"]
        },
        "ranks" : {
            "type" : "integer"
        },
        "description" : {
            "type" : "string"
        }
    },
    "required" : [ "name", "dist_type" ]
}
```

TODO: formally constrain the combinations of "par" key types for specification.

If a prior-type parameter (i.e. "dist_type" is "UNFORM" or "NORMAL" (or "GAUSSIAN")): must specify "par1" and "par2". Future: allow these to instead be distribution-specific parameter names (e.g. "min" and "max" for "UNIFORM"; "mean" and "sd" for "NORMAL").

If a posterior parameter, can instead use "ranks" - this is "how many" posterior particles to use. It will effectively be translated to a sequence from 0 to ranks-1. Generically, any specific sample ids from an empirical posterior can be used (and therefore, specified in the configuration file), but practically you'll most frequently be using the 0 through N-1 (for a total of N), so "ranks" provides an easier way to specify that.

TODO: enable setting "ranks" at the top level, since it must be identical across all POSTERIOR parameters?

For PSEUDO parameters, you can specify a `par1` (start) / `par2` (end) / `step` (optional, defaults to 1) OR a `values` array. If provided start/end/step, the `values` array is effectively created from that, by creating a `par1, par1 + step, par1 + 2*step,  ..., par1 + n*step` array, where the last element is the only value equal-to-or-greater-than `par2`.

The `par1` must be less than or equal to `par2`, and `step` must be positive. If `par1 == par2`, the parameter sequence has only one value. This is considered a special case, but does have practical uses, e.g. when doing a batch of test runs using a subset of scenario settings.

That expression mirrors how the values are created (versus e.g. repeated `+=` operations); if that presents floating-point arithmetic issues, you may need to manually specify values.

The metric objects have the following JSON schema structure:

```JSON
{
    "type" : "object",
    "properties" : {
        "name" : {
            "type" : "string"
        },
        "value" : {
            "type" : ["integer", "number"]
        },
        "short_name" : {
            "type" : "string"
        },
        "num_type" : {
            "type" : "string",
            "enum" : ["INT", "FLOAT"],
            "default" : "FLOAT"
        },
        "description" : {
            "type" : "string"
        }
    },
    "required" : [ "name", "value" ]
}
```