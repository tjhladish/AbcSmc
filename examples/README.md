# AbcSmc Examples

This directory contains examples for the different use cases / modes for the `AbcSmc` library.

# Example Simulator

Those examples are all motivated in terms of a particular core simulation: rolling $n$ dice, each with $m$ faces (labelled 1 to $m$), and observing the sum and standard deviation of the faces for a particular roll.

This problem is a convolution of identical, indepently drawn, discrete uniform deviates. Per [Caiado and Rathie](https://www.researchgate.net/profile/Camila-Caiado/publication/228457326_Polynomial_coefficients_and_distribution_of_the_sum_of_discrete_uniform_variables/links/00b7d53501e1ea6fbc000000/Polynomial-coefficients-and-distribution-of-the-sum-of-discrete-uniform-variables.pdf): the theoretical mean of this distribution is $\mu = \frac{n(m+1)}{2}$. The variance in terms of a particular roll (not the variance in sums of many rolls), is $\sigma^2 = \frac{m^2-1}{12}$. We can invert these equations to solve for $m$ and $n$ in terms of $\mu$ and $\sigma$:

$$
m = \sqrt{12\sigma^2 + 1} \\
n = \frac{2\mu}{m+1}
$$

# Technical Note

This is an plain language description of technical code in this directory. The relationships between various source elements, executables, and outputs are also formally captured in the `Makefile` in this directory.

# Use Cases

`AbcSmc` supports two stages of the typical dynamical modeling analysis workflow: fitting and projection (including scenario analysis). In fitting mode, given a model, parameters, and observed metrics, `AbcSmc` will estimate the parameters that best reproduce the observed measures with that model. In projection mode, `AbcSmc` will use posterior parameter samples and scenario input grids with the model to compute projected outcomes.

In general, these two stages will happen together using the "same" model. Emphasis on "same" here to indicate that, for example, the fitting version of the model might be used to create one set of outputs, while in projection a different set might be desired for the research question, leading to different instantiations of a simulator function, but around the same core dynamical model.

In addition to the different modes, `AbcSmc` supports a variety of ways to work with a concrete simulator function.

## Options to Specify Simulator

`AbcSmc` supports three different ways to use a simulator:

 1. The simulator is compiled into an executable along with the library (`main_sql.cpp` --> `abc_sql`). (fastest running option, requires more recompilation)
 2. ... a shared library (`dice.cpp` --> `dice.so`) and loaded dynamically (`main_dynamic.cpp` --> `abc_dynamic`). (next fastest, provides more simulator flexibility)
 3. ... a separate executable (`dice_game.cpp`) and invoked via forking and interacting with that executable via stdout and stdin (`main_exec.cpp` --> `abc_exec`). (slowest, can provide more simulator flexibility)

This `Makefile` demonstrates how to compile and run each of these options.