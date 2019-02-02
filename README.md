# Potential functions for first-order methods

This code comes jointly with the following reference:

> [1] Adrien Taylor, and Francis Bach. "Stochastic first-order methods: non-asymptotic and computer-aided analyses via potential functions."

Date:    February 4, 2019

#### Authors

- [**Adrien Taylor**](http://www.di.ens.fr/~ataylor/)
- [**Francis Bach**](https://www.di.ens.fr/~fbach/)

**Note:** This code requires [YALMIP](https://yalmip.github.io/) along with a suitable SDP solver (e.g., Sedumi, SDPT3, Mosek).


## Organization of the code

### Deterministic first-order methods
- [`A_GradientDescent`](1_Deterministic_SmoothConvex/A_GradientDescent.m) Code for reproducing the result of Section 3.2.1, Theorem 3 and Figure 1 (gradient method).
- [`B_ProximalGradientDescent`](1_Deterministic_SmoothConvex/B_ProximalGradientDescent.m) Code for reproducing the result presented in Appendix C.2, Theorem 9 and Figure 1-like results (proximal gradient method).
- [`C_StepsizeSelection_FirstAcceleratedMethod`](1_Deterministic_SmoothConvex/C_StepsizeSelection_FirstAcceleratedMethod.m) Code for reproducing the first result of Appendix C.3, Theorem 10 and Figure 3 (first accelerated gradient method).
- [`D_StepsizeSelection_SecondAcceleratedMethod`](1_Deterministic_SmoothConvex/D_StepsizeSelection_SecondAcceleratedMethod.m) Code for reproducing the second result of Appendix C.3, Theorem 11 and Figure 4 (second accelerated gradient method).

### Stochastic first-order methods: unbiased oracles with bounded variance
- [`A_StochasticGradientDescent`](2_Stochastic_BoundedVariance/A_StochasticGradientDescent.m) Code for reproducing the first result of Section 3.2.2, Theorem 5 and Figure 2 (stochastic gradient descent).
- [`B_StochasticGradientDescentWithAveraging`](2_Stochastic_BoundedVariance/B_StochasticGradientDescentWithAveraging.m) Code for reproducing the second result of Section 3.2.2, Theorem 6  (stochastic gradient descent with averaging).
- [`C_StochasticGradientDescentWithPrimalAveraging`](2_Stochastic_BoundedVariance/C_StochasticGradientDescentWithPrimalAveraging.m) Code for reproducing the third result of Section 3.2.2, Theorem 7 (stochastic gradient descent with primal averaging).
- [`D_StochasticGradientDescent_EvaluationAtAveragedPoint`](2_Stochastic_BoundedVariance/D_StochasticGradientDescent_EvaluationAtAveragedPoint.m) Code for reproducing the result of Appendix D.4, Theorem 12 (stochastic gradient descent with evaluation at the averaged iterate).

### Stochastic first-order methods: unbiased oracles arising from sampling in expectations of smooth convex functions

#### Over-parametrized models
- [`A_ParameterSelection`](3_Stochastic_Overparametrized/A_ParameterSelection.m) Code for reproducing the result of Section 4 and Appendix E, Theorem 8 and Figure 5 (obtaining stochastic gradient descent with primal averaging, using the parameter selection technique).

#### Bounded variance at optimum
- [`A_ParameterSelection`](5_Stochastic_VarianceAtOptimum/A_ParameterSelection.m) Code for recovering the result of Appendix G, Theorem 15 (stochastic gradient descent with primal averaging from the parameter selection technique).


### Stochastic first-order methods: unbiased oracles under weak growth conditions
- - [`A_StochasticGradientDescentWithPrimalAveraging`](5_Stochastic_WeakGrowth/A_StochasticGradientDescentWithPrimalAveraging.m) Code for verifying the result of Appendix F, Theorem 15 (stochastic gradient descent with primal averaging, parameter selection technique).

