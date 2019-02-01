# Potential functions for first-order methods

This code comes jointly with the following reference:

> [1] Adrien Taylor, and Francis Bach. "Stochastic first-order methods: non-asymptotic and computer-aided analyses via potential functions."

Date:    February 4, 2019

#### Authors

- [**Adrien Taylor**](http://www.di.ens.fr/~ataylor/)
- [**Francis Bach**](https://www.di.ens.fr/~fbach/)


#### Setup

**Note:** This code requires [YALMIP](https://yalmip.github.io/) along with a suitable SDP solver (e.g., Sedumi, SDPT3, Mosek).


## Organization of the code

### Deterministic first-order methods
- [`1_Deterministic_SmoothConvex/AA_GradientDescent.m`](AA_GradientDescent.m) Code for reproducing the result of Section 3.2.1, Theorem 3 (gradient method).
- [``](.m) Code for reproducing the result presented in Appendix C.2, Theorem 9 (proximal gradient method).
- [``](.m) Code for reproducing the first result of Appendix C.3, Theorem 10 (first accelerated gradient method).
- [``](.m) Code for reproducing the second result of Appendix C.3, Theorem 11 (second accelerated gradient method).

### Stochastic first-order methods: unbiased oracles with bounded variance
- [``](.m) Code for reproducing the first result of Section 3.2.2, Theorem 5 (stochastic gradient descent).
- [``](.m) Code for reproducing the second result of Section 3.2.2, Theorem 6  (stochastic gradient descent with averaging).
- [``](.m) Code for reproducing the third exaresultmple of Section 3.2.2, Theorem 7 (stochastic gradient descent with primal averaging).
- [``](.m) Code for reproducing the result of Appendix D.4, Theorem 12 (stochastic gradient descent with evaluation at the averaged iterate).
- [``](.m) Code for reproducing the first result of Appendix D.6, Theorem 13 (stochastic gradient descent with momentum).

### Stochastic first-order methods: unbiased oracles arising from sampling in expectations of smooth convex functions

#### Over-parametrized models
- [``](.m) Code for reproducing the result of Section 4 and Appendix E, Theorem 8 (stochastic gradient descent with primal averaging, using the parameter selection technique).

#### Bounded variance at optimum
- [``](.m) Code for reproducing the result of Appendix G, Theorem 15 (stochastic gradient descent with primal averaging, parameter selection technique).


### Stochastic first-order methods: unbiased oracles under weak growth conditions
- [``](.m) Code for reproducing the result of Appendix F, Theorem 15 (stochastic gradient descent with primal averaging, parameter selection technique).

