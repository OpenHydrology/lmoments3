lmoments3 Library documentation
===============================

Version |version|

This library was designed to use L-moments to calculate optimal parameters for a number of distributions. This library
extends a number of :mod:`scipy` distributions and provides some additional distributions frequently used in Extreme
Value Analyses.

========================= ================ ============= ===============================================================
Name                      `lmoments3` name `scipy` name  Parameters
========================= ================ ============= ===============================================================
Exponential               `exp`            `expon`       `loc`, `scale`
Gamma                     `gam`            `gamma`       `a`, `loc`, `scale` (The location parameter is not calculated
                                                         using L-moments and assumed to be zero.)
Generalised Extreme Value `gev`            `genextreme`  `c`, `loc`, `scale`
Generalised Logistic      `glo`            n/a           `k`, `loc`, `scale`
Generalised Normal        `gno`            n/a           `k`, `loc`, `scale`
Generalised Pareto        `gpa`            `genpareto`   `c`, `loc`, `scale`
Gumbel                    `gum`            `gumbel_r`    `loc`, `scale`
Kappa                     `kap`            n/a           `k`, `h`, `loc`, `scale`
Normal                    `nor`            `norm`        `loc`, `scale`
Pearson III               `pe3`            `pearson3`    `skew`, `loc`, `scale`
Wakeby                    `wak`            n/a           `beta`, `gamma`, `delta`, `loc`, `scale`
Weibull                   `wei`            `weibull_min` `c`, `loc`, `scale`
========================= ================ ============= ===============================================================

All distributions in the table above are included in the :mod:`lmoments3.distr` module.

L-moment estimation from sample data
------------------------------------

The primary purpose of this library is to estimate L-moments from a sample dataset.

The function :func:`lmoments3.lmom_ratios(data, nmom)` takes an input list or `numpy` array `data` and the number of
L-moments to estimate from the dataset.

Example:

>>> import lmoments3 as lm
>>> data = [2.0, 3.0, 4.0, 2.4, 5.5, 1.2, 5.4, 2.2, 7.1, 1.3, 1.5]
>>> lm.lmom_ratios(data, nmom=5)
[3.2363636363636363, 1.1418181818181818, 0.27388535031847133, 0.023354564755838598, -0.042462845010615709]

This returns the first five sample L-moments, in the structured as l1, l2, t3, t4, t5. Where t3..5 = l3..5 / l2.

Fitting distribution functions to sample data
---------------------------------------------

Sample data can be fitted directly to statistical distribution functions using the :mod:`lmoments3` library.

For example, using the gamma distribution:

>>> import lmoments3 as lm
>>> from lmoments3 import distr
>>> data = [2.0, 3.0, 4.0, 2.4, 5.5, 1.2, 5.4, 2.2, 7.1, 1.3, 1.5]
>>> paras = distr.gam.lmom_fit(data)
>>> paras
OrderedDict([('a', 2.295206110128833), ('loc', 0), ('scale', 1.4100535991436407)])

This returns the distribution's parameters as an :class:`OrderedDict` in the same order as a standard `scipy` list of
distribution function parameters. The distribution parameters can be used, for example, like this:

>>> fitted_gam = distr.gam(**paras)
>>> median = fitted_gam.ppf(0.5)
>>> median
2.7804212925067344

For full details of distribution function methods, see the
`scipy.stats documentation <http://docs.scipy.org/doc/scipy/reference/stats.html>`_. Some useful methods include:

 - `pdf`: Probability density function
 - `cdf`: Cumulative distribution function
 - `ppf`: Inverse cumulative distribution function (also known as quantile function or percentage point function)
 - `rvs`: Random numbers generator

Computing L-moments from distribution parameters
------------------------------------------------

The :mod:`lmoments3` package provides two additional methods to compute the L-moments (λ1..n) or L-moment ratios
(λ1, λ2, τ3..n) for a distribution with given parameters.

Example:

>>> distr.gam.lmom(nmom=3, **paras)
[3.2363636363636363, 1.1418181181569327, 0.24963415541016151]
>>> distr.gam.lmom_ratios(nmom=4, **paras)
[3.2363636363636363, 1.1418181181569327, 0.21862865148182167, 0.13877337951549581]

Or using the frozen distribution:

>>> moments = fitted_gam.lmom(nmom=3)
>>> ratios = fitted_gam.lmom_ratios(nmom=4)

Modified implementation of negative log likelihood function
-----------------------------------------------------------

:meth:`nnlf(data, *args, **kwds)`

Calculates the Negative Log Likelihood. Provide data to calculate the negeative log likelihood. If no distribution
parameters are provided, the `scipy` defaults of `loc=0` and `scale=1` are used.

Example: Calculate the Negative Log Likelihood of a Gamma distribution fitted to `data`:

>>> from lmoments3 import distr
>>> paras = distr.gam.lmom_fit(data)
>>> distr.gam.nnlf(data, **paras)
21.283995091031549

Example:  Calculate the Negative Log Likelihood of a Gamma distribution with parameters 2.5 and 1.0 when fitted to
`data`:

>>> from lmoments3 import distr
>>> from collections import OrderedDict
>>> distr.gam.nnlf(data, a=2.5, scale=1)
22.166452544264637

Other statistical methods
-------------------------

The :mod:`lmoments3.stats` module provides some additional statistical parametes to evaluate fitting of data to
distribution function.

:func:`AIC(data, distr_name, distr_paras)`

Calculate the Akaike Information Criterion (AIC) using the chosen dataset and distribution.

Example: Calculate the Akaike Information Criterion for the weibull distribution using the input dataset `data`:

>>> from lmoments3 import stats, distr
>>> paras = {'loc': 0.67, 'scale': 2.71, 'c': 1.18}
>>> stats.AIC(data, 'wei', paras)
47.500528639652515

Functions :func:`AICc` and :func:`BIC` have a similar structure and calculate the corrected Akaike Information Criterion
and the Bayesian Information Criterion respectively.
