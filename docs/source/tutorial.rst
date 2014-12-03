Tutorial
========

Estimate L-moments from sample data
-----------------------------------

The function :meth:`lmoments3.lmom_ratios` is used to calculate a specified number of L-moment ratios for a sample
dataset (`list` or 1-dimensional numpy `array`).

Example:

>>> import lmoments3 as lm
>>> data = [2.0, 3.0, 4.0, 2.4, 5.5, 1.2, 5.4, 2.2, 7.1, 1.3, 1.5]
>>> lm.lmom_ratios(data, nmom=5)
[3.2363636363636363, 1.1418181818181818, 0.27388535031847133, 0.023354564755838598, -0.042462845010615709]

This returns the first two L-moments followed by three L-moment ratios like this: l1, l2, t3, t4, t5. Where t3..5 =
l3..5 / l2.

Fitting distributions to sample data
------------------------------------

Distributions can be fitted to sample data using the :meth:`lmoments3.distr.LmomDistrMixin.lmom_fit` method.

For example, using the gamma distribution:

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

Using distributions
-------------------

:mod:`lmoments3` distributions can be used like any other :class:`scipy.stats.rv_continous` distribution. For full
details of all methods, see the `scipy.stats documentation`_. Some useful methods include:

- `pdf`: Probability density function
- `cdf`: Cumulative distribution function
- `ppf`: Inverse cumulative distribution function (also known as quantile function or percentage point function)
- `rvs`: Random numbers generator

The class :class:`lmoments3.distr.LmomDistrMixin` extends the distribution and provides the following additional
methods:

- `lmom_fit`: Fit distributions to sample data
- `lmom`: Distribution's L-moments
- `lmom_ratios`: Distribution's L-moment ratios

Computing distribution L-moments
--------------------------------

The :mod:`lmoments3` package provides methods to compute the L-moments (λ1..n) and L-moment ratios (λ1, λ2, τ3..n) for a
distribution with given parameters.

.. note::

   When fitting a distribution to (sample) data, the higher order **distribution's** L-moments may be
   different from the **sample** L-moments. For example for a two-parameter distribution, τ3 (the distribution's third
   L-moment) may be different from t3 (the sample's corresponding L-moment).

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

.. _scipy.stats documentation: http://docs.scipy.org/doc/scipy/reference/stats.html