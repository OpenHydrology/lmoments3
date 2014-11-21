lmoments3 Library documentation
===============================

This library was designed to use L-moments to calculate optimal parameters for a number of distributions. This library
extends a number of :mod:`scipy` distributions and provides some additional distributions frequently used in Extreme
Value Analyses.

========================= ================ ============= ========================================
Name                      `lmoments3` name `scipy` name  Parameters
========================= ================ ============= ========================================
Exponential               `exp`            `expon`       `loc`, `scale`
Gamma                     `gam`            `gamma`       `a`, `loc`, `scale`
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
========================= ================ ============= ========================================

All distributions in the table above are included in the :mod:`lmoments3.distr` module.

L-moment estimation from sample data
------------------------------------

The primary purpose of this library is to estimate L-moments from a sample dataset.

The method :meth:`lmoments3.samlmu(x, nmom)` takes an input list or `numpy` array `x` and the number of L-moments to
estimate from the dataset.

Example:

>>> import lmoments3 as lm
>>> data = [2.0, 3.0, 4.0, 2.4, 5.5, 1.2, 5.4, 2.2, 7.1, 1.3, 1.5]
>>> lm.samlmu(data, nmom=5)
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
`scipy.stats documentation <http://docs.scipy.org/doc/scipy/reference/stats.html>`_.

Description of functions
------------------------

:meth:`pel(x, nmom)`

Parameter Estimates.  This takes the L-Moments calculated by :meth:`samlmu()`
and predicts the parameter estimates for that function.

Example: Find Wakeby distribution that best fits dataset `data`::

    import lmoments3 as lm
    para = lm.pelwak(lm.samlmu(data,5))


:meth:`qua(f, para)`

Quantile Estimates.  This takes the parameter estimates for a
distribution, and a given quantile value `f` to calculate the quantile for the
given function.

Example: Find the Upper Quantile (75%) of the Kappa distribution that
best fits dataset `data`::

    import lmoments3 as lm
    para = lm.pelkap(lm.samlmu(data, 5))
    UQ = lm.quakap(0.75, para)


:meth:`lmr(para, nmom)`

L-Moment Ratios.  This takes the parameter estimates for a distribution
and calculates nmom L-Moment ratios.

Example: Find 4 lmoment ratios for the Gumbel distribution that
best fits dataset `data`::

    import lmoments3 as lm
    para = lm.pelgum(lm.samlmu(data, 5))
    LMR = lm.lmrgum(para, 4)


:meth:`cdf(x, para)`

Cumulative Distribution Function.  This takes the parameter estimates
for a distribution and calculates the quantile for a given value `x`.

Example: Find the quantile of the datapoint 6.4  for the Weibull
Distribution that best fits the dataset `data`::

    import lmoments3 as lm
    para = lm.pelwei(lm.samlmu(data, 5))
    quantile = lm.cdfwei(6.4, para)


:meth:`pdf(x, para)`

Probability Distribution Function.  This takes the parameter estimates
for a distribution and calculates the p-value for a given value `x`.

Example: Find the p-value of the datapoint 6.4 for the Weibull
Distribution that best fits the dataset `data`::

    import lmoments3 as lm
    para = lm.pelwei(lm.samlmu(data, 5))
    quantile = lm.pdfwei(6.4, para)

:meth:`lmom(para)`

L-Moment Estimation from Parameters.  This function takes the input
parameters for a given distribution, and attempt to calculate the
L-Moments that would correspond to this distribution.

Example: Estimate the L-Moments of the Weibull Distribution that has
parameters (2.5, 1.5, 0.5)::

    import lmoments3 as lm
    Lmoments = lm.lmomwei([2.5, 1.5, 0.5])

:meth:`rand(n, para)`

Random Number Generator for a given function.  This takes a curve fit
and returns a list of random numbers generated from that distribution

Example: Generate 10 numbers from the weibull distribution that makes
up the dataset `data`::

    import lmoments3 as lm
    weibullfit = lm.pelwei(lm.samlmu(data))
    randnums = lm.randwei(10, weibullfit)

:meth:`NlogL(data, dist, peldist)`

Calculates the Negative Log Likelihood for use in AIC/AICc/BIC calculations.
Provide data, and distribution to calculate NlogL.  Can also provide curve
fitting parameters, but if they aren't provided, then the function will
generate them via pelxxx(samlmu(data)).  To specify the distribution, use the
three letters typically assigned.

Example: Calculate the Negative Log Likelihood of a Gamma distribution
fitted to `data`::

    import lmoments3 as lm
    NLL = lm.NlogL(data, 'GAM')

Example:  Calculate the Negative Log Likelihood of a Gamma distribution
with parameters `[2.5, 1.0]` when fitted to `data`::

    import lmoments3 as lm
    NLL = lm.NlogL(data, 'GAM', [2.5, 1.0])

:meth:`AIC(data, dist, *distfit)`

Calculate the Akaike Information Criterion (AIC) using the chosen dataset
and distribution

Example:  Calculate the Akaike Information Criterion for the weibull
distribution using the input dataset `data`::

    import lmoments3 as lm
    Akaike = lm.AIC(data, 'WEI')


Contents:

.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

