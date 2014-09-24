lmoments3 Library documentation
===============================

This library was designed to use L-moments to predict optimal parameters
for a number of distributions.  Distributions supported in this file are
listed below, with their distribution suffix:

* Exponential (EXP)
* Gamma (GAM)
* Generalised Extreme Value (GEV)
* Generalised Logistic (GLO)
* Generalised Normal (GNO)
* Generalised Pareto (GPA)
* Gumbel (GUM)
* Kappa (KAP)
* Normal (NOR)
* Pearson III (PE3)
* Wakeby (WAK)
* Weibull (WEI)

The primary function in this file is the :meth:`samlmu(x, nmom)` function, which takes
an input dataset `x` and input of the number of moments to produce the log
moments of that dataset.

For Instance, given a list `Data`, if 5 L-moments are needed, the function
would be called by `lm.samlmu(Data,5)`

In this file contains four different functions for using each distribution.
Each function can be called by the prefix `FUN` with the suffix `DIS`.

Description of functions
------------------------

:meth:`pel(x, nmom)`

Parameter Estimates.  This takes the L-Moments calculated by :meth:`samlmu()`
and predicts the parameter estimates for that function.

Example: Find Wakeby distribution that best fits dataset `data`::

    import lmoments3 as lm as lm
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

