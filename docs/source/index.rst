lmoments3 Library documentation
===============================

   "In statistics, L-moments are a sequence of statistics used to summarize the shape of a probability distribution.
   They are L-statistics (linear combinations of order statistics, hence the "L" for "linear") analogous to conventional
   moments, and can be used to calculate quantities analogous to standard deviation, skewness and kurtosis, termed the
   L-scale, L-skewness and L-kurtosis respectively (the L-mean is identical to the conventional mean). Standardised
   L-moments are called L-moment ratios and are analogous to standardized moments. Just as for conventional moments, a
   theoretical distribution has a set of population L-moments. Sample L-moments can be defined for a sample from the
   population, and can be used as estimators of the population L-moments."

   -- `Wikipedia`_

The :mod:`lmoments3` library can be used to:

1. Calculate L-moments from sample data.
2. Fit probability distributions to sample data using L-moments.

The library extends the `scipy.stats module`_ with additional methods. Supported distributions are described in
:doc:`/distributions`.

Contents:

.. toctree::
   :numbered:
   :maxdepth: 2

   tutorial
   distributions
   reference

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. _Wikipedia: https://en.wikipedia.org/wiki/L-moment
.. _scipy.stats module: http://docs.scipy.org/doc/scipy/reference/stats.html