Probability distributions
=========================

The probablity distributions supported by :mod:`lmoments3` are summarised in the table below.


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

All distributions in the table above are included in the :mod:`lmoments3.distr` module. They can be used like any other
:class:`scipy.stats.rv_continous` distribution (`documentation`_).

Basic usage is as follows:

>>> from lmoments3 import distr
>>> general_exp_distr = distr.exp

Or a distribution with parameters "frozen":

>>> exp_distr = distr.exp(loc=0, scale=2)

.. _documentation: http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.html#scipy.stats.rv_continuous