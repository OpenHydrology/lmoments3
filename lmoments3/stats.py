# -*- coding: utf-8 -*-

# lmoments3 library
# Copyright (C) 2012, 2014  J. R. M. Hosking, William Asquith,
# Sam Gillespie, Pierre Gérard-Marchant, Florenz A. P. Hollebrandse
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import scipy as sp
from . import distr
import lmoments3 as lm


def neg_log_lik(x, distr_name, distr_paras={}):
    """
    Calculate the negative log likelihood of a dataset with a distribution
    :param x: fitted dataset
    :type x: array-like
    :param distr_name: 3-letter code of distribution function, e.g. `exp`
    :type distr_name: str
    :param distr_paras: the distribution function's parameters (optional). If not provided, the data will be fitted
                        using the L-moments method.
    :type distr_paras: dict
    :return: negative log likelihood
    :type: float
    """
    x = sp.asarray(x)
    distr_name = distr_name.lower()  # Ignore case

    # Fit x to estimate distribution function parameters, if not provided
    if not distr_paras:
        #lmoms = lm.samlmu(x)
        #pel_f = getattr(lm, 'pel' + distr_name)
        #distr_paras = pel_f(lmoms)
        distr_paras = getattr(distr, distr_name).lmom_fit(x)

    distr_f = getattr(distr, distr_name)  # scipy rv_continous class
    nll = distr_f.nnlf(theta=list(distr_paras.values()), x=x)
    return nll


def distr_n_params(distr_name):
    distr_f = getattr(distr, distr_name)
    return distr_f.numargs + 2  # Include location and scale in addition to shape parameters


def AIC(x, distr_name, distr_paras={}):
    distr_name = distr_name.lower()  # Ignore case

    NLL = neg_log_lik(x, distr_name, distr_paras)
    k = distr_n_params(distr_name)
    AIC = 2 * k + 2 * NLL
    return AIC


def AICc(x, distr_name, distr_paras={}):
    distr_name = distr_name.lower()  # Ignore case

    AICbase = AIC(x, distr_name, distr_paras)
    k = distr_n_params(distr_name)
    diff = 2 * k * (k + 1) / (len(x) - k - 1)
    AICc = AICbase + diff
    return AICc


def BIC(x, distr_name, distr_paras={}):
    distr_name = distr_name.lower()  # Ignore case

    NLL = neg_log_lik(x, distr_name, distr_paras)
    k = distr_n_params(distr_name)
    BIC = k * sp.log(len(x)) + 2 * NLL
    return BIC
