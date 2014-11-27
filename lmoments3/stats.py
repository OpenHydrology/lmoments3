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


def neg_log_lik(data, distr_name, distr_paras={}):
    """
    Calculate the negative log likelihood of a dataset with a distribution
    :param data: fitted dataset
    :tydatae x: array-like
    :param distr_name: 3-letter code of distribution function, e.g. `exp`
    :type distr_name: str
    :param distr_paras: the distribution function's parameters (optional). If not provided, the data will be fitted
                        using the L-moments method.
    :type distr_paras: dict
    :return: negative log likelihood
    :type: float
    """
    data = sp.asarray(data)
    distr_name = distr_name.lower()  # Ignore case

    # Fit data to estimate distribution function parameters, if not provided
    if not distr_paras:
        distr_paras = getattr(distr, distr_name).lmom_fit(data)

    distr_f = getattr(distr, distr_name)  # scipy rv_continous class
    nll = distr_f.nnlf(theta=list(distr_paras.values()), x=data)
    return nll


def distr_n_params(distr_name):
    distr_f = getattr(distr, distr_name)
    return distr_f.numargs + 2  # Include location and scale in addition to shape parameters


def AIC(data, distr_name, distr_paras={}):
    distr_name = distr_name.lower()  # Ignore case

    NLL = neg_log_lik(data, distr_name, distr_paras)
    k = distr_n_params(distr_name)
    AIC = 2 * k + 2 * NLL
    return AIC


def AICc(data, distr_name, distr_paras={}):
    distr_name = distr_name.lower()  # Ignore case

    AICbase = AIC(data, distr_name, distr_paras)
    k = distr_n_params(distr_name)
    diff = 2 * k * (k + 1) / (len(data) - k - 1)
    AICc = AICbase + diff
    return AICc


def BIC(data, distr_name, distr_paras={}):
    distr_name = distr_name.lower()  # Ignore case

    NLL = neg_log_lik(data, distr_name, distr_paras)
    k = distr_n_params(distr_name)
    BIC = k * sp.log(len(data)) + 2 * NLL
    return BIC
