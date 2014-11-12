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


import numpy as np
import scipy as sp


class KappaGen(sp.stats.rv_continuous):
    """
    The CDF is given by

    .. math::
       F(x; a, b) = \left[1-h\{1-kx\}^{1/k}\\right]^{1/h}
    """

    def _argcheck(self, k, h):
        k = np.asarray(k)
        h = np.asarray(h)
        # Upper bound
        self.b = np.where(k <= 0, np.inf, 1. / k)
        # Lower bound
        self.a = np.where(h > 0,
                          np.where(k == 0, 0., (1 - h ** (-k)) / k),
                          np.where(k < 0, 1. / k, -np.inf))
        return (k == k) | (h == h)

    def _cdf(self, x, k, h):
        y = np.where(k == 0, np.exp(-x), (1 - k * x) ** (1. / k))
        return np.where(h == 0, np.exp(-y), (1. - h * y) ** (1. / h))

    def _pdf(self, x, k, h):
        y = (1 - k * x) ** (1. / k - 1.)
        y *= self._cdf(x, k, h) ** (1. - h)
        return y

    def _ppf(self, q, k, h):
        y = np.where(h == 0, -np.log(q), (1. - q ** h) / h)
        y = np.where(k == 0, -np.log(y), (1. - y ** k) / k)
        return y


Kappa = KappaGen(name='kappa', shapes='k,h')