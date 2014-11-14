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
import math


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


class WakebyGen(sp.stats.rv_continuous):
    """
    The Wakeby distribution is defined by the transformation:
    (x-xi)/a = (1/b).[1 - (1-U)^b] - (c/d).[1 - (1-U)^(-d)]

    """

    def _argcheck(self, b, c, d):
        b = np.asarray(b)
        c = np.asarray(c)
        d = np.asarray(d)
        check = np.where(b + d > 0,
                         np.where(c == 0, d == 0, True),
                         (b == c) & (c == d) & (d == 0))
        np.putmask(check, c > 0, d > 0)
        np.putmask(check, c < 0, False)
        return check

    def _ppf(self, q, b, c, d):
        z = -np.log(1. - q)
        u = np.where(b == 0, z, (1. - np.exp(-b * z)) / b)
        v = np.where(d == 0, z, (1. - np.exp(d * z)) / d)
        return u - c * v

    def _cdf(self, x, b, c, d):
        if hasattr(x, '__iter__'):
            return np.array([self._cdfwak(_, parameters)
                             for (_, parameters) in zip(x, zip(b, c, d))])
        else:
            return self._cdfwak(x, (b, c, d))

    def _cdfwak(self, x, para):
        # Only for a single value of x!

        EPS = 1e-8
        MAXIT = 20
        ZINCMX = 3
        ZMULT = 0.2
        UFL = -170
        XI = 0  # sp.stats.rv_continuous deals with scaling
        A = 1   # sp.stats.rv_continuous deals with scaling
        B, C, D = para

        CDFWAK = 0
        if x <= XI:
            return CDFWAK

        # Test for _special cases
        if B == 0 and C == 0 and D == 0:
            Z = (x - XI) / A
            CDFWAK = 1
            if -Z >= UFL:
                CDFWAK = 1 - math.exp(-Z)
            return CDFWAK

        if C == 0:
            CDFWAK = 1
            if x >= (XI + A / B):
                return (CDFWAK)
            Z = -math.log(1 - (x - XI) * B / A) / B
            if -Z >= UFL:
                CDFWAK = 1 - math.exp(-Z)
            return CDFWAK

        if A == 0:
            Z = math.log(1 + (x - XI) * D / C) / D
            if -Z >= UFL:
                CDFWAK = 1 - math.exp(-Z)
            return CDFWAK

        CDFWAK = 1
        if D < 0 and x >= (XI + A / B - C / D):
            return CDFWAK

        Z = 0.7
        if x < self._ppf(0.1, *para):
            Z = 0
        if x < self._ppf(0.99, *para):
            pass
        else:
            if D < 0:
                Z = math.log((x - XI - A / B) * D / C + 1) / D
            if D == 0:
                Z = (x - XI - A / B) / C
            if D > 0:
                Z = math.log((x - XI) * D / C + 1) / D

        for IT in range(1, MAXIT + 1):
            EB = 0
            BZ = -B * Z
            if BZ >= UFL:
                EB = math.exp(BZ)
            GB = Z

            if abs(B) > EPS:
                GB = (1 - EB) / B
            ED = math.exp(D * Z)
            GD = -Z

            if abs(D) > EPS:
                GD = (1 - ED) / D

            XEST = XI + A * GB - C * GD
            FUNC = x - XEST
            DERIV1 = A * EB + C * ED
            DERIV2 = -A * B * EB + C * D * ED
            TEMP = DERIV1 + 0.5 * FUNC * DERIV2 / DERIV1

            if TEMP <= 0:
                TEMP = DERIV1
            ZINC = FUNC / TEMP

            if ZINC > ZINCMX:
                ZINC = ZINCMX

            ZNEW = Z + ZINC

            if ZNEW <= 0:
                Z = Z * ZMULT
            else:
                Z = ZNEW
                if abs(ZINC) <= EPS:
                    CDFWAK = 1
                    if -Z >= UFL:
                        CDFWAK = 1 - math.exp(-Z)
                    return CDFWAK

    def _pdf(self, x, b, c, d):
        t = (1. - self._cdf(x, b, c, d))
        f = t ** (d + 1) / (t ** (b + d) + c)
        return f


Wakeby = WakebyGen(name='wakeby', shapes='beta, gamma, delta')