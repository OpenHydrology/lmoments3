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

from collections import OrderedDict
import numpy as np
import scipy as sp
import scipy.stats
import scipy.stats.distributions
from scipy import special
import math
import lmoments3 as lm


class LmomDistrMixin(object):
    """
    Mixin class to add L-moment methods to :class:`scipy.stats.rv_continous` distribution functions. Distributions using
    the mixin should override the methods :meth:`._lmom_fit`, :meth:`.lmom`, and :meth:`.lmom_ratios`.
    """

    def lmom_fit(self, data=[], lmom_ratios=[]):
        """
        Fit the distribution function to the given data or given L-moments.

        :param data: Data to use in calculating the distribution parameters
        :type data: array_like
        :param lmom_ratios: L-moments (ratios) l1, l2, t3, t4, .. to use in calculating the distribution parameters
        :type lmom_ratios: array_like
        :returns: Distribution parameters in `scipy` order, e.g. scale, loc, shape
        :rtype: :class:`OrderedDict`
        """
        n_min = self.numargs + 2
        if len(data) > 0:
            if len(data) <= n_min:
                raise ValueError("At least {} data points must be provided.".format(n_min))
            lmom_ratios = lm.samlmu(data, nmom=n_min)
        elif not lmom_ratios:
            raise Exception("Either `data` or `lmom_ratios` must be provided.")
        elif len(lmom_ratios) < n_min:
            raise ValueError("At least {} number of L-moments must be provided.".format(n_min))

        return self._lmom_fit(lmom_ratios)

    def _lmom_fit(self, lmom_ratios):
        raise NotImplementedError

    def lmom(self, *args, **kwds):
        """
        Compute the distribution's L-moments, e.g. l1, l2, l3, l4, ..

        :param args: Distribution parameters in order of shape(s), loc, scale
        :type args: float
        :param kwds: Distribution parameters as named arguments. See :attr:`rv_continous.shapes` for names of shape
                     parameters
        :type kwds: float
        :returns: List of L-moments
        :rtype: list
        """
        raise NotImplementedError

    def lmom_ratios(self, *args, **kwds):
        """
        Compute the distribution's L-moment ratios, e.g. l1, l2, t3, t4, ..

        :param args: Distribution parameters in order of shape(s), loc, scale
        :type args: float
        :param kwds: Distribution parameters as named arguments. See :attr:`rv_continous.shapes` for names of shape
                     parameters
        :type kwds: float
        :returns: List of L-moment ratios
        :rtype: list
        """
        raise NotImplementedError


"""
The following distributions are **not** available in :mod:`scipy.stats`.
"""


class GenlogisticGen(LmomDistrMixin, scipy.stats.rv_continuous):
    """
    The CDF is given by

    .. math::
       F(x;k) = \\frac{1}{1 + \left[1 - kx\\right]^{1/k}}
    """

    def _argcheck(self, k):
        return (k == k)

    def _cdf(self, x, k):
        u = np.where(k == 0, np.exp(-x), (1. - k * x) ** (1. / k))
        return 1. / (1. + u)

    def _pdf(self, x, k):
        u = np.where(k == 0, np.exp(-x), (1. - k * x) ** (1. / k))
        return u ** (1. - k) / (1. + u) ** 2

    def _ppf(self, q, k):
        F = q / (1. - q)
        return np.where(k == 0, np.log(F), (1 - F ** (-k)) / k)

    def _lmom_fit(self, lmom_ratios):
        SMALL = 1e-6
        G = -lmom_ratios[2]
        if lmom_ratios[1] <= 0 or abs(G) >= 1:
            raise ValueError("L-Moments invalid")

        if abs(G) <= SMALL:
            G = 0
            para1 = lmom_ratios[0]
            A = lmom_ratios[1]
        else:
            GG = G * math.pi / math.sin(G * math.pi)
            A = lmom_ratios[1] / GG
            para1 = lmom_ratios[0] - A * (1 - GG) / G

        para = OrderedDict([('k', G),
                            ('loc', para1),
                            ('scale', A)])
        return para


glo = GenlogisticGen(name='glogistic', shapes='k')


class GennormGen(LmomDistrMixin, scipy.stats.rv_continuous):
    """
    The CDF is given by

    .. math::
       F(x) = \Phi{\left[ -k^{-1} \log\{1 - kx\} \\right]}
    """

    def _argcheck(self, k):
        return (k == k)

    def _cdf(self, x, k):
        y = np.where(k == 0, x, -np.log(1. - k * x) / k)
        return 0.5 * (1 + sp.special.erf(y * np.sqrt(0.5)))

    def _pdf(self, x, k):
        u = np.where(k == 0, x, -np.log(1. - k * x) / k)
        return np.exp(k * u - u * u / 2.) / np.sqrt(2 * np.pi)

    def _ppf(self, q, k):
        u = sp.special.ndtri(q)  # Normal distribution's ppf
        return np.where(k == 0, u, (1. - np.exp(-k * u)) / k)

    def _lmom_fit(self, lmom_ratios):
        SMALL = 1e-8
        A0 = 0.20466534e+01
        A1 = -0.36544371e+01
        A2 = 0.18396733e+01
        A3 = -0.20360244e+00
        B1 = -0.20182173e+01
        B2 = 0.12420401e+01
        B3 = -0.21741801e+00

        T3 = lmom_ratios[2]
        if lmom_ratios[1] <= 0 or abs(T3) >= 1:
            raise ValueError("L-Moments invalid")

        if abs(T3) >= 0.95:
            U, A, G = 0, -1, 0
        elif abs(T3) <= SMALL:
            U, A, G = lmom_ratios[0], lmom_ratios[1] * math.sqrt(math.pi), 0
        else:
            TT = T3 ** 2
            G = -T3 \
                * (A0 + TT * (A1 + TT * (A2 + TT * A3))) \
                / (1 + TT * (B1 + TT * (B2 + TT * B3)))
            E = math.exp(0.5 * G ** 2)
            A = lmom_ratios[1] * G / (E * special.erf(0.5 * G))
            U = lmom_ratios[0] + A * (E - 1) / G
        para = OrderedDict([('k', G),
                            ('loc', U),
                            ('scale', A)])
        return para


gno = GennormGen(name='gennorm', shapes='k')


class KappaGen(LmomDistrMixin, scipy.stats.rv_continuous):
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

    def _lmom_fit(self, lmom_ratios):
        EPS = 1e-6
        MAXIT = 20
        MAXSR = 10
        HSTART = 1.001
        BIG = 10
        OFLEXP = 170
        OFLGAM = 53

        T3 = lmom_ratios[2]
        T4 = lmom_ratios[3]
        if lmom_ratios[1] <= 0 or abs(T3) >= 1 or abs(T4) >= 1 or T4 <= (5 * T3 * T3 - 1) / 4 or \
                        T4 >= (5 * T3 * T3 + 1) / 6:
            raise ValueError("L-Moments invalid")

        G = (1 - 3 * T3) / (1 + T3)
        H = HSTART
        Z = G + H * 0.725
        Xdist = BIG

        # Newton-Raphson Iteration
        for it in range(1, MAXIT + 1):
            for i in range(1, MAXSR + 1):
                if G > OFLGAM:
                    raise Exception("Failed to converge")
                if H > 0:
                    U1 = sp.exp(special.gammaln(1 / H) - special.gammaln(1 / H + 1 + G))
                    U2 = sp.exp(special.gammaln(2 / H) - special.gammaln(2 / H + 1 + G))
                    U3 = sp.exp(special.gammaln(3 / H) - special.gammaln(3 / H + 1 + G))
                    U4 = sp.exp(special.gammaln(4 / H) - special.gammaln(4 / H + 1 + G))
                else:
                    U1 = sp.exp(special.gammaln(-1 / H - G) - special.gammaln(-1 / H + 1))
                    U2 = sp.exp(special.gammaln(-2 / H - G) - special.gammaln(-2 / H + 1))
                    U3 = sp.exp(special.gammaln(-3 / H - G) - special.gammaln(-3 / H + 1))
                    U4 = sp.exp(special.gammaln(-4 / H - G) - special.gammaln(-4 / H + 1))

                ALAM2 = U1 - 2 * U2
                ALAM3 = -U1 + 6 * U2 - 6 * U3
                ALAM4 = U1 - 12 * U2 + 30 * U3 - 20 * U4
                if ALAM2 == 0:
                    raise Exception("Failed to converge")
                TAU3 = ALAM3 / ALAM2
                TAU4 = ALAM4 / ALAM2
                E1 = TAU3 - T3
                E2 = TAU4 - T4

                DIST = max(abs(E1), abs(E2))
                if DIST < Xdist:
                    Success = 1
                    break
                else:
                    DEL1 = 0.5 * DEL1
                    DEL2 = 0.5 * DEL2
                    G = XG - DEL1
                    H = XH - DEL2

            if Success == 0:
                raise Exception("Failed to converge")

            # Test for convergence
            if DIST < EPS:
                TEMP = special.gammaln(1 + G)
                if TEMP > OFLEXP:
                    raise Exception("Failed to converge")
                GAM = sp.exp(TEMP)
                TEMP = (1 + G) * sp.log(abs(H))
                if TEMP > OFLEXP:
                    raise Exception("Failed to converge")

                HH = sp.exp(TEMP)
                scale = lmom_ratios[1] * G * HH / (ALAM2 * GAM)
                loc = lmom_ratios[0] - scale / G * (1 - GAM * U1 / HH)
                return OrderedDict([('k', G),
                                    ('h', H),
                                    ('loc', loc),
                                    ('scale', scale)])
            else:
                XG = G
                XH = H
                XZ = Z
                Xdist = DIST
                RHH = 1 / (H ** 2)
                if H > 0:
                    U1G = -U1 * special.psi(1 / H + 1 + G)
                    U2G = -U2 * special.psi(2 / H + 1 + G)
                    U3G = -U3 * special.psi(3 / H + 1 + G)
                    U4G = -U4 * special.psi(4 / H + 1 + G)
                    U1H = RHH * (-U1G - U1 * special.psi(1 / H))
                    U2H = 2 * RHH * (-U2G - U2 * special.psi(2 / H))
                    U3H = 3 * RHH * (-U3G - U3 * special.psi(3 / H))
                    U4H = 4 * RHH * (-U4G - U4 * special.psi(4 / H))
                else:
                    U1G = -U1 * special.psi(-1 / H - G)
                    U2G = -U2 * special.psi(-2 / H - G)
                    U3G = -U3 * special.psi(-3 / H - G)
                    U4G = -U4 * special.psi(-4 / H - G)
                    U1H = RHH * (-U1G - U1 * special.psi(-1 / H + 1))
                    U2H = 2 * RHH * (-U2G - U2 * special.psi(-2 / H + 1))
                    U3H = 3 * RHH * (-U3G - U3 * special.psi(-3 / H + 1))
                    U4H = 4 * RHH * (-U4G - U4 * special.psi(-4 / H + 1))

                DL2G = U1G - 2 * U2G
                DL2H = U1H - 2 * U2H
                DL3G = -U1G + 6 * U2G - 6 * U3G
                DL3H = -U1H + 6 * U2H - 6 * U3H
                DL4G = U1G - 12 * U2G + 30 * U3G - 20 * U4G
                DL4H = U1H - 12 * U2H + 30 * U3H - 20 * U4H
                D11 = (DL3G - TAU3 * DL2G) / ALAM2
                D12 = (DL3H - TAU3 * DL2H) / ALAM2
                D21 = (DL4G - TAU4 * DL2G) / ALAM2
                D22 = (DL4H - TAU4 * DL2H) / ALAM2
                DET = D11 * D22 - D12 * D21
                H11 = D22 / DET
                H12 = -D12 / DET
                H21 = -D21 / DET
                H22 = D11 / DET
                DEL1 = E1 * H11 + E2 * H12
                DEL2 = E1 * H21 + E2 * H22

                # TAKE NEXT N-R STEP
                G = XG - DEL1
                H = XH - DEL2
                Z = G + H * 0.725

                # REDUCE STEP IF G AND H ARE OUTSIDE THE PARAMETER _spACE

                FACTOR = 1
                if G <= -1:
                    FACTOR = 0.8 * (XG + 1) / DEL1
                if H <= -1:
                    FACTOR = min(FACTOR, 0.8 * (XH + 1) / DEL2)
                if Z <= -1:
                    FACTOR = min(FACTOR, 0.8 * (XZ + 1) / (XZ - Z))
                if H <= 0 and G * H <= -1:
                    FACTOR = min(FACTOR, 0.8 * (XG * XH + 1) / (XG * XH - G * H))

                if FACTOR == 1:
                    pass
                else:
                    DEL1 = DEL1 * FACTOR
                    DEL2 = DEL2 * FACTOR
                    G = XG - DEL1
                    H = XH - DEL2
                    Z = G + H * 0.725


kap = KappaGen(name='kappa', shapes='k, h')


class WakebyGen(LmomDistrMixin, scipy.stats.rv_continuous):
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
            if hasattr(b, '__iter__'):
                # Assume x, b, c, d are arrays with matching length
                result = np.array([self._cdfwak(_, parameters)
                                   for (_, parameters) in zip(x, zip(b, c, d))])
            else:
                # Only x is an array, paras are scalars
                result = np.array([self._cdfwak(_, [b, c, d])
                                   for _ in x])
        else:
            result = self._cdfwak(x, (b, c, d))
        return result

    def _cdfwak(self, x, para):
        # Only for a single value of x!

        EPS = 1e-8
        MAXIT = 20
        ZINCMX = 3
        ZMULT = 0.2
        UFL = -170
        XI = 0  # stats.rv_continuous deals with scaling
        A = 1  # stats.rv_continuous deals with scaling
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

    def _lmom_fit(self, lmom_ratios):
        if lmom_ratios[1] <= 0 or abs(lmom_ratios[2]) >= 1 or abs(lmom_ratios[3]) >= 1 or abs(lmom_ratios[4]) >= 1:
            raise ValueError("Invalid L-Moments")

        ALAM1 = lmom_ratios[0]
        ALAM2 = lmom_ratios[1]
        ALAM3 = lmom_ratios[2] * ALAM2
        ALAM4 = lmom_ratios[3] * ALAM2
        ALAM5 = lmom_ratios[4] * ALAM2

        XN1 = 3 * ALAM2 - 25 * ALAM3 + 32 * ALAM4
        XN2 = -3 * ALAM2 + 5 * ALAM3 + 8 * ALAM4
        XN3 = 3 * ALAM2 + 5 * ALAM3 + 2 * ALAM4
        XC1 = 7 * ALAM2 - 85 * ALAM3 + 203 * ALAM4 - 125 * ALAM5
        XC2 = -7 * ALAM2 + 25 * ALAM3 + 7 * ALAM4 - 25 * ALAM5
        XC3 = 7 * ALAM2 + 5 * ALAM3 - 7 * ALAM4 - 5 * ALAM5

        XA = XN2 * XC3 - XC2 * XN3
        XB = XN1 * XC3 - XC1 * XN3
        XC = XN1 * XC2 - XC1 * XN2
        DISC = XB * XB - 4 * XA * XC
        skip20 = False
        if DISC >= 0:
            DISC = math.sqrt(DISC)
            ROOT1 = 0.5 * (-XB + DISC) / XA
            ROOT2 = 0.5 * (-XB - DISC) / XA
            B = max(ROOT1, ROOT2)
            D = -min(ROOT1, ROOT2)
            if D < 1:
                A = (1 + B) * (2 + B) * (3 + B) / (4 * (B + D)) * ((1 + D) * ALAM2 - (3 - D) * ALAM3)
                C = -(1 - D) * (2 - D) * (3 - D) / (4 * (B + D)) * ((1 - B) * ALAM2 - (3 + B) * ALAM3)
                XI = ALAM1 - A / (1 + B) - C / (1 - D)
                skip20 = bool(C >= 0 and A + C >= 0)

        if not skip20:
            D = -(1 - 3 * lmom_ratios[2]) / (1 + lmom_ratios[2])
            C = (1 - D) * (2 - D) * lmom_ratios[1]
            B = 0
            A = 0
            XI = lmom_ratios[0] - C / (1 - D)
            if D <= 0:
                A = C
                B = -D
                C = 0
                D = 0

        para = OrderedDict([('beta', B),
                            ('gamma', C),
                            ('delta', D),
                            ('loc', XI),
                            ('scale', A)])
        return para


wak = WakebyGen(name='wakeby', shapes='beta, gamma, delta')

"""
The following distributions are available in `scipy.stats` and are redefined here with an `LmomDistrMixin` to extend the
scipy distribution with L-moment methods.
"""


class GenParetoGen(LmomDistrMixin, scipy.stats.distributions.genpareto_gen):
    def _lmom_fit(self, lmom_ratios):
        T3 = lmom_ratios[2]
        if lmom_ratios[1] <= 0 or abs(T3) >= 1:
            raise ValueError("L-Moments invalid")

        G = (1 - 3 * T3) / (1 + T3)

        # CHANGE: shape parameter `c` has been negated from original `lmoments` package to be compatible with scipy's
        # GPA distribution function.
        PARA3 = -G
        PARA2 = (1 + G) * (2 + G) * lmom_ratios[1]
        PARA1 = lmom_ratios[0] - PARA2 / (1 + G)

        para = OrderedDict([('c', PARA3),
                            ('loc', PARA1),
                            ('scale', PARA2)])
        return para


gpa = GenParetoGen(a=0.0, name='genpareto', shapes='c')


class ExponGen(LmomDistrMixin, scipy.stats.distributions.expon_gen):
    def _lmom_fit(self, lmom_ratios):
        if lmom_ratios[1] <= 0:
            raise ValueError("L-Moments invalid")

        para = OrderedDict([('loc', lmom_ratios[0] - 2 * lmom_ratios[1]),
                            ('scale', 2 * lmom_ratios[1])])
        return para


exp = ExponGen(a=0.0, name='expon')


class GammaGen(LmomDistrMixin, scipy.stats.distributions.gamma_gen):
    def _lmom_fit(self, lmom_ratios):
        A1 = -0.3080
        A2 = -0.05812
        A3 = 0.01765
        B1 = 0.7213
        B2 = -0.5947
        B3 = -2.1817
        B4 = 1.2113

        if lmom_ratios[0] <= lmom_ratios[1] or lmom_ratios[1] <= 0:
            raise ValueError("L-Moments invalid")

        CV = lmom_ratios[1] / lmom_ratios[0]
        if CV >= 0.5:
            T = 1 - CV
            ALPHA = T * (B1 + T * B2) / (1 + T * (B3 + T * B4))
        else:
            T = math.pi * CV ** 2
            ALPHA = (1 + A1 * T) / (T * (1 + T * (A2 + T * A3)))

        para = OrderedDict([('a', ALPHA),
                            ('loc', 0),
                            ('scale', lmom_ratios[0] / ALPHA)])
        return para


gam = GammaGen(a=0.0, name='gamma', shapes='a')


class GenextremeGen(LmomDistrMixin, scipy.stats.distributions.genextreme_gen):
    def _lmom_fit(self, lmom_ratios):
        SMALL = 1e-5
        eps = 1e-6
        maxit = 20
        EU = 0.57721566
        DL2 = math.log(2)
        DL3 = math.log(3)
        A0 = 0.28377530
        A1 = -1.21096399
        A2 = -2.50728214
        A3 = -1.13455566
        A4 = -0.07138022
        B1 = 2.06189696
        B2 = 1.31912239
        B3 = 0.25077104
        C1 = 1.59921491
        C2 = -0.48832213
        C3 = 0.01573152
        D1 = -0.64363929
        D2 = 0.08985247

        T3 = lmom_ratios[2]
        if lmom_ratios[1] <= 0 or abs(T3) >= 1:
            raise ValueError("L-Moments Invalid")

        if T3 <= 0:
            G = (A0 + T3 * (A1 + T3 * (A2 + T3 * (A3 + T3 * A4)))) / (1 + T3 * (B1 + T3 * (B2 + T3 * B3)))

            if T3 >= -0.8:
                para3 = G
                GAM = math.exp(special.gammaln(1 + G))
                para2 = lmom_ratios[1] * G / (GAM * (1 - 2 ** (-G)))
                para1 = lmom_ratios[0] - para2 * (1 - GAM) / G
                para = [para1, para2, para3]
                return para
            elif T3 <= -0.97:
                G = 1 - math.log(1 + T3) / DL2

            T0 = (T3 + 3) * 0.5
            for IT in range(1, maxit):
                X2 = 2 ** (-G)
                X3 = 3 ** (-G)
                XX2 = 1 - X2
                XX3 = 1 - X3
                T = XX3 / XX2
                DERIV = (XX2 * X3 * DL3 - XX3 * X2 * DL2) / (XX2 ** 2)
                GOLD = G
                G = G - (T - T0) / DERIV

                if abs(G - GOLD) <= eps * G:
                    para3 = G
                    GAM = math.exp(special.gammaln(1 + G))
                    para2 = lmom_ratios[1] * G / (GAM * (1 - 2 ** (-G)))
                    para1 = lmom_ratios[0] - para2 * (1 - GAM) / G
                    para = OrderedDict([('c', para3),
                                        ('loc', para1),
                                        ('scale', para2)])
                    return para
            raise Exception("Iteration has not converged")
        else:
            Z = 1 - T3
            G = (-1 + Z * (C1 + Z * (C2 + Z * C3))) / (1 + Z * (D1 + Z * D2))
            if abs(G) < SMALL:
                para2 = lmom_ratios[1] / DL2
                para1 = lmom_ratios[0] - EU * para2
                para = OrderedDict([('c', 0),
                                    ('loc', para1),
                                    ('scale', para2)])
            else:
                para3 = G
                GAM = math.exp(special.gammaln(1 + G))
                para2 = lmom_ratios[1] * G / (GAM * (1 - 2 ** (-G)))
                para1 = lmom_ratios[0] - para2 * (1 - GAM) / G
                para = OrderedDict([('c', para3),
                                    ('loc', para1),
                                    ('scale', para2)])
            return para


gev = GenextremeGen(name='genextreme', shapes='c')


class GumbelGen(LmomDistrMixin, scipy.stats.distributions.gumbel_r_gen):
    def _lmom_fit(self, lmom_ratios):
        EU = 0.577215664901532861
        if lmom_ratios[1] <= 0:
            raise ValueError("L-Moments Invalid")

        para2 = lmom_ratios[1] / math.log(2)
        para1 = lmom_ratios[0] - EU * para2
        para = OrderedDict([('loc', para1),
                            ('scale', para2)])
        return para


gum = GumbelGen(name='gumbel_r')


class NormGen(LmomDistrMixin, scipy.stats.distributions.norm_gen):
    def _lmom_fit(self, lmom_ratios):
        if lmom_ratios[1] <= 0:
            raise ValueError("L-Moments invalid")

        para = OrderedDict([('loc', lmom_ratios[0]),
                            ('scale', lmom_ratios[1] * math.sqrt(math.pi))])
        return para


nor = NormGen(name='norm')


class Pearson3Gen(LmomDistrMixin, scipy.stats.distributions.pearson3_gen):
    def _lmom_fit(self, lmom_ratios):
        Small = 1e-6
        # Constants used in Minimax Approx:
        C1 = 0.2906
        C2 = 0.1882
        C3 = 0.0442
        D1 = 0.36067
        D2 = -0.59567
        D3 = 0.25361
        D4 = -2.78861
        D5 = 2.56096
        D6 = -0.77045

        T3 = abs(lmom_ratios[2])
        if lmom_ratios[1] <= 0 or T3 >= 1:
            raise ValueError("L-Moments invalid")

        if T3 <= Small:
            loc = lmom_ratios[0]
            scale = lmom_ratios[1] * math.sqrt(math.pi)
            skew = 0
        else:
            if T3 >= (1.0 / 3):
                T = 1 - T3
                Alpha = T * (D1 + T * (D2 + T * D3)) / (1 + T * (D4 + T * (D5 + T * D6)))
            else:
                T = 3 * math.pi * T3 * T3
                Alpha = (1 + C1 * T) / (T * (1 + T * (C2 + T * C3)))

            RTALPH = math.sqrt(Alpha)
            BETA = math.sqrt(math.pi) * lmom_ratios[1] * math.exp(special.gammaln(Alpha) - special.gammaln(Alpha + 0.5))
            loc = lmom_ratios[0]
            scale = BETA * RTALPH
            skew = 2 / RTALPH
            if lmom_ratios[2] < 0:
                skew *= -1

        return OrderedDict([('skew', skew),
                            ('loc', loc),
                            ('scale', scale)])


pe3 = Pearson3Gen(name="pearson3", shapes='skew')


class FrechetRGen(LmomDistrMixin, scipy.stats.distributions.frechet_r_gen):
    def _lmom_fit(self, lmom_ratios):
        if lmom_ratios[1] <= 0 or lmom_ratios[2] >= 1 or lmom_ratios[2] <= -lm.lmrgum([0, 1], nmom=3)[2]:
            raise ValueError("L-Moments invalid")

        pg = gev.lmom_fit(lmom_ratios=[-lmom_ratios[0], lmom_ratios[1], -lmom_ratios[2]])
        delta = 1 / pg[2]
        beta = pg[1] / pg[2]
        para = OrderedDict([('c', delta),
                            ('loc', -pg[0] - beta),
                            ('scale', beta)])
        return para


wei = FrechetRGen(a=0.0, name='weibull_min', shapes='c')
