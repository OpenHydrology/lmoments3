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

exp = scipy.stats.expon
gam = scipy.stats.gamma
gev = scipy.stats.genextreme
gum = scipy.stats.gumbel_r
nor = scipy.stats.norm
pe3 = scipy.stats.pearson3
wei = scipy.stats.weibull_min


class LmomDistrMixin(object):
    def lmom_fit(self, data):
        """
        Fit the distribution function to the given data.

        :param data: Data to use in calculating the distribution parameters
        :type data: array_like
        :returns: Distribution parameters in `scipy` order, e.g. scale, loc, shape
        :rtype: :class:`OrderedDict`
        """
        lmom = lm.samlmu(data, nmom=self.numargs + 2)
        return self._lmom_fit(lmom)

    def _lmom_fit(self, lmom):
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


class GlogisticGen(LmomDistrMixin, scipy.stats.rv_continuous):
    """
    The CDF is given by

    .. math::
       F(x;k) = \\frac{1}{1 + \left[1 - kx\\right]^{1/k}}
    """

    numargs = 1

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

    def _lmom_fit(self, lmom):
        SMALL = 1e-6
        G = -lmom[2]
        if lmom[1] <= 0 or abs(G) >= 1:
            raise ValueError("L-Moments invalid")

        if abs(G) <= SMALL:
            G = 0
            para1 = lmom[0]
            A = lmom[1]
        else:
            GG = G * math.pi / math.sin(G * math.pi)
            A = lmom[1] / GG
            para1 = lmom[0] - A * (1 - GG) / G

        para = OrderedDict([('k', G),
                            ('loc', para1),
                            ('scale', A)])
        return para


glo = GlogisticGen(name='glogistic', shapes='k')


class GennormGen(LmomDistrMixin, scipy.stats.rv_continuous):
    """
    The CDF is given by

    .. math::
       F(x) = \Phi{\left[ -k^{-1} \log\{1 - kx\} \\right]}
    """

    numargs = 1

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

    def _lmom_fit(self, lmom):
        SMALL = 1e-8
        A0 = 0.20466534e+01
        A1 = -0.36544371e+01
        A2 = 0.18396733e+01
        A3 = -0.20360244e+00
        B1 = -0.20182173e+01
        B2 = 0.12420401e+01
        B3 = -0.21741801e+00

        T3 = lmom[2]
        if lmom[1] <= 0 or abs(T3) >= 1:
            raise ValueError("L-Moments invalid")

        if abs(T3) >= 0.95:
            U, A, G = 0, -1, 0
        elif abs(T3) <= SMALL:
            U, A, G = lmom[0], lmom[1] * math.sqrt(math.pi), 0
        else:
            TT = T3 ** 2
            G = -T3 \
                * (A0 + TT * (A1 + TT * (A2 + TT * A3))) \
                / (1 + TT * (B1 + TT * (B2 + TT * B3)))
            E = math.exp(0.5 * G ** 2)
            A = lmom[1] * G / (E * special.erf(0.5 * G))
            U = lmom[0] + A * (E - 1) / G
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

    def _lmom_fit(self, lmom):
        EPS = 1e-6
        MAXIT = 20
        MAXSR = 10
        HSTART = 1.001
        BIG = 10
        OFLEXP = 170
        OFLGAM = 53

        T3 = lmom[2]
        T4 = lmom[3]
        if lmom[1] <= 0 or abs(T3) >= 1 or abs(T4) >= 1 or T4 <= (5 * T3 * T3 - 1) / 4 or T4 >= (5 * T3 * T3 + 1) / 6:
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
                scale = lmom[1] * G * HH / (ALAM2 * GAM)
                loc = lmom[0] - scale / G * (1 - GAM * U1 / HH)
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

    def _lmom_fit(self, lmom):
        if lmom[1] <= 0 or abs(lmom[2]) >= 1 or abs(lmom[3]) >= 1 or abs(lmom[4]) >= 1:
            raise ValueError("Invalid L-Moments")

        ALAM1 = lmom[0]
        ALAM2 = lmom[1]
        ALAM3 = lmom[2] * ALAM2
        ALAM4 = lmom[3] * ALAM2
        ALAM5 = lmom[4] * ALAM2

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
            D = -(1 - 3 * lmom[2]) / (1 + lmom[2])
            C = (1 - D) * (2 - D) * lmom[1]
            B = 0
            A = 0
            XI = lmom[0] - C / (1 - D)
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


class GenParetoGen(LmomDistrMixin, scipy.stats.distributions.genpareto_gen):
    def _lmom_fit(self, lmom):
        T3 = lmom[2]
        if lmom[1] <= 0 or abs(T3) >= 1:
            raise ValueError("L-Moments invalid")

        G = (1 - 3 * T3) / (1 + T3)

        # CHANGE: shape parameter `c` has been negated from original `lmoments` package to be compatible with scipy's GPA
        # distribution function.
        PARA3 = -G
        PARA2 = (1 + G) * (2 + G) * lmom[1]
        PARA1 = lmom[0] - PARA2 / (1 + G)

        para = OrderedDict([('c', PARA3),
                            ('loc', PARA1),
                            ('scale', PARA2)])
        return para


gpa = GenParetoGen(a=0.0, name='genpareto')
