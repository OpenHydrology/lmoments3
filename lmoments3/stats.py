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


def NlogL(data, distr_name, distr_paras={}):
    distr_name = distr_name.lower()  # Ignore case

    # Fit data to estimate distribution function parameters, if not provided
    if not distr_paras:
        lmoms = lm.samlmu(data)
        pel_f = getattr(lm, 'pel' + distr_name)
        distr_paras = pel_f(lmoms)

    distr_f = getattr(distr, distr_name)  # scipy rv_continous class
    L = distr_f(**distr_paras).pdf(data)
    NLL = -sum(sp.log(L))
    return NLL

    # if distr_paras:
    #     peldist = distr_paras
    #     if distr_name == "EXP":
    #         pdf = pdfexp
    #     elif distr_name == "GAM":
    #         pdf = pdfgam
    #     elif distr_name == "GEV":
    #         pdf = pdfgev
    #     elif distr_name == "GLO":
    #         pdf = pdfglo
    #     elif distr_name == "GNO":
    #         pdf = pdfgno
    #     elif distr_name == "GPA":
    #         pdf = pdfgpa
    #     elif distr_name == "GUM":
    #         pdf = pdfgum
    #     elif distr_name == "KAP":
    #         pdf = pdfkap
    #     elif distr_name == "NOR":
    #         pdf = pdfnor
    #     elif distr_name == "PE3":
    #         pdf = pdfpe3
    #     elif distr_name == "WAK":
    #         pdf = pdfwak
    #     elif distr_name == "WEI":
    #         pdf = pdfwei
    # else:
    #     if distr_name == "EXP":
    #         pdf = pdfexp
    #         peldist = pelexp(samlmu(data))
    #     elif distr_name == "GAM":
    #         pdf = pdfgam
    #         peldist = pelgam(samlmu(data))
    #     elif distr_name == "GEV":
    #         pdf = pdfgev
    #         peldist = pelgev(samlmu(data))
    #     elif distr_name == "GLO":
    #         pdf = pdfglo
    #         peldist = pelglo(samlmu(data))
    #     elif distr_name == "GNO":
    #         pdf = pdfgno
    #         peldist = pelgno(samlmu(data))
    #     elif distr_name == "GPA":
    #         pdf = pdfgpa
    #         peldist = pelgpa(samlmu(data))
    #     elif distr_name == "GUM":
    #         pdf = pdfgum
    #         peldist = pelgum(samlmu(data))
    #     elif distr_name == "KAP":
    #         pdf = pdfkap
    #         peldist = pelkap(samlmu(data))
    #     elif distr_name == "NOR":
    #         pdf = pdfnor
    #         peldist = pelnor(samlmu(data))
    #     elif distr_name == "PE3":
    #         pdf = pdfpe3
    #         peldist = pelpe3(samlmu(data))
    #     elif distr_name == "WAK":
    #         pdf = pdfwak
    #         peldist = pelwak(samlmu(data))
    #     elif distr_name == "WEI":
    #         pdf = pdfwei
    #         peldist = pelwei(samlmu(data))
    #
    # L = pdf(data, peldist)


def NumParam(dist):
    if dist == "EXP":
        return (2)
    elif dist == "GAM":
        return (2)
    elif dist == "GEV":
        return (3)
    elif dist == "GLO":
        return (3)
    elif dist == "GNO":
        return (3)
    elif dist == "GPA":
        return (3)
    elif dist == "GUM":
        return (2)
    elif dist == "KAP":
        return (4)
    elif dist == "NOR":
        return (2)
    elif dist == "PE3":
        return (3)
    elif dist == "WAK":
        return (5)
    elif dist == "WEI":
        return (3)


def AIC(data, dist, *args):
    if len(args) >= 2:
        print('Invalid Number of Arguments')
        return ()
    elif len(args) == 1:
        peldist = args[0]
        NLL = NlogL(data, dist, peldist)
        k = len(peldist)
        AIC = 2 * k + 2 * NLL
        return (AIC)
    else:
        NLL = NlogL(data, dist)
        k = NumParam(dist)
        AIC = 2 * k + 2 * NLL
        return (AIC)


def AICc(data, dist, *args):
    if len(args) == 0:
        AICbase = AIC(data, dist)
    else:
        AICbase = AIC(data, dist, *args)
    k = NumParam(dist)
    diff = 2 * k * (k + 1) / (len(data) - k - 1)
    AICc = AICbase + diff
    return (AICc)


def BIC(data, dist, *args):
    if len(args) >= 2:
        print('Invalid Number of Arguments')
        return ()
    elif len(args) == 1:
        peldist = args[0]
        NLL = NlogL(data, dist, peldist)
        k = len(peldist)
        BIC = k * sp.log(len(data)) + 2 * NLL
        return (BIC)
    else:
        NLL = NlogL(data, dist)
        k = NumParam(dist)
        BIC = k * sp.log(len(data)) + 2 * NLL
        return (BIC)
