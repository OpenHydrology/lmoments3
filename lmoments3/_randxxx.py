# -*- coding: utf-8 -*-

# lmoments3 library
# Copyright (C) 2012, 2014  J. R. M. Hosking, William Asquith, Sam Gillespie, Pierre Gérard-Marchant,
# Florenz A. P. Hollebrandse
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

"""
Random functions
"""

from ._quaxxx import *
from numpy.random import random


def randexp(n, para):
    return quaexp(random(n), para)


def randgam(n, para):
    return quagam(random(n), para)


def randgev(n, para):
    return quagev(random(n), para)


def randglo(n, para):
    return quaglo(random(n), para)


def randgno(n, para):
    return quagno(random(n), para)


def randgpa(n, para):
    return quagpa(random(n), para)


def randgum(n, para):
    return quagum(random(n), para)


def randkap(n, para):
    return quakap(random(n), para)


def randnor(n, para):
    return quanor(random(n), para)


def randpe3(n, para):
    return quape3(random(n), para)


def randwak(n, para):
    return quawak(random(n), para)


def randwei(n, para):
    return quawei(random(n), para)
