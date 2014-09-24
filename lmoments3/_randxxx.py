"""
Random functions
"""

import scipy as _sp
from ._quaxxx import *


def randexp(n, para):
    return quaexp(_sp.rand(n), para)


def randgam(n, para):
    return quagam(_sp.rand(n), para)


def randgev(n, para):
    return quagev(_sp.rand(n), para)


def randglo(n, para):
    return quaglo(_sp.rand(n), para)


def randgno(n, para):
    return quagno(_sp.rand(n), para)


def randgpa(n, para):
    return quagpa(_sp.rand(n), para)


def randgum(n, para):
    return quagum(_sp.rand(n), para)


def randkap(n, para):
    return quakap(_sp.rand(n), para)


def randnor(n, para):
    return quanor(_sp.rand(n), para)


def randpe3(n, para):
    return quape3(_sp.rand(n), para)


def randwak(n, para):
    return quawak(_sp.rand(n), para)


def randwei(n, para):
    return quawei(_sp.rand(n), para)
