import scipy as _sp
import scipy.special as _spsp
import scipy.stats as _spst
import math as _math
import sys as _sys
import random
import lmoments
#############################################################
##RANDOM FUNCTIONS
#############################################################

def randexp(n,para):
    return(lmoments.quaexp(_sp.rand(n),para))

def randgam(n,para):
    return(lmoments.quagam(_sp.rand(n),para))

def randgev(n,para):
    return(lmoments.quagev(_sp.rand(n),para))

def randglo(n,para):
    return(lmoments.quaglo(_sp.rand(n),para))

def randgno(n,para):
    return(lmoments.quagno(_sp.rand(n),para))

def randgpa(n,para):
    return(lmoments.quagpa(_sp.rand(n),para))

def randgum(n,para):
    return(lmoments.quagum(_sp.rand(n),para))

def randkap(n,para):
    return(lmoments.quakap(_sp.rand(n),para))

def randnor(n,para):
    return(lmoments.quanor(_sp.rand(n),para))

def randpe3(n,para):
    return(lmoments.quape3(_sp.rand(n),para))

def randwak(n,para):
    return(lmoments.quawak(_sp.rand(n),para))

def randwei(n,para):
    return(lmoments.quawei(_sp.rand(n),para))
