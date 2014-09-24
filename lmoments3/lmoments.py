"""
This library was designed to use L-moments to predict optimal parameters
for a number of distributions.  Distributions supported in this file are
listed below, with their distribution suffix:
    *Exponential (EXP)
    *Gamma (GAM)
    *Generalised Extreme Value (GEV)
    *Generalised Logistic (GLO)
    *Generalised Normal (GNO)
    *Generalised Pareto (GPA)
    *Gumbel (GUM)
    *Kappa (KAP)
    *Normal (NOR)
    *Pearson III (PE3)
    *Wakeby (WAK)
    *Weibull (WEI)

The primary function in this file is the samlmu(x,nmom) function, which takes
an input dataset x and  input of the number of moments to produce the log
moments of that dataset.

For Instance, given a list "Data", if 5 l-moments are needed, the function
would be called by lmoments.samlmu(Data,5)

In this file contains four different functions for using each distribution.
Each function can be called by the prefix FUN with the suffix DIS.

*PEL: (x,nmom):
      Parameter Estimates.  This takes the L-Moments calculated by samlmu()
      and predicts the parameter estimates for that function.
      
      EXAMPLE: Find Wakeby distribution that best fits dataset DATA:

          import lmoments
          para = lmoments.pelwak(lmoments.samlmu(DATA,5))

*QUA: (f,para)
      Quantile Estimates.  This takes the parameter estimates for a
      distribution, and a given Quantile value to calculate the quantile for the
      given function.

      EXAMPLE: Find the Upper Quantile (75%) of the Kappa distribution that
      best fits dataset DATA:

          import lmoments
          para = lmoments.pelkap(lmoments.samlmu(DATA,5))
          UQ = lmoments.quakap(0.75,para)

*LMR: (para,nmom):
      L-Moment Ratios.  This takes the parameter estimates for a distribution
      and calculates nmom L-Moment ratios.

      EXAMPLE: Find 4 lmoment ratios for the Gumbel distribution that
      best fits dataset DATA:

          import lmoments
          para = lmoments.pelgum(lmoments.samlmu(DATA,5))
          LMR = lmoments.lmrgum(para,4)

*CDF: (x,para):
      Cumulative Distribution Function.  This takes the parameter estimates
      for a distribution and calculates the quantile for a given value x.

      EXAMPLE: Find the quantile of the datapoint 6.4 for the Weibull
      Distribution that best fits the dataset DATA:

          import lmoments
          para = lmoments.pelwei(lmoments.samlmu(DATA,5))
          quantile = lmoments.cdfwei(6.4,para)


*PDF: (x,para):
      Probability Distribution Function.  This takes the parameter estimates
      for a distribution and calculates the p value for a given value x.

      EXAMPLE: Find the p-value of the datapoint 6.4 for the Weibull
      Distribution that best fits the dataset DATA:

          import lmoments
          para = lmoments.pelwei(lmoments.samlmu(DATA,5))
          quantile = lmoments.pdfwei(6.4,para)

*LMOM: (para):
      L-Moment Estimation from Parameters.  This function takes the input
      parameters for a given distribution, and attempt to calculate the
      L-Moments that would correspond to this distribution.

      EXAMPLE: Estimate the L-Moments of the Weibull Distribution that has
      parameters (2.5,1.5,0.5)

          import lmoments
          Lmoments = lmoments.lmomwei([2.5,1.5,0.5])

#RAND:  (n,para)
      Random Number Generator for a given function.  This takes a curve fit
      and returns a list of random numbers generated from that distribution

      EXAMPLE: Generate 10 numbers from the weibull distribution that makes
      up the dataset "data"

          import lmoments
          weibullfit = lmoments.pelwei(lmoments.samlmu(data))
          randnums = lmoments.randwei(10,weibullfit)

*NlogL: (data,dist,peldist):
      Calculates the Negative Log Likelihood for use in AIC/AICc/BIC calculations.
      Provide data, and distribution to calculate NlogL.  Can also provide curve
      fitting parameters, but if they aren't provided, then the function will
      generate them via pelxxx(samlmu(data)).  To specify the distribution, use the
      three letters typically assigned.

      EXAMPLE: Calculate the Negative Log Likelihood of a Gamma distribution
      fitted to Data.

          import lmoments
          NLL = lmoments.NlogL(data,"GAM")

      EXAMPLE:  Calculate the Negative Log Likelihood of a Gamma distribution
      with parameters [2.5,1.0] when fitted to Data.

          import lmoments
          NLL = lmoments.NlogL(data,"GAM",[2.5,1.0])

*AIC: (data,dist,*distfit):
      Calculate the Akaike Information Criterion (AIC) using the chosen dataset
      and distribution

      EXAMPLE:  Calculate the Akaike Information Criterion for the weibull
      distribution using the input dataset data:

          import lmoments
          Akaike = AIC(data,"WEI")

This file contains a Python implimentation of the lmoments.f library created by
J. R. M. HOSKING.                                                                                 

The base Fortran code is copyright of the IBM Corperation, and the licensing
information is shown below:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IBM software disclaimer

LMOMENTS: Fortran routines for use with the method of L-moments
Permission to use, copy, modify and distribute this software for any purpose
and without fee is hereby granted, provided that this copyright and permission
notice appear on all copies of the software. The name of the IBM Corporation
may not be used in any advertising or publicity pertaining to the use of the
software. IBM makes no warranty or representations about the suitability of the
software for any purpose. It is provided "AS IS" without any express or implied
warranty, including the implied warranties of merchantability, fitness for a
particular purpose and non-infringement. IBM shall not be liable for any direct,
indirect, _special or consequential damages resulting from the loss of use,
data or projects, whether in an action of contract or tort, arising out of or
in connection with the use or performance of this software.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Additional code from the R library "lmomco" has been converted into Python.
This library was developed by WILLIAM ASQUITH, and was released under the GPL-3
License. Copyright (C) 2012 WILLIAM ASQUITH

The Python translation was conducted by:
    Sam Gillespie
    Numerical Analyst
    C&R Consulting
    Townsville Australia
    September 2013

For more information, or to report bugs, contact:
    sam@candrconsulting.com.au

Licensing for Python Translation:
####################################################
    Copyright (C) 2014 Sam Gillespie

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.Version 0.1.0:
####################################################

"""

import scipy as _sp
import scipy.special as _spsp
import scipy.stats as _spst
import math as _math
import sys as _sys
from ._cdfxxx import *
from ._lmrxxx import *
from ._pelxxx import *
from ._quaxxx import *
from ._pdfxxx import *
from ._lmomxxx import *
from ._randxxx import *

################################################################
##L-MOMENT CALCULATION FUNCTION samlmu
################################################################

def _is_numeric(obj):
    try:
        obj+obj, obj-obj, obj*obj, obj**obj, obj/obj
    except ZeroDivisionError:
        return True
    except Exception:
        return False
    else:
        try:
            a = len(obj)
            if a == 1:
                return True
            return False
        except:
            return True

def _comb(N,k):
    if (k > N) or (N < 0) or (k < 0):
        return 0
    val = 1
    for j in range(min(k, N-k)):
        val = (val*(N-j))//(j+1)
    return val

def samlmu(x,nmom=5):
    if nmom <= 5:
        var = _samlmusmall(x,nmom)
        return(var)
    else:
        var = _samlmularge(x,nmom)
        return(var)

##LARGE can be used to calculate samlmu when nmom > 5, less efficient
##than samlmusmall.
def _samlmularge(x,nmom=5):
    checkx = []
    for i in x:
        if _is_numeric(i):
            checkx.append(i)
    x = checkx
    if nmom <= 0:
        return("Invalid number of Sample L-Moments")
    
    x = sorted(x)
    n = len(x)   

    if n < nmom:
        return("Insufficient length of data for specified nmoments")
    ##Calculate first order
    ##Pretty efficient, no loops
    coefl1 = 1.0/_comb(n,1)
    suml1 = sum(x)
    l = [coefl1*suml1]

    if nmom == 1:
        return(l[0])

    #Setup comb table, where comb[i][x] refers to comb(x,i)
    comb = []
    for i in range(1,nmom):
        comb.append([])
        for j in range(n):
            comb[-1].append(_comb(j,i))

    for mom in range(2,nmom+1):
##        print(mom)
        coefl = 1.0/mom * 1.0/_comb(n,mom)
        xtrans = []
        for i in range(0,n):
            coeftemp = []
            for j in range(0,mom):
                coeftemp.append(1)

            for j in range(0,mom-1):
                coeftemp[j] = coeftemp[j]*comb[mom-j-2][i]

            for j in range(1,mom):
                coeftemp[j] = coeftemp[j]*comb[j-1][n-i-1]

            for j in range(0,mom):
                coeftemp[j] = coeftemp[j]*_comb(mom-1,j)

            for j in range(0,int(0.5*mom)):
                coeftemp[j*2+1] = -coeftemp[j*2+1] 
            coeftemp = sum(coeftemp)
            xtrans.append(x[i]*coeftemp)

        if mom > 2:
            l.append(coefl*sum(xtrans)/l[1])
        else:
            l.append(coefl*sum(xtrans))
    return(l)

def _samlmusmall(x,nmom=5):
    checkx = []
    for i in x:
        if _is_numeric(i):
            checkx.append(i)
    x = checkx
    if nmom <= 0:
        return("Invalid number of Sample L-Moments")
    
    x = sorted(x)
    n = len(x)   

    if n < nmom:
        return("Insufficient length of data for specified nmoments")
    ##Calculate first order
    ##Pretty efficient, no loops
    coefl1 = 1.0/_comb(n,1)
    suml1 = sum(x)
    l1 = coefl1*suml1

    if nmom == 1:
        ret = l1
        return(ret)

    ##Calculate Second order

    #comb terms appear elsewhere, this will decrease calc time
    #for nmom > 2, and shouldn't decrease time for nmom == 2
    #comb(x,1) = x
    #for i in range(1,n+1):
    ##        comb1.append(_comb(i-1,1))
    ##        comb2.append(_comb(n-i,1))
    #Can be simplifed to comb1 = range(0,n)
    
    comb1 = range(0,n)
    comb2 = range(n-1,-1,-1)
    
    coefl2 = 0.5 * 1.0/_comb(n,2)
    xtrans = []
    for i in range(0,n):
        coeftemp = comb1[i]-comb2[i]
        xtrans.append(coeftemp*x[i])

        
    l2 = coefl2 * sum(xtrans)

    if nmom  ==2:
        ret = [l1,l2]
        return(ret)

    ##Calculate Third order
    #comb terms appear elsewhere, this will decrease calc time
    #for nmom > 2, and shouldn't decrease time for nmom == 2
    #comb3 = comb(i-1,2)
    #comb4 = comb3.reverse()
    comb3 = []
    comb4 = []
    for i in range(0,n):
        combtemp = _comb(i,2)
        comb3.append(combtemp)
        comb4.insert(0,combtemp)
   
    
    coefl3 = 1.0/3 * 1.0/_comb(n,3)
    xtrans = []
    for i in range(0,n):
        coeftemp = (comb3[i]-
                    2*comb1[i]*comb2[i] +
                    comb4[i])
        xtrans.append(coeftemp*x[i])

    l3 = coefl3 *sum(xtrans) /l2

    if nmom  ==3:
        ret = [l1,l2,l3]
        return(ret)

    ##Calculate Fourth order
    #comb5 = comb(i-1,3)
    #comb6 = comb(n-i,3)
    comb5 = []
    comb6 = []
    for i in range(0,n):
        combtemp = _comb(i,3)
        comb5.append(combtemp)
        comb6.insert(0,combtemp)
    
    coefl4 = 1.0/4 * 1.0/_comb(n,4)
    xtrans = []
    for i in range(0,n):
        coeftemp = (comb5[i]-
                    3*comb3[i]*comb2[i] +
                    3*comb1[i]*comb4[i] -
                    comb6[i])
        xtrans.append(coeftemp*x[i])

    l4 = coefl4 *sum(xtrans)/l2

    if nmom  ==4:
        ret = [l1,l2,l3,l4]
        return(ret)

    ##Calculate Fifth order
    comb7 = []
    comb8 = []
    for i in range(0,n):
        combtemp = _comb(i,4)
        comb7.append(combtemp)
        comb8.insert(0,combtemp)
        
    coefl5 = 1.0/5 * 1.0/_comb(n,5)
    xtrans = []
    for i in range(0,n):
        coeftemp = (comb7[i]-
                    4*comb5[i]*comb2[i] +
                    6*comb3[i]*comb4[i] -
                    4*comb1[i]*comb6[i] +
                    comb8[i])
        xtrans.append(coeftemp*x[i])

    l5 = coefl5 *sum(xtrans)/l2

    if nmom ==5:
        ret = [l1,l2,l3,l4,l5]
        return(ret)

##############################################################
#MODEL SELECTION AND INFORMATION STATISTICS: AIC
##############################################################
def NlogL(data,dist,*args):

    if len(args) >= 2:
        print("Invalid number of arguments")
    elif len(args) == 1:
        peldist = args[0]
        if dist == "EXP":
            pdf = pdfexp
        elif dist == "GAM":
            pdf = pdfgam
        elif dist == "GEV":
            pdf = pdfgev
        elif dist == "GLO":
            pdf = pdfglo
        elif dist == "GNO":
            pdf = pdfgno
        elif dist == "GPA":
            pdf = pdfgpa
        elif dist == "GUM":
            pdf = pdfgum
        elif dist == "KAP":
            pdf = pdfkap        
        elif dist == "NOR":
            pdf = pdfnor
        elif dist == "PE3":
            pdf = pdfpe3
        elif dist == "WAK":
            pdf = pdfwak
        elif dist == "WEI":
            pdf = pdfwei
    else:
        if dist == "EXP":
            pdf = pdfexp
            peldist = pelexp(samlmu(data))
        elif dist == "GAM":
            pdf = pdfgam
            peldist = pelgam(samlmu(data))
        elif dist == "GEV":
            pdf = pdfgev
            peldist = pelgev(samlmu(data))
        elif dist == "GLO":
            pdf = pdfglo
            peldist = pelglo(samlmu(data))
        elif dist == "GNO":
            pdf = pdfgno
            peldist = pelgno(samlmu(data))
        elif dist == "GPA":
            pdf = pdfgpa
            peldist = pelgpa(samlmu(data))
        elif dist == "GUM":
            pdf = pdfgum
            peldist = pelgum(samlmu(data))
        elif dist == "KAP":
            pdf = pdfkap
            peldist = pelkap(samlmu(data))
        elif dist == "NOR":
            pdf = pdfnor
            peldist = pelnor(samlmu(data))
        elif dist == "PE3":
            pdf = pdfpe3
            peldist = pelpe3(samlmu(data))
        elif dist == "WAK":
            pdf = pdfwak
            peldist = pelwak(samlmu(data))
        elif dist == "WEI":
            pdf = pdfwei
            peldist = pelwei(samlmu(data))

    L = pdf(data,peldist)
    NLL =-sum(_sp.log(L))
    return(NLL)

##############################################################
def NumParam(dist):
    if dist == "EXP":
        return(2)
    elif dist == "GAM":
        return(2)
    elif dist == "GEV":
        return(3)
    elif dist == "GLO":
        return(3)
    elif dist == "GNO":
        return(3)
    elif dist == "GPA":
        return(3)
    elif dist == "GUM":
        return(2)
    elif dist == "KAP":
        return(4)
    elif dist == "NOR":
        return(2)
    elif dist == "PE3":
        return(3)
    elif dist == "WAK":
        return(5)
    elif dist == "WEI":
        return(3)

##############################################################
def AIC(data,dist,*args):
    if len(args) >= 2:
        print('Invalid Number of Arguments')
        return()
    elif len(args) == 1:
        peldist = args[0]
        NLL = NlogL(data,dist,peldist)
        k = len(peldist)
        AIC = 2*k+2*NLL
        return(AIC)
    else:
        NLL = NlogL(data,dist)
        k = NumParam(dist)
        AIC = 2*k+2*NLL
        return(AIC)

##############################################################
def AICc(data,dist,*args):
    if len(args) == 0:
        AICbase = AIC(data,dist)
    else:
        AICbase = AIC(data,dist,*args)
    k = NumParam(dist)
    diff = 2*k*(k+1)/(len(data)-k-1)
    AICc = AICbase + diff
    return(AICc)

##############################################################
def BIC(data,dist,*args):
    if len(args) >= 2:
        print('Invalid Number of Arguments')
        return()
    elif len(args) == 1:
        peldist = args[0]
        NLL = NlogL(data,dist,peldist)
        k = len(peldist)
        BIC = k*_sp.log(len(data))+2*NLL
        return(BIC)
    else:
        NLL = NlogL(data,dist)
        k = NumParam(dist)
        BIC = k*_sp.log(len(data))+2*NLL
        return(BIC)
    
