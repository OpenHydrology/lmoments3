import scipy as _sp
import scipy.special as _spsp
import scipy.stats as _spst
import math as _math
import sys as _sys
from _otherfunct import is_numeric as _is_numeric
import lmoments
##############################################################
#PDF FUNCTIONS
##############################################################
def pdfexp(x,para):
    [U,A] = para
    f = []
    if _is_numeric(x)==True:
        x = [x]
    else:
        x = list(x)

    for i in x:
        Y = (i-U)/A
        if Y <= 0:
            f.append(0)
        else:
            f.append(_sp.exp(-Y)/A)
    if len(f) == 1:
        f = f[0]
        
    return(f)

#############################################################

def pdfgam(x,para):
    [K, theta] = para
    f = []
    if _is_numeric(x)==True:
        x = [x]
    else:
        x = list(x)

    for i in x:
        if i <= 0:
            f.append(0)
        else:
            f.append(1/(_spsp.gamma(K)*theta**K)*i**(K-1)*_sp.exp(-i/theta))
         
    if len(f) == 1:
        f = f[0]
    return(f)

#############################################################

def pdfgev(x,para):
    SMALL =1e-15 
    [XI,A,K] = para
    f = []
    if _is_numeric(x)==True:
        x = [x]
    else:
        x = list(x)
        
    for i in x:
        Y = (i-XI)/A
        if K == 0:
            f.append(A**(-1)*_sp.exp(-(1-K)*Y-_sp.exp(-Y)))
            continue
        ARG = 1-K*Y
        if ARG > SMALL:
            Y = -_sp.log(ARG)/K
            f.append(A**(-1) * _sp.exp(-(1-K)*Y-_sp.exp(-Y)))
            continue
        if K < 0:
            f.append(0)
            continue
        f.append(1)
        
    if len(f) == 1:
        f = f[0]
        
    return(f)

#############################################################

def pdfglo(x,para):
    SMALL = 1e-15
    [XI,A,K] = para
    f = []
    if _is_numeric(x)==True:
        x = [x]
    else:
        x = list(x)
        
    for i in x:
        Y = (i-XI)/A
        if K ==0:
            f.append(A**(-1)* _sp.exp(-(1-K)*Y)/(1+_sp.exp(-Y))**2)
            continue
        ARG = 1-K*Y
        if ARG > SMALL:
            Y = -_sp.log(ARG)/K
            f.append(A**(-1)*_sp.exp(-(1-K)*Y)/(1+_sp.exp(-Y))**2)
            continue
        if K < 0:
            f.append(0)
            continue
        if K > 0:
            f.append(1)
            continue

    if len(f) == 1:
        f = f[0]
        
    return(f)

#############################################################
   
def pdfgno(x,para):
    SMALL = 1e-15
    [XI,A,K] = para
    f = []

    if _is_numeric(x)==True:
        x = [x]
    else:
        x = list(x)
        
    for i in x:
        Y = (i-XI)/A
        if K != 0:
            ARG = 1-K*Y
            if ARG > SMALL:
                Y = -_sp.log(ARG)/K
            else:
                if K < 0:
                    f.append(0)
                    continue
                else:
                    f.append(1)
                    continue

        f.append(_sp.exp(K*Y-Y**2/2)/(A*_sp.sqrt(2*_sp.pi)))

    if len(f) == 1:
        f = f[0]
        
    return(f)

#############################################################

def pdfgpa(x,para):
    SMALL = 1e-15
    [XI,A,K] = para
    f = []
    
    if _is_numeric(x)==True:
        x = [x]
    else:
        x = list(x)

    for i in x:
        Y = (i-XI)/A
        if Y <= 0:
            f.append(0)
            continue
        if K == 0:
            f.append(A**(-1)*_sp.exp(-(1-K)*Y))
        else:
            ARG = 1-K*Y
            if ARG > SMALL:
                Y = -_sp.log(ARG)/K
                f.append(A**(-1)*_sp.exp(-(1-K)*Y))
            else:
                f.append(1)

    if len(f) == 1:
        f = f[0]
        
    return(f)

#############################################################

def pdfgum(x,para):
    [U,A] = para
    f = []

    if _is_numeric(x)==True:
        x = [x]
    else:
        x = list(x)
        
    for i in x:
        Y = -(i-U)/A
        f.append(A**(-1)*_sp.exp(Y)*_sp.exp(-_sp.exp(Y)))
    if len(f) == 1:
        f = f[0]
    return(f)

#############################################################

def pdfkap(x,para):
    SMALL = 1e-15
    [XI,A,K,H] = para
    Fs = lmoments.cdfkap(x, para)
    if _is_numeric(x)==True:
        Fs = [Fs]
        
    f = []

    if _is_numeric(x)==True:
        x = [x]
    else:
        x = list(x)
        
    for i in range(0,len(x)):
        Y = (x[i]-XI)/A
        if K != 0:
            Y = 1-K*Y
            if Y <= SMALL:
                f.append(0)
            else:
                Y = (1-1/K)*_sp.log(Y)

        Y = _sp.exp(-Y)
        f.append(Y/A * Fs[i]**(1-H))
    if len(f) == 1:
        f = f[0]
    return(f)

#############################################################

def pdfnor(x,para):
    f = []
    
    if _is_numeric(x)==True:
        x = [x]
    else:
        x = list(x)
        
    for i in x:
        f.append(_spst.norm(para[0],para[1]).pdf(i))
    if len(f) == 1:
        f = f[0]
    return(f)

#############################################################

def pdfpe3(x,para):    
    [mu, sigma, gamma] = para
    f = []
    
    if _is_numeric(x)==True:
        x = [x]
    else:
        x = list(x)
        
    alpha = 4/gamma**2
    tmp = _spsp.gamma(alpha)
    if gamma == 0 or _math.isinf(tmp):
        for i in x:
            f.append(_spst.norm(0,1).pdf(((i-mu)/sigma)))
        return(f)

    beta = 0.5 * sigma * abs(gamma)
    xi = mu-2*sigma/gamma
    if gamma > 0:
        for i in x:
            Y = i-xi
            f.append((Y)**(alpha-1) * _sp.exp(-Y/beta)/(beta**alpha*tmp))
        if len(f) == 1:
            f = f[0]
        return(f)
    else:
        for i in x:
            Y = xi - i
            f.append((Y)**(alpha-1) * _sp.exp(-Y/beta)/(beta**alpha*tmp))
        if len(f) == 1:
            f = f[0]
        return(f)

#############################################################
    
def pdfwak(x,para):
    [XI,A,B,C,D] = para
    f = []
    
    if _is_numeric(x)==True:
        x = [x]
    else:
        x = list(x)
        
    for i in x:
        Fc = 1-lmoments.cdfwak(i,para)
        if Fc == 0:
            Fc = 1e-10
        tmp = A*Fc**(B - 1) + C * Fc**(-D - 1)
        f.append(1/tmp)
    if len(f) == 1:
        f = f[0]
    return(f)

#############################################################

def pdfwei(x,para):
    f = []
    
    if _is_numeric(x)==True:
        x = [x]
    else:
        x = list(x)
        
    [ZETA,B,D] = para
    K = 1/D
    A = B/D
    XI = ZETA -B
    
    for i in x:
        f.append(lmoments.pdfgev(-i,[XI,A,K]))

    if len(f) == 1:
        f = f[0]
        
    return(f)
        
