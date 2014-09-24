import scipy as _sp
import scipy.special as _spsp
import scipy.stats as _spst
import math as _math
import sys as _sys
from lmoments import *
from _otherfunct import is_numeric as _is_numeric
#######################################################
#CDF FUNCTIONS
#######################################################

def cdfexp(x,para):
    if _is_numeric(x) == False:
        x = _sp.array(x)

    U = para[0]
    A = para[1]
    if A <= 0:
        cdfexp = 0
        print("Parameters Invalid")
        return(cdfexp)
    else:
        Y = (x-U)/A
        if U <= 0:
            cdfexp = 0
            print("Parameters Invalid")
            return(cdfexp)
        else:
            cdfexp = 1-_sp.exp(-Y)
            if _is_numeric(x)==True:
                if cdfexp >= 0:
                    return(cdfexp)
                else:
                    return(0)
            else:
                for i in range(0,len(cdfexp)):
                    if cdfexp[i] < 0:
                        cdfexp[i] = 0

            return(cdfexp)
            
#############################################################
            
def cdfgam(x,para):

    CDFGAM=0
    Alpha=para[0]
    Beta=para[1]
    if Alpha <= 0 or Beta <= 0:
        print("Parameters Invalid")
        return
    if _is_numeric(x)==True:
        if x <= 0:
            print("x Parameter Invalid")
            return
    else:
        for i in x:
            if i <= 0:
                print("One X Parameter in list is Invalid")
                return
        x = _sp.array(x)
    CDFGAM = _spsp.gammainc(Alpha,x/Beta)
    return(CDFGAM)

#############################################################

def cdfgev(x,para):
    SMALL = 1e-15
    U=para[0]
    A=para[1]
    G=para[2]
    if A <= 0:
        print("Parameters Invalid")
        return

    if _is_numeric(x)==False:
        x = _sp.array(x)

        
    Y = (x-U)/A
    if G==0:
        CDFGEV = _sp.exp(-_sp.exp(-Y))
    else:
        Arg = 1-G*Y

        if _is_numeric(x)==True:
            if Arg > SMALL:
                Y = -_sp.log(Arg)/G
                CDFGEV = _sp.exp(-_sp.exp(-Y))
            elif G<0:
                CDFGEV = 0
            else:
                CDFGEV = 1
        else:
            CDFGEV = []
            for i in range(0,len(Y)):
                if Arg[i] > SMALL:
                    Y[i] = -_sp.log(Arg[i])/G
                    CDFGEV.append(_sp.exp(-_sp.exp(-Y[i])))
                elif G<0:
                    CDFGEV.append(0)
                else:
                    CDFGEV.append(1)

            CDFGEV = list(CDFGEV)

    return(CDFGEV)

#############################################################

def cdfglo(x,para):
    SMALL = 1e-15
    U=para[0]
    A=para[1]
    G=para[2]

    if _is_numeric(x)==False:
        x = _sp.array(x)
        
    if A <= 0:
        print("Parameters Invalid")
        return
    Y = (x-U)/A
    if G==0:
        CDFGLO=1/(1+_sp.exp(-Y))
    else:
        Arg = 1-G*Y
        if _is_numeric(x)==True:
            if Arg > SMALL:
                Y = -_sp.log(Arg)/G
                CDFGLO=1/(1+_sp.exp(-Y))
            elif G<0:
                CDFGLO = 0
            else: 
                CDFGLO = 1
        else:
            CDFGLO = []
            for i in range(0,len(Y)):
                if Arg[i] > SMALL:
                    Y[i] = -_sp.log(Arg[i])/G
                    CDFGLO.append(1/(1+_sp.exp(-Y[i])))
                elif G<0:
                    CDFGLO = 0
                else:
                    CDFGLO = 1

            CDFGLO = list(CDFGLO)
            
    return(CDFGLO)

#############################################################

def cdfgno(x,para):
    SMALL = 1e-15
    U=para[0]
    A=para[1]
    G=para[2]
    
    if _is_numeric(x)==False:
        x = _sp.array(x)
        
    if A <= 0:
        print("Parameters Invalid")
        return
    Y = (x-U)/A
    if G==0:
        CDFGNO = 0.5+0.5*_sp.erg(Y*1/_sp.sqrt(2))
    else:
        Arg = 1-G*Y
        if _is_numeric(x)==True:
            if Arg > SMALL:
                Y = -_sp.log(Arg)/G
                CDFGNO = 0.5+0.5*_spsp.erf(Y*1/_sp.sqrt(2))
            elif G<0:
                CDFGNO = 0
            else:
                CDFGNO = 1

        else:
            CDFGNO = []
            for i in range(0,len(Y)):
                if Arg[i] > SMALL:
                    Y[i] = -_sp.log(Arg[i])/G
                    CDFGNO.append(0.5+0.5*_spsp.erf(Y[i]*1/_sp.sqrt(2)))
                elif G<0:
                    CDFGNO.append(0)
                else:
                    CDFGNO.append(1)
            CDFGNO = list(CDFGNO)
    return(CDFGNO)

#############################################################

def cdfgpa(x,para):
    SMALL = 1e-15
    U=para[0]
    A=para[1]
    G=para[2]
    
    if _is_numeric(x)==False:
        x = _sp.array(x)
        
    CDFGPA = 0
    if A <= 0:
        print("Parameters Invalid")
        return
    Y = (x-U)/A

    if _is_numeric(x)==True:
        if Y <= 0:
            print("Parameters Invalid")
            return
    else:
        for i in Y:
            if i <= 0:
                print("Parameters Invalid")
                return
    
    if G==0:
        CDFGPA=1-_sp.exp(-Y)
    else:
        Arg = 1-G*Y
        if _is_numeric(x)==True:
            if Arg > SMALL:
                Y = -_sp.log(Arg)/G
                CDFGPA=1-_sp.exp(-Y)
            else:
                CDFGPA = 1
        else:
            CDFGPA = []
            for i in range(0,len(Y)):
                if Arg[i] > SMALL:
                    Y[i] = -_sp.log(Arg[i])/G
                    CDFGPA.append(1-_sp.exp(-Y[i]))
                else:
                    CDFGPA.append(1)

            CDFGPA = list(CDFGPA)
    return(CDFGPA)

#############################################################

def cdfgum(x,para):
    U = para[0]
    A = para[1]

    if _is_numeric(x)==False:
        x = _sp.array(x)
        
    if A <= 0:
        print("Parameters Invalid")
        return
    else:
        Y = (x-U)/A
        CDFGUM = _sp.exp(-_sp.exp(-Y))
        return(CDFGUM)

#############################################################
    
def cdfkap(x,para):
    para = map(float,para)
    if _is_numeric(x)==True:
        pass
    elif len(x) > 1:
        x = _sp.array(x)
        
    SMALL = 1e-15
    [U,A,G,H] = para
    if A <= 0:
        print("Invalid Parameters")
        return
    if _is_numeric(x)==True:
        Y = [(x-U)/A]
    else:
        Y = list((x-U)/A)
        
    CDFKAP = []
    for i in Y:
        if G != 0:
            ARG = 1-G*i
            if ARG < SMALL:
                if G < 0:
                    CDFKAP.append(0)
                    continue
                if G > 0:
                    CDFKAP.append(1)
                    continue    
            i = -_sp.log(ARG)/G

        i = _sp.exp(-i)
        if H == 0:
            CDFKAP.append(_sp.exp(-i))
        else:
            ARG = 1-H*i
            if ARG > SMALL:
                i = -_sp.log(ARG)/H
                CDFKAP.append(_sp.exp(-i))
                continue
            else:
                CDFKAP.append(0)
                continue
    if len(CDFKAP) == 1:
        CDFKAP = CDFKAP[0]
    return(CDFKAP)
#############################################################
        
def cdfnor(x,para):
    if _is_numeric(x)==False:
        x = _sp.array(x)
        
    if para[1] < 0:
        print("Invalid Parameters")
    cdfnor = 0.5+0.5*_spsp.erf((x-para[0])/para[1]*1.0/_sp.sqrt(2))
    return(cdfnor)

#############################################################

def cdfpe3(x,para):
    SMALL = 1e-6
    
    if _is_numeric(x)==False:
        x = _sp.array(x)
        
    CDFPE3 = 0
    if para[1]<= 0:
        print("Parameters Invalid")
        return
    else:
        Gamma = para[2]
        if abs(Gamma) <= SMALL:
            Z = (x-para[0])/para[1]
            
            if type(Z) == list:
                Z = _sp.array(Z)
                
            CDFPE3 = 0.5+0.5*_spsp.erf(Z*1/_sp.sqrt(2))
            return(CDFPE3)
        else:
            Alpha = 4/(Gamma**2)
            Z = 2*(x-para[0])/(para[1]*Gamma)+Alpha
            if _is_numeric(x)==True:
                if Z > 0:
                    CDFPE3 = _spsp.gammainc(Alpha,Z)
                if Gamma < 0:
                    CDFPE3 = 1-CDFPE3
                return(CDFPE3)
            else:
                CDFPE3 = []
                for i in Z:
                    if i > 0:
                        CDFPE3.append(_spsp.gammainc(Alpha,i))
                    if Gamma < 0:
                        CDFPE3[-1] = 1-CDFPE3[-1]
                return(CDFPE3)

        
#############################################################

##THIS IS REALLY INEFFICIENT, LOOK INTO OPTIMISING TO REDUCE
##FUNCTION CALL
def cdfwak(x,para):
    if _is_numeric(x):
        return(_cdfwak2(x,para))
    else:
        save = []
        for i in x:
            save.append(_cdfwak2(i,para))
        return(save)
        
def _cdfwak2(x,para):
    
    EPS = 1e-8
    MAXIT = 20
    ZINCMX =3
    ZMULT = 0.2
    UFL = -170
    XI = para[0]
    A = para[1]
    B = para[2]
    C = para[3]
    D = para[4]

    if B+D <= 0 and (B!=0 or C!=0 or D!= 0):
        print("Invalid Parameters")
        return
    if A == 0 and B!= 0:
        print("Invalid Parameters")
        return
    if C == 0 and D != 0:
        print("Invalid Parameters")
        return
    if C < 0 or A+C < 0:
        print("Invalid Parameters")
        return
    if A == 0 and C == 0:
        print("Invalid Parameters")
        return

    CDFWAK = 0
    if x <= XI:
        return(CDFWAK)

    #Test for _special cases
    if B == 0 and C == 0 and D == 0:
        
        Z = (x-XI)/A
        CDFWAK = 1
        if -Z >= UFL:
            CDFWAK = 1-_sp.exp(-Z)
        return(CDFWAK)
    
    if C == 0:
        CDFWAK = 1
        if x >= (XI+A/B):
            return(CDFWAK)
        Z = -_sp.log(1-(x-XI)*B/A)/B
        if -Z >= UFL:
            CDFWAK = 1-_sp.exp(-Z)
        return(CDFWAK)

    
    if A == 0:
        Z = _sp.log(1+(x-XI)*D/C)/D
        if -Z >= UFL:
            CDFWAK = 1-_sp.exp(-Z)
        return(CDFWAK)


    CDFWAK=1
    if D <0 and x >= (XI+A/B-C/D):
        return(CDFWAK)

    Z=0.7
    if x < quawak(0.1,para):
        Z = 0
    if x < quawak(0.99,para):
        pass
    else:
        if D < 0:
            Z = _sp.log((x-XI-A/B)*D/C+1)/D
        if D == 0:
            Z = (x-XI-A/B)/C
        if D > 0:
            Z = _sp.log((x-XI)*D/C+1)/D
            
    for IT in range(1,MAXIT+1):
        EB = 0
        BZ = -B*Z
        if BZ >= UFL:
            EB = _sp.exp(BZ)
        GB = Z
        
        if abs(B)>EPS:
            GB = (1-EB)/B
        ED = _sp.exp(D*Z)
        GD = -Z

        if abs(D)>EPS:
            GD = (1-ED)/D

        XEST =XI +A*GB-C*GD
        FUNC = x-XEST
        DERIV1 = A*EB+C*ED
        DERIV2 = -A*B*EB+C*D*ED
        TEMP = DERIV1+0.5*FUNC*DERIV2/DERIV1

        if TEMP <= 0:
            TEMP = DERIV1
        ZINC = FUNC/TEMP

        if ZINC > ZINCMX:
            ZINC = ZINCMX

        ZNEW = Z+ZINC

        if ZNEW <= 0:
            Z = Z*ZMULT
        else:
            Z = ZNEW
            if abs(ZINC) <= EPS:
                CDFWAK = 1
                if -Z >= UFL:
                    CDFWAK = 1-_sp.exp(-Z)
                return(CDFWAK)

#############################################################   
def cdfwei(x,para):
    U = para[0]
    A = para[1]
    G = para[2]
    if _is_numeric(x)==False:
        x = _sp.array(x)
        
    if len(para) < 3:
        print("Invalid number of parameters")
        return
    elif para[1] <= 0 or para[2] <= 0:
        print("Invalid Parameters")
        return
    else:
        cdfwei = 1-_sp.exp(-((x-para[0])/para[1])**para[2])
        return(cdfwei)
    
