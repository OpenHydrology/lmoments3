import scipy as _sp
import scipy.special as _spsp
import scipy.stats as _spst
import math as _math
import sys as _sys
from _otherfunct import is_numeric as _is_numeric
import lmoments

#############################################################
###PEL FUNCTIONS
#############################################################


def pelexp(xmom):
    if xmom[1] <= 0:
        print("L-Moments Invalid")
        return
    else:
        para = [xmom[0]-2*xmom[1],2*xmom[1]]
        return(para)

#############################################################

def pelgam(xmom):
    A1 = -0.3080
    A2 = -0.05812
    A3 = 0.01765
    B1 = 0.7213
    B2 = -0.5947
    B3 = -2.1817
    B4 = 1.2113
    
    if xmom[0] <= xmom[1] or xmom[1]<= 0:
        print("L-Moments Invalid")
        return
    CV = xmom[1]/xmom[0]
    if CV >= 0.5:
        T = 1-CV
        ALPHA =T*(B1+T*B2)/(1+T*(B3+T*B4))
    else:
        T=_sp.pi*CV**2
        ALPHA=(1+A1*T)/(T*(1+T*(A2+T*A3)))
        
    para = [ALPHA,xmom[0]/ALPHA]
    return(para)

#############################################################

def pelgev(xmom):
    SMALL = 1e-5
    eps = 1e-6
    maxit = 20
    EU =0.57721566
    DL2 = _sp.log(2)
    DL3 = _sp.log(3)
    A0 =  0.28377530
    A1 = -1.21096399
    A2 = -2.50728214
    A3 = -1.13455566
    A4 = -0.07138022
    B1 =  2.06189696 
    B2 =  1.31912239 
    B3 =  0.25077104
    C1 =  1.59921491
    C2 = -0.48832213
    C3 =  0.01573152
    D1 = -0.64363929
    D2 =  0.08985247

    T3 = xmom[2]
    if xmom[1]<= 0 or abs(T3)>= 1:
        print("L-Moments Invalid")
        return
    if T3<= 0:
        G=(A0+T3*(A1+T3*(A2+T3*(A3+T3*A4))))/(1+T3*(B1+T3*(B2+T3*B3)))
        if T3>= -0.8:
            para3 = G
            GAM = _sp.exp(_sp.special.gammaln(1+G))
            para2=xmom[1]*G/(GAM*(1-2**(-G)))
            para1=xmom[0]-para2*(1-GAM)/G
            para = [para1,para2,para3]
            return(para)

        if T3 <= -0.97:
            G = 1-_sp.log(1+T3)/DL2
            
        T0=(T3+3)*0.5
        for IT in range(1,maxit):
            X2=2**(-G)
            X3=3**(-G)
            XX2=1-X2
            XX3=1-X3
            T=XX3/XX2
            DERIV=(XX2*X3*DL3-XX3*X2*DL2)/(XX2**2)
            GOLD=G
            G=G-(T-T0)/DERIV
            if abs(G-GOLD) <= eps*G:
                para3 = G
                GAM = _sp.exp(_sp.special.gammaln(1+G))
                para2=xmom[1]*G/(GAM*(1-2**(-G)))
                para1=xmom[0]-para2*(1-GAM)/G
                para = [para1,para2,para3]
                return(para)
            
        print("Iteration has not converged")

    Z=1-T3
    G=(-1+Z*(C1+Z*(C2+Z*C3)))/(1+Z*(D1+Z*D2))
    if abs(G)<SMALL:
        para2 = xmom[1]/DL2
        para1 = xmom[0]-EU*para2
        para = [para1,para2,0]
        return(para)
    else:
        para3 = G
        GAM = _sp.exp(_sp.special.gammaln(1+G))
        para2=xmom[1]*G/(GAM*(1-2**(-G)))
        para1=xmom[0]-para2*(1-GAM)/G
        para = [para1,para2,para3]
        return(para)
 
#############################################################

def pelglo(xmom):
    SMALL = 1e-6
    
    G=-xmom[2]
    if xmom[1]<= 0 or abs(G)>= 1:
        print("L-Moments Invalid")
        return

    if abs(G)<= SMALL:
        para = [xmom[0],xmom[1],0]
        return(para)

    GG = G*_sp.pi/_sp.sin(G*_sp.pi)
    A = xmom[1]/GG
    para1 = xmom[0]-A*(1-GG)/G
    para = [para1,A,G]
    return(para)

#############################################################

def pelgno(xmom):
    A0 =  0.20466534e+01
    A1 = -0.36544371e+01
    A2 =  0.18396733e+01
    A3 = -0.20360244e+00
    B1 = -0.20182173e+01
    B2 =  0.12420401e+01
    B3 = -0.21741801e+00
    SMALL = 1e-8

    T3=xmom[2]
    if xmom[1] <= 0 or abs(T3) >= 1:
        print("L-Moments Invalid")
        return
    if abs(T3)>= 0.95:
        para = [0,-1,0]
        return(para)

    if abs(T3)<= SMALL:
        para =[xmom[0],xmom[1]*_sp.sqrt(_sp.pi),0] 

    TT=T3**2
    G=-T3*(A0+TT*(A1+TT*(A2+TT*A3)))/(1+TT*(B1+TT*(B2+TT*B3)))
    E=_sp.exp(0.5*G**2)
    A=xmom[1]*G/(E*_sp.special.erf(0.5*G))
    U=xmom[0]+A*(E-1)/G
    para = [U,A,G]
    return(para)

#############################################################

def pelgpa(xmom):
    T3=xmom[2]
    if xmom[1]<= 0:
        print("L-Moments Invalid")
        return
    if abs(T3)>= 1:
        print("L-Moments Invalid")
        return

    G=(1-3*T3)/(1+T3)
    
    PARA3=G
    PARA2=(1+G)*(2+G)*xmom[1]
    PARA1=xmom[0]-PARA2/(1+G)
    para = [PARA1,PARA2,PARA3]
    return(para)

#############################################################

def pelgum(xmom):
    EU = 0.577215664901532861
    if xmom[1] <= 0:
        print("L-Moments Invalid")
        return
    else:
        para2 = xmom[1]/_sp.log(2)
        para1 = xmom[0]-EU*para2
        para = [para1, para2]
        return(para)

#############################################################

def pelkap(xmom):
    EPS = 1e-6
    MAXIT = 20
    MAXSR = 10
    HSTART = 1.001
    BIG = 10
    OFLEXP = 170
    OFLGAM = 53

    T3 = xmom[2]
    T4 = xmom[3]
    para = [0]*4
    if xmom[1] <= 0:
        print("L-Moments Invalid")
        return
    if abs(T3) >= 1 or abs(T4) >=  1:
        print("L-Moments Invalid")
        return

    if T4 <= (5*T3*T3-1)/4:
        print("L-Moments Invalid")
        return

    if T4 >= (5*T3*T3+1)/6:
        print("L-Moments Invalid")
        return

    G = (1-3*T3)/(1+T3)
    H = HSTART
    Z = G+H*0.725
    Xdist = BIG

    #Newton-Raphson Iteration
    for it in range(1,MAXIT+1):
        for i in range(1,MAXSR+1):
            if G > OFLGAM:
                print("Failed to converge")
                return
            if H > 0:
                U1 = _sp.exp(_spsp.gammaln(1/H)-_spsp.gammaln(1/H+1+G))
                U2 = _sp.exp(_spsp.gammaln(2/H)-_spsp.gammaln(2/H+1+G))
                U3 = _sp.exp(_spsp.gammaln(3/H)-_spsp.gammaln(3/H+1+G))
                U4 = _sp.exp(_spsp.gammaln(4/H)-_spsp.gammaln(4/H+1+G))
            else:
                U1 = _sp.exp(_spsp.gammaln(-1/H-G)-_spsp.gammaln(-1/H+1))
                U2 = _sp.exp(_spsp.gammaln(-2/H-G)-_spsp.gammaln(-2/H+1))
                U3 = _sp.exp(_spsp.gammaln(-3/H-G)-_spsp.gammaln(-3/H+1))
                U4 = _sp.exp(_spsp.gammaln(-4/H-G)-_spsp.gammaln(-4/H+1))

            ALAM2 =  U1-2*U2
            ALAM3 = -U1+6*U2-6*U3
            ALAM4 =  U1-12*U2+30*U3-20*U4
            if ALAM2 == 0:
                print("Failed to Converge")
                return
            TAU3 = ALAM3/ALAM2
            TAU4 = ALAM4/ALAM2
            E1 = TAU3-T3
            E2 = TAU4-T4

            DIST = max(abs(E1),abs(E2))
            if DIST < Xdist:
                Success = 1
                break
            else:
                DEL1 = 0.5*DEL1
                DEL2 = 0.5*DEL2
                G = XG-DEL1
                H = XH-DEL2
                
        if Success == 0:
            print("Failed to converge")
            return

        #Test for convergence
        if DIST < EPS:
            para[3]=H
            para[2]=G
            TEMP = _spsp.gammaln(1+G)
            if TEMP > OFLEXP:
                print("Failed to converge")
                return
            GAM = _sp.exp(TEMP)
            TEMP = (1+G)*_sp.log(abs(H))
            if  TEMP > OFLEXP:
                print("Failed to converge")
                return

            HH = _sp.exp(TEMP)
            para[1] = xmom[1]*G*HH/(ALAM2*GAM)
            para[0] = xmom[0]-para[1]/G*(1-GAM*U1/HH)
            return(para)
        else:
            XG=G
            XH=H
            XZ=Z
            Xdist=DIST
            RHH=1/(H**2)
            if H > 0:
                U1G=-U1*_spsp.psi(1/H+1+G)
                U2G=-U2*_spsp.psi(2/H+1+G)
                U3G=-U3*_spsp.psi(3/H+1+G)
                U4G=-U4*_spsp.psi(4/H+1+G)
                U1H=  RHH*(-U1G-U1*_spsp.psi(1/H))
                U2H=2*RHH*(-U2G-U2*_spsp.psi(2/H))
                U3H=3*RHH*(-U3G-U3*_spsp.psi(3/H))
                U4H=4*RHH*(-U4G-U4*_spsp.psi(4/H))
            else:
                U1G=-U1*_spsp.psi(-1/H-G)
                U2G=-U2*_spsp.psi(-2/H-G)
                U3G=-U3*_spsp.psi(-3/H-G)
                U4G=-U4*_spsp.psi(-4/H-G)
                U1H=  RHH*(-U1G-U1*_spsp.psi(-1/H+1))
                U2H=2*RHH*(-U2G-U2*_spsp.psi(-2/H+1))
                U3H=3*RHH*(-U3G-U3*_spsp.psi(-3/H+1))
                U4H=4*RHH*(-U4G-U4*_spsp.psi(-4/H+1))

            DL2G=U1G-2*U2G
            DL2H=U1H-2*U2H
            DL3G=-U1G+6*U2G-6*U3G
            DL3H=-U1H+6*U2H-6*U3H
            DL4G=U1G-12*U2G+30*U3G-20*U4G
            DL4H=U1H-12*U2H+30*U3H-20*U4H
            D11=(DL3G-TAU3*DL2G)/ALAM2
            D12=(DL3H-TAU3*DL2H)/ALAM2
            D21=(DL4G-TAU4*DL2G)/ALAM2
            D22=(DL4H-TAU4*DL2H)/ALAM2
            DET=D11*D22-D12*D21
            H11= D22/DET
            H12=-D12/DET
            H21=-D21/DET
            H22= D11/DET
            DEL1=E1*H11+E2*H12
            DEL2=E1*H21+E2*H22
            
##          TAKE NEXT N-R STEP
            G=XG-DEL1
            H=XH-DEL2
            Z=G+H*0.725

##          REDUCE STEP IF G AND H ARE OUTSIDE THE PARAMETER _spACE

            FACTOR=1
            if G <= -1:
                FACTOR = 0.8*(XG+1)/DEL1
            if H <= -1:
                FACTOR = min(FACTOR,0.8*(XH+1)/DEL2)
            if Z <= -1:
                FACTOR = min(FACTOR,0.8*(XZ+1)/(XZ-Z))
            if H <= 0 and G*H<= -1:
                FACTOR = min(FACTOR,0.8*(XG*XH+1)/(XG*XH-G*H))

            if FACTOR == 1:
                pass
            else:
                DEL1 = DEL1*FACTOR
                DEL2 = DEL2*FACTOR
                G = XG-DEL1
                H = XH-DEL2
                Z = G+H*0.725
    
#############################################################

def pelnor(xmom):
    if xmom[1] <= 0:
        print("L-Moments Invalid")
        return
    else:
        para = [xmom[0],xmom[1]*_sp.sqrt(_sp.pi)]
        return(para)

#############################################################
    
def pelpe3(xmom):
    Small = 1e-6
    #Constants used in Minimax Approx:

    C1 = 0.2906
    C2 = 0.1882
    C3 = 0.0442
    D1 = 0.36067
    D2 = -0.59567
    D3 = 0.25361
    D4 = -2.78861
    D5 = 2.56096
    D6 = -0.77045

    T3=abs(xmom[2])
    if xmom[1] <= 0 or T3 >= 1:
        para = [0]*3
        print("L-Moments Invalid")
        return(para)

    if T3<= Small:
        para = []
        para.append(xmom[0])
        para.append(xmom[1]*_sp.sqrt(_sp.pi))
        para.append(0)
        return(para)

    if T3 >= (1.0/3):
        T = 1-T3
        Alpha = T*(D1+T*(D2+T*D3))/(1+T*(D4+T*(D5+T*D6)))
    else:
        T=3*_sp.pi*T3*T3
        Alpha=(1+C1*T)/(T*(1+T*(C2+T*C3)))
                    
    RTALPH=_sp.sqrt(Alpha)
    BETA=_sp.sqrt(_sp.pi)*xmom[1]*_sp.exp(_spsp.gammaln(Alpha)-_spsp.gammaln(Alpha+0.5))
    para = []
    para.append(xmom[0])
    para.append(BETA*RTALPH)
    para.append(2/RTALPH)
    if xmom[2] < 0:
        para[2]=-para[2]

    return(para)

#############################################################

def pelwak(xmom):
    if xmom[1] <= 0:
        print("Invalid L-Moments")
        return()
    if abs(xmom[2]) >= 1 or abs(xmom[3]) >= 1 or abs(xmom[4]) >= 1:
        print("Invalid L-Moments")
        return()

    ALAM1 = xmom[0]
    ALAM2 = xmom[1]
    ALAM3 = xmom[2]*ALAM2
    ALAM4 = xmom[3]*ALAM2
    ALAM5 = xmom[4]*ALAM2

    XN1= 3*ALAM2-25*ALAM3 +32*ALAM4
    XN2=-3*ALAM2 +5*ALAM3  +8*ALAM4
    XN3= 3*ALAM2 +5*ALAM3  +2*ALAM4
    XC1= 7*ALAM2-85*ALAM3+203*ALAM4-125*ALAM5
    XC2=-7*ALAM2+25*ALAM3  +7*ALAM4 -25*ALAM5
    XC3= 7*ALAM2 +5*ALAM3  -7*ALAM4  -5*ALAM5

    XA=XN2*XC3-XC2*XN3
    XB=XN1*XC3-XC1*XN3
    XC=XN1*XC2-XC1*XN2
    DISC=XB*XB-4*XA*XC
    skip20 = 0
    if DISC < 0:
        pass
    else:       
        DISC=_sp.sqrt(DISC)
        ROOT1=0.5*(-XB+DISC)/XA
        ROOT2=0.5*(-XB-DISC)/XA
        B= max(ROOT1,ROOT2)
        D=-min(ROOT1,ROOT2)
        if D >= 1:
            pass
        else:          
            A=(1+B)*(2+B)*(3+B)/(4*(B+D))*((1+D)*ALAM2-(3-D)*ALAM3)
            C=-(1-D)*(2-D)*(3-D)/(4*(B+D))*((1-B)*ALAM2-(3+B)*ALAM3)
            XI=ALAM1-A/(1+B)-C/(1-D)
            if C >= 0 and A+C>= 0:
                skip20 = 1
                
    if skip20 == 0:
        IFAIL=1
        D=-(1-3*xmom[2])/(1+xmom[2])
        C=(1-D)*(2-D)*xmom[1]
        B=0
        A=0
        XI=xmom[0]-C/(1-D)
        if D <= 0:
            A = C
            B=-D
            C = 0
            D = 0

    para =[XI,A,B,C,D]
    return(para)


#############################################################
def pelwei(xmom):
    if len(xmom) < 3:
        print("Insufficient L-Moments: Need 3")
        return
    if xmom[1] <= 0 or xmom[2] >= 1 or xmom[2] <= -lmoments.lmrgum([0,1],3)[2]:
        print("L-Moments Invalid")
        return
    pg = pelgev([-xmom[0],xmom[1],-xmom[2]])
    delta = 1/pg[2]
    beta = pg[1]/pg[2]
    out = [-pg[0]-beta,beta,delta]
    return(out)
