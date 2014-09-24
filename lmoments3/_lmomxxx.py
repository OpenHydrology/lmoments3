import scipy as _sp
import scipy.special as _spsp
import scipy.stats as _spst
import math as _math
import sys as _sys
from _otherfunct import is_numeric as _is_numeric

##############################################################
#LMOM FUNCTIONS
##############################################################
def lmomexp(para):
    A = para[1]
    L1 = para[0]+A
    L2 = 0.5*A
    TAU3 = 2.0/6
    TAU4 = 2.0/12
    TAU5 = 2.0/20
    L3 = TAU3*L2
    L4 = TAU4*L2
    L5 = TAU5*L2

    return([L1,L2,L3,L4,L5])

##############################################################
def lmomgam(para):
    CONST = 0.564189583547756
    A0 = 0.32573501
    A1 = 0.1686915
    A2 = 0.078327243
    A3 = -0.0029120539
    B1 = 0.46697102
    B2 = 0.24255406
    C0 = 0.12260172
    C1 = 0.05373013
    C2 = 0.043384378
    C3 = 0.011101277
    D1 = 0.18324466
    D2 = 0.20166036
    E1 = 2.3807576
    E2 = 1.5931792
    E3 = 0.11618371
    F1 = 5.1533299
    F2 = 7.142526
    F3 = 1.9745056
    G1 = 2.1235833
    G2 = 4.1670213
    G3 = 3.1925299
    H1 = 9.0551443
    H2 = 26.649995
    H3 = 26.193668
    ALPHA = para[0]
    BETA = para[1]
    L1 = ALPHA*BETA
    L2 = BETA*CONST*_sp.exp(_spsp.gammaln(ALPHA+0.5)-_spsp.gammaln(ALPHA))
    if ALPHA < 1:
        Z = ALPHA
        TAU3 = (((E3*Z+E2)*Z+E1)*Z+1)/(((F3*Z+F2)*Z+F1)*Z+1)
        TAU4 = (((G3*Z+G2)*Z+G1)*Z+1)/(((H3*Z+H2)*Z+H1)*Z+1)
    else:
        Z = 1/ALPHA
        TAU3 = _sp.sqrt(Z)*(((A3*Z+A2)*Z+A1)*Z+A0)/((B2*Z+B1)*Z+1)
        TAU4 = (((C3*Z+C2)*Z+C1)*Z+C0)/((D2*Z+D1)*Z+1)
    L3 = TAU3 * L2
    L4 = TAU4 * L2

    return([L1,L2,L3,L4])



##############################################################
def lmomgev(para):
    
    ZMOM = [0.577215664901533, 0.693147180559945, 0.169925001442312, 
        0.150374992788438, 0.0558683500577583]
    SMALL = 1e-06

    [U,A,G] = para
    if abs(G) <= SMALL:
        L1 = U
        L2 = A*ZMOM[1]
        L3 = L2*ZMOM[2]
        L4 = L2*ZMOM[3]
        L5 = L2*ZMOM[4]
        return([L1,L2,L3,L4,L5])
   
    else:
        GAM = _sp.exp(_spsp.gammaln(1 + G))
        L1 = U+A*(1-GAM)/G
        XX2 = 1-2**(-G)
        L2 = A * XX2 * GAM/G
        Z0 = 1
        save = []
        for J in range(3,6):
            BETA = (1-J**(-G))/XX2
            Z0 = Z0*(4*J-6)/J
            Z = Z0*3*(J-1)/(J+1)
            SUM = Z0*BETA - Z
            if J > 3:
                for I in range(2,J-1):
                    Z = Z*(I+I+1)*(J-I)/((I+I-1)*(J+I))
                    SUM = SUM - Z *save[I-2]

            save.append(SUM)

        L3 = save[0]*L2
        L4 = save[1]*L2
        L5 = save[2]*L2
        return([L1,L2,L3,L4,L5])
            
##############################################################
def lmomglo(para):
    SMALL = 1e-04
    C1 = 1.64493406684823
    C2 = 1.89406565899449
    [XI,A,K] = para
    KK = K*K
    ALAM1 = -K * (C1 + KK * C2)
    ALAM2 = 1 + KK * (C1 + KK * C2)
    if abs(K) > SMALL:
        ALAM2 = K * _sp.pi/_sp.sin(K*_sp.pi)
        ALAM1 = (1 - ALAM2)/K
    L1 = XI + A * ALAM1
    L2 = A * ALAM2
    TAU3 = -K
    TAU4 = (1 + 5 * K**2)/6
    LCV = L2/L1
    L3 = TAU3 * L2
    L4 = TAU4 * L2
    return([L1,L2,L3,L4])
    
##############################################################
def lmomgno(para):

    SUM = [0,0,0,0,0]
    EST = [0,0,0,0,0]
    ESTX = [0,0,0,0,0]
    ZMOM = [0, 0.564189583547756, 0, 0.122601719540891, 0]
    RRT2 = 1/_sp.sqrt(2)
    RRTPI = 1/_sp.sqrt(_sp.pi)
    RANGE = 5
    EPS = 1e-08
    MAXIT = 10
    [U,A,G] = para
    if abs(G) <= EPS:
        L1 = U
        L2 = A * ZMOM[1]
        TAU3 = ZMOM[2]
        TAU4 = ZMOM[3]
        TAU5 = ZMOM[4]
        LCV = 2/L1
        L3 = TAU3 * L2
        L4 = TAU4 * L2
        L5 = TAU5 * L2
        return([L1,L2,L3,L4,L5])
    
    EGG = _sp.exp(0.5 * G**2)
    ALAM1 = (1 - EGG)/G
    L1 = U + A * ALAM1
    ALAM2 = EGG * _spsp.erf(0.5 * G)/G
    L2 = A * ALAM2
    CC = -G * RRT2
    XMIN = CC - RANGE
    XMAX = CC + RANGE
    N = 16
    XINC = (XMAX - XMIN)/N
    for i in range(1,N):
        X = XMIN + i * XINC
        E = _sp.exp(-((X - CC)**2))
        D = _spsp.erf(X)
        P1 = 1
        P = D
        for m in range(3,6):
            C1 = m + m - 3
            C2 = m - 2
            C3 = m - 1
            P2 = P1
            P1 = P
            P = (C1 * D * P1 - C2 * P2)/C3
            SUM[m-1] = SUM[m-1] + E * P

    EST[2] = SUM[2] * XINC
    EST[3] = SUM[3] * XINC
    EST[4] = SUM[4] * XINC
    for it in range(1, MAXIT+1):
        ESTX[2] = EST[2]
        ESTX[3] = EST[3]
        ESTX[4] = EST[4]
        N = N * 2
        XINC = (XMAX - XMIN)/N
        for i in range(1, N, 2):
            X = XMIN + i * XINC
            E = _sp.exp(-((X - CC)**2))
            D = _spsp.erf(X)
            P1 = 1
            P = D
            for m in range(3, 6):
                C1 = m + m - 3
                C2 = m - 2
                C3 = m - 1
                P2 = P1
                P1 = P
                P = (C1 * D * P1 - C2 * P2)/C3
                SUM[m-1] = SUM[m-1] + E * P

        NOTCGD = 0
        for m in range(5, 2, -1):
            EST[m-1] = SUM[m-1] * XINC
            if abs(EST[m-1] - ESTX[m-1]) > EPS * abs(EST[m-1]):
                NOTCGD = m
        
        if NOTCGD == 0:
            break
    
    if NOTCGD != 0:
        print("ITERATION HAS NOT CONVERGED. ONLY THE FIRST "+str(NOTCGD - 1)+" L-MOMENTS ARE RELIABLE")
    
    CONST = -_sp.exp(CC * CC) * RRTPI/(ALAM2 * G)
    TAU3 = CONST * EST[2]
    TAU4 = CONST * EST[3]
    TAU5 = CONST * EST[4]
    LCV = L2/L1
    L3 = TAU3 * L2
    L4 = TAU4 * L2
    L5 = TAU5 * L2
    return([L1,L2,L3,L4,L5])
    
##############################################################
def lmomgpa(para):
    [XI,A,K] = para
    Y = 1/(1 + K)
    L1 = XI + A*Y
    Y = Y/(2+K)
    L2 = A*Y
    x = [0,0,0,0,0]
    Y = 1
    for m in range(3, 6):
        AM = m - 2
        Y = Y * (AM - K)/(m + K)
        x[m-1] = Y
    
    TAU3 = x[2]
    TAU4 = x[3]
    TAU5 = x[4]
    LCV = L2/L1
    L3 = TAU3 * L2
    L4 = TAU4 * L2
    L5 = TAU5 * L2
    return([L1,L2,L3,L4,L5])

    
##############################################################
def lmomgum(para):
    ZMOM = [0.577215664901533, 0.693147180559945, 0.169925001442312, 
        0.150374992788438, 0.0558683500577583]
    A = para[1]
    L1 = para[0] + A * ZMOM[0]
    L2 = A * ZMOM[1]
    TAU3 = ZMOM[2]
    TAU4 = ZMOM[3]
    TAU5 = ZMOM[4]
    LCV = L2/L1
    L3 = TAU3*L2
    L4 = TAU4*L2
    L5 = TAU5*L2
    return([L1,L2,L3,L4,L5])
    
##############################################################
def lmomkap(para):

    def kapicase1(U,A,G,H):
        BETA = [0,0,0,0,0]
        OFL = _sp.log(_sys.float_info.max)
        DLGAM = _spsp.gammaln(1+G)
        for R in range(1,6):
            ARG = DLGAM +_spsp.gammaln(-R/H-G)-_spsp.gammaln(-R/H)-G*_sp.log(-H)
            if abs(ARG) > OFL:
                print("Calculations of L-Moments have Failed")
                return()
            BETA[R-1] = _sp.exp(ARG)
        return(BETA)

    def kapicase2(U,A,G,H):
        BETA = [0,0,0,0,0]
        DLGAM = _spsp.gammaln(1+G)
        for R in range(1,6):
            BETA[R-1] = _sp.exp(DLGAM-G*_sp.log(R))*(1-0.5*H*G*(1+G)/R)
        return(BETA)

    def kapicase3(U,A,G,H):
        BETA = [0,0,0,0,0]
        OFL = _sp.log(_sys.float_info.max)
        DLGAM = _spsp.gammaln(1+G)
        for R in range(1, 6):
            ARG = DLGAM + _spsp.gammaln(1+R/H)-_spsp.gammaln(1+G+R/H)-G*_sp.log(H)
            if abs(ARG) > OFL:
                print("Calculations of L-moments have broken down")
                return()
            BETA[R-1] = _sp.exp(ARG)
        return(BETA)
    

    def kapicase4(U,A,G,H):
        BETA = [0,0,0,0,0]
        EU = 0.577215664901533
        for R in range(1,6): 
            BETA[R-1] = EU + _sp.log(-H) + _spsp.psi(-R/H)
        return(BETA)

    def kapicase5(U,A,G,H):
        BETA = [0,0,0,0,0]
        EU = 0.577215664901533
        for R in range(1,6):
            BETA[R-1] = EU + _sp.log(R)
        return(BETA)

    def kapicase6(U,A,G,H):
        BETA = [0,0,0,0,0]
        EU = 0.577215664901533
        for R in range(1,6):
            BETA[R-1] = EU + _sp.log(H)+_spsp.psi(1+R/H)
        return(BETA)
    
    SMALL = 1e-08
    [U,A,G,H] = para
    ICASE = 1
    if H > 0:
        ICASE = 3
    if abs(H) < SMALL:
        ICASE = 2
    if G == 0:
        ICASE = ICASE + 3
    if ICASE == 1:
        BETA = kapicase1(U, A, G, H)
    if ICASE == 2:
        BETA = kapicase2(U, A, G, H)
    if ICASE == 3:
        BETA = kapicase3(U, A, G, H)
    if ICASE == 4:
        BETA = kapicase4(U, A, G, H)
    if ICASE == 5:
        BETA = kapicase5(U, A, G, H)
    if ICASE == 6:
        BETA = kapicase6(U, A, G, H)
    if G == 0:
        L1 = U+A*BETA[0]
    else:
        L1 = U+A*(1 - BETA[0])/G
    ALAM2 = BETA[1] - BETA[0]
    if G == 0:
        L2 = A*ALAM2
    else:
        L2 = A*ALAM2/(-G)
    
    LCV = L2/L1
    Z0 = 1
    x = [0,0,0,0,0]
    for J in range(3, 6):
        Z0 = Z0 * (4 * J - 6)/J
        Z = Z0 * 3 * (J - 1)/(J + 1)
        SUM = Z0 * (BETA[J-1] - BETA[0])/ALAM2 - Z
        if J == 3:
            x[J-1] = SUM
        else:
            for I in range(2, J - 1):
                Z = Z*(I+I+1)*(J-I)/((I+I-1)*(J+I))
                SUM = SUM - Z * x[I]
            
            x[J-1] = SUM
        
    
    TAU3 = x[2]
    TAU4 = x[3]
    TAU5 = x[4]
    L3 = TAU3 * LCV
    L4 = TAU4 * LCV
    L5 = TAU5 * LCV
    return([L1,L2,L3,L4,L5])

##############################################################
def lmomnor(para):
    ZMOM = [0, 0.564189583547756, 0, 0.122601719540891, 0]
    RRT2 = 1/_sp.sqrt(2)
    RRTPI = 1/_sp.sqrt(_sp.pi)
    L1 = para[0]
    L2 = para[1] * ZMOM[1]
    TAU3 = ZMOM[2]
    TAU4 = ZMOM[3]
    TAU5 = ZMOM[4]
    LCV = L2/L1
    L3 = TAU3 * L2
    L4 = TAU4 * L2
    L5 = TAU5 * L2
    return([L1,L2,L3,L4,L5])

##############################################################
def lmompe3(para):
    SMALL = 1e-06
    CONST = 1/_sp.sqrt(_sp.pi)
    A0 = 1/_sp.sqrt(3*_sp.pi)
    A1 = 0.1686915
    A2 = 0.078327243
    A3 = -0.0029120539
    B1 = 0.46697102
    B2 = 0.24255406
    C0 = 0.12260172
    C1 = 0.05373013
    C2 = 0.043384378
    C3 = 0.011101277
    D1 = 0.18324466
    D2 = 0.20166036
    E1 = 2.3807576
    E2 = 1.5931792
    E3 = 0.11618371
    F1 = 5.1533299
    F2 = 7.142526
    F3 = 1.9745056
    G1 = 2.1235833
    G2 = 4.1670213
    G3 = 3.1925299
    H1 = 9.0551443
    H2 = 26.649995
    H3 = 26.193668
    SD = para[1]
    L1 = para[0]
    GAMMA = para[2]
    if abs(GAMMA) < SMALL:
        L2 = CONST * para[1]
        TAU3 = 0
        TAU4 = C0
        L3 = L2 * TAU3
        L4 = L2 * TAU4
    
    else:
        ALPHA = 4/(GAMMA*GAMMA)
        BETA = abs(0.5*SD*GAMMA)
        ALAM2 = CONST*_sp.exp(_spsp.gammaln(ALPHA+0.5) - _spsp.gammaln(ALPHA))
        L2 = ALAM2*BETA
        if ALPHA < 1:
            Z = ALPHA
            TAU3 = (((E3 * Z + E2) * Z + E1) * Z + 1)/(((F3 * 
                Z + F2) * Z + F1) * Z + 1)
            if GAMMA < 0:
                TAU3 = -TAU3
            TAU4 = (((G3 * Z + G2) * Z + G1) * Z + 1)/(((H3 * 
                Z + H2) * Z + H1) * Z + 1)
            L3 = L2 * TAU3
            L4 = L2 * TAU4
        
        else:
            Z = 1/ALPHA
            TAU3 = _sp.sqrt(Z) * (((A3 * Z + A2) * Z + A1) * Z + 
                A0)/((B2 * Z + B1) * Z + 1)
            if GAMMA < 0:
                TAU3 = -TAU3
            TAU4 = (((C3 * Z + C2) * Z + C1) * Z + C0)/((D2 * 
                Z + D1) * Z + 1)
            L3 = L2 * TAU3
            L4 = L2 * TAU4


    LCV = L2/L1
    return([L1,L2,L3,L4])


##############################################################
def lmomwak(para):
    [XI,A,B,C,D] = para

    Y = A/(1+B)
    Z = C/(1-D)
    L1 = XI+Y+Z
    Y = Y/(2 + B)
    Z = Z/(2 - D)
    ALAM2 = Y + Z
    L2 = ALAM2
    x = [0,0,0,0,0]
    for M in range(3, 6):
        Y = Y * (M - 2 - B)/(M + B)
        Z = Z * (M - 2 + D)/(M - D)
        x[M-1] = (Y + Z)/ALAM2

    TAU3 = x[2]
    TAU4 = x[3]
    TAU5 = x[4]
    LCV = L2/L1
    L3 = TAU3 * L2
    L4 = TAU4 * L2
    L5 = TAU5 * L2
    return([L1,L2,L3,L4,L5])

##############################################################
def lmomwei(para):
    [ZETA,B,D] = para
    K = 1/D
    A = B/D
    XI = ZETA - B
    z = lmomgev([XI,A,K])
    z[0] = -z[0]
    z[2] = -z[2]
    z[4] = -z[4]
    return(z)
