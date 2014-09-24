
import scipy as _sp
import scipy.special as _spsp
import scipy.stats as _spst
import math as _math
import sys as _sys
from _otherfunct import is_numeric as _is_numeric

#############################################################
#LMR FUNCTIONS
#############################################################

def lmrexp(para,nmom):
    A=para[1]
    if A <= 0:
        print("Invalid Parameters")
        return
    if nmom > 20:
        print("Parameter nmom too large")
    xmom = []
    xmom.append(para[0]+A)
    if nmom == 1:
        return(xmom)

    xmom.append(0.5*A)
    if nmom ==2:
        return(xmom)

    for i in range(3,nmom+1):
        xmom.append(2/float(i*(i-1)))

    return(xmom)
    
#############################################################

def lmrgam(para,nmom):
    A0 = 0.32573501
    [A1,A2,A3] = [0.16869150, 0.078327243,-0.0029120539]
    [B1,B2] = [0.46697102, 0.24255406]
    C0 = 0.12260172
    [C1,C2,C3] = [0.053730130, 0.043384378, 0.011101277]
    [D1,D2]    = [0.18324466, 0.20166036]
    [E1,E2,E3] = [2.3807576, 1.5931792, 0.11618371]
    [F1,F2,F3] = [5.1533299, 7.1425260, 1.9745056]
    [G1,G2,G3] = [2.1235833, 4.1670213, 3.1925299]
    [H1,H2,H3] = [9.0551443, 26.649995, 26.193668]

    Alpha = para[0]
    Beta = para[1]
    if Alpha <= 0 or Beta <= 0:
        print("Invalid Parameters")
        return
    if nmom > 4:
        print("Parameter nmom too large")
        return
    
    xmom = []
    xmom.append(Alpha*Beta)
    if nmom == 1:
        return(xmom)

    xmom.append(Beta*1/_sp.sqrt(_sp.pi)*_sp.exp(_spsp.gammaln(Alpha+0.5)-_spsp.gammaln(Alpha)))
    if nmom == 2:
        return(xmom)

    if Alpha < 1:
        Z= Alpha
        xmom.append((((E3*Z+E2)*Z+E1)*Z+1)/(((F3*Z+F2)*Z+F1)*Z+1))
        if nmom == 3:
            return(xmom)
        xmom.append((((C3*Z+C2)*Z+C1)*Z+C0)/((D2*Z+D1)*Z+1))
        if nmom == 4:
            return(xmom)
    else:
        Z=1/Alpha
        xmom.append(_sp.sqrt(Z)*(((A3*Z+A2)*Z+A1)*Z+A0)/((B2*Z+B1)*Z+1))
        if nmom == 3:
            return(xmom)
        
        xmom.append((((C3*Z+C2)*Z+C1)*Z+C0)/((D2*Z+D1)*Z+1))
        if nmom == 4:
            return(xmom)

#############################################################

def lmrgev(para,nmom):

    ZMOM=[0.577215664901532861, 0.693147180559945309,
        0.169925001442312363,0.150374992788438185,
        0.558683500577583138e-1,0.581100239999710876e-1,
        0.276242584297309125e-1,0.305563766579053126e-1,
        0.164650282258328802e-1,0.187846624298170912e-1,
        0.109328215063027148e-1,0.126973126676329530e-1,
        0.778982818057231804e-2,0.914836179621999726e-2,
        0.583332389328363588e-2,0.690104287590348154e-2,
        0.453267970180679549e-2,0.538916811326595459e-2,
        0.362407767772390e-2,0.432387608605538096e-2]
    SMALL = 1e-6
    U = para[0]
    A = para[1]
    G = para[2]
    if A<= 0 or G <= -1:
        print("Invalid Parameters")
        return
    if nmom > 20:
        print("Parameter nmom too large")
        return

    if abs(G)>SMALL:
        GAM = _sp.exp(_spsp.gammaln(1+G))
        xmom = [U+A*(1-GAM)/G]
        if nmom == 1:
            return(xmom)

        XX2 = 1-2**(-G)
        xmom.append(A*XX2*GAM/G)
        if nmom == 2:
            return(xmom)
 
        Z0=1
        for j in range(2,nmom):
            DJ=j+1
            BETA = (1-DJ**(-G))/XX2
            Z0 = Z0*(4*DJ-6)/DJ
            Z = Z0*3*(DJ-1)/(DJ+1)
            SUM = Z0*BETA-Z
            if j == 2:
                xmom.append(SUM)
            else:
                for i in range(1,j-1):
                    DI = i+1
                    Z = Z*(DI+DI+1)*(DJ-DI)/((DI+DI-1)*(DJ+DI))
                    SUM = SUM-Z*xmom[i+1]
                xmom.append(SUM)
        return(xmom)
    
    else:
        xmom = [U]
        if nmom == 1:
            return(xmom)

        xmom.append(A*ZMOM[1])
        if nmom == 2:
            return(xmom)

        for i in range(2,nmom):
            xmom.append(zmom[i-1])

        return(xmom)
  
#############################################################
    
def lmrglo(para,nmom):
    SMALL = 1e-4
    C1 = _sp.pi**2/6
    C2 = 7*_sp.pi**4/360


    Z = [[0],[0]]
    Z.append([1])
    Z.append([0.166666666666666667,  0.833333333333333333])
    Z.append([0.416666666666666667,  0.583333333333333333])
    Z.append([0.666666666666666667e-1,  0.583333333333333333,
              0.350000000000000000])
    Z.append([0.233333333333333333,  0.583333333333333333,
              0.183333333333333333])

    Z.append([0.357142857142857143e-1,  0.420833333333333333,
              0.458333333333333333,  0.851190476190476190e-1])

    Z.append([0.150992063492063492,  0.515625000000000000,
              0.297916666666666667,  0.354662698412698413e-1])

    Z.append([0.222222222222222222e-1,  0.318893298059964727,
              0.479976851851851852,  0.165509259259259259,
              0.133983686067019400e-1])

    Z.append([0.106507936507936508,  0.447663139329805996,
              0.360810185185185185,  0.803902116402116402e-1,
              0.462852733686067019e-2])

    Z.append([0.151515151515151515e-1,  0.251316137566137566,
              0.469695216049382716,  0.227650462962962963,
              0.347139550264550265e-1,  0.147271324354657688e-2])

    Z.append([0.795695045695045695e-1,  0.389765946502057613,
              0.392917309670781893,  0.123813106261022928,
              0.134998713991769547e-1,  0.434261597456041900e-3])

    Z.append([0.109890109890109890e-1,  0.204132996632996633,
              0.447736625514403292,  0.273053442827748383,
              0.591917438271604938e-1,  0.477687757201646091e-2,
              0.119302636663747775e-3])

    Z.append([0.619345205059490774e-1,  0.342031759392870504,
              0.407013705173427396,  0.162189192806752331,
              0.252492100235155791e-1,  0.155093427662872107e-2,
              0.306778208563922850e-4])

    Z.append([0.833333333333333333e-2,  0.169768364902293474,
              0.422191282868366202,  0.305427172894620811,
              0.840827939972285210e-1,  0.972435791446208113e-2,
              0.465280282988616322e-3,  0.741380670696146887e-5])

    Z.append([0.497166028416028416e-1,  0.302765838589871328,
              0.410473300089185506,  0.194839026503251764,
              0.386598063704648526e-1,  0.341399407642897226e-2,
              0.129741617371825705e-3,  0.168991182291033482e-5])
             
    Z.append([0.653594771241830065e-2,  0.143874847595085690,
              0.396432853710259464,  0.328084180720899471,
              0.107971393165194318,  0.159653369932077769e-1,
              0.110127737569143819e-2,  0.337982364582066963e-4,
              0.364490785333601627e-6])

    Z.append([0.408784570549276431e-1,  0.270244290725441519,
              0.407599524514551521,  0.222111426489320008,
              0.528463884629533398e-1,  0.598298239272872761e-2,
              0.328593965565898436e-3,  0.826179113422830354e-5,
              0.746033771150646605e-7])

    Z.append([0.526315789473684211e-2,  0.123817655753054913,
              0.371859291444794917,  0.343568747670189607,
              0.130198662812524058,  0.231474364899477023e-1,
              0.205192519479869981e-2,  0.912058258107571930e-4,
              0.190238611643414884e-5,  0.145280260697757497e-7])

    U = para[0]
    A = para[1]
    G = para[2]

    if A <= 0 or abs(G) >= 1:
        print("Invalid Parameters")
        return

    if nmom > 20:
        print("Parameter nmom too large")
        return
    GG = G*G
    ALAM1 = -G*(C1+GG*C2)
    ALAM2 = 1+GG*(C1+GG*C2)
    if abs(G) > SMALL:
        ALAM2=G*_sp.pi/_sp.sin(G*_sp.pi)
        ALAM1=(1-ALAM2)/G

    xmom = [U+A*ALAM1]
    if nmom == 1:
        return(xmom)
             
    xmom.append(A*ALAM2)
    if nmom == 2:
        return(xmom)

    for M in range(3,nmom+1):
        kmax = M/2
        SUMM=Z[M-1][kmax-1]
        for K in range(kmax-1,0,-1):
            SUMM = SUMM*GG+Z[M-1][K-1]
        if M != M/2*2:
            SUMM = -G*SUMM
        xmom.append(SUMM)

    return(xmom)

#############################################################

def lmrgno(para,nmom):

    ZMOM = [0,   0.564189583547756287, 0,   0.122601719540890947,
            0,   0.436611538950024944e-1,0, 0.218431360332508776e-1,
            0,   0.129635015801507746e-1,0, 0.852962124191705402e-2,
            0,   0.601389015179323333e-2,0, 0.445558258647650150e-2,
            0,   0.342643243578076985e-2,0, 0.271267963048139365e-2]


    RRT2 = 1/_sp.sqrt(2)
    RRTPI = 1/_sp.sqrt(_sp.pi)
    
    RANGE = 5
    EPS = 1e-8
    MAXIT = 10

    U = para[0]
    A = para[1]
    G = para[2]
    if A <= 0:
        print("Invalid Parameters")
        return
    if nmom > 20:
        print("Parameter nmom too large")
        return

    

    if abs(G)<=EPS:
        xmom = [U]
        if nmom == 1:
            return(xmom)

        xmom.append(A*ZMOM[1])
        if nmom == 2:
            return(xmom)

        for i in range(3,nmom+1):
            xmom.append(zmom[i-1])

        return(xmom)


    EGG = _sp.exp(0.5*G**2)
    ALAM1 = (1-EGG)/G
    xmom = [U+A*ALAM1]
    if nmom == 1:
        return(xmom)
    
    ALAM2=EGG*_spsp.erf(0.5*G)/G
    xmom.append(A*ALAM2)
    if nmom == 2:
        return(xmom)
  
    CC=-G*RRT2
    XMIN=CC-RANGE
    XMAX=CC+RANGE
    SUMM = [0]*nmom
    
    N=16
    XINC=(XMAX-XMIN)/N

    for i in range(1,N):
        X = XMIN+i*XINC
        E = _sp.exp(-((X-CC)**2))
        D = _spsp.erf(X)
        P1 = 1
        P = D
        for M in range(3,nmom+1):
            C1=M+M-3
            C2=M-2
            C3=M-1
            P2=P1
            P1=P
            P=(C1*D*P1-C2*P2)/C3
            SUMM[M-1] = SUMM[M-1]+E*P

    EST = []
    for i in SUMM:
        EST.append(i*XINC)


    for IT in range(1,MAXIT+1):
        ESTX = EST
        N=N*2
        XINC=(XMAX-XMIN)/N
        for i in range(1,N-1,2):
            X = XMIN+i*XINC
            E = _sp.exp(-((X-CC)**2))
            D = _spsp.erf(X)
            P1 = 1
            P = D
            for M in range(3,nmom+1):
                C1=M+M-3
                C2=M-2
                C3=M-1
                P2=P1
                P1=P
                P=(C1*D*P1-C2*P2)/C3
                SUMM[M-1] = SUMM[M-1]+E*P

        NOTCGD = 0
        for M in range(nmom,2,-1):
            EST[M-1] = SUMM[M-1]*XINC
            if abs(EST[M-1]-ESTX[M-1]) > EPS*abs(EST[M-1]):
                NOTCGD = M

        if NOTCGD == 0:
            CONST = -_sp.exp(CC**2)*RRTPI/(ALAM2*G)
            
            for M in range(3,nmom+1):
                xmom.append(CONST*EST[M-1])
            return(xmom)
        else:
            print("Did Not Converge")
            return
        
#############################################################
                
def lmrgpa(para,nmom):
    U = para[0]
    A = para[1]
    G = para[2]
    if A <=0 or G < -1:
        print("Invalid Parameters")
        return
    if nmom > 20:
        print("Parameter nmom too large")
        return

    Y = 1/(1+G)
    xmom = [U+A*Y]
    if nmom == 1:
        return(xmom)
    
    Y = Y/(2+G)
    xmom.append(A*Y)
    if nmom == 2:
        return(xmom)
    
    Y = 1
    for i in range(3,nmom+1):
        AM = i-2
        Y = Y*(AM-G)/(i+G)
        xmom.append(Y)
    return(xmom)

#############################################################

def lmrgum(para,nmom):
    ZMOM = [0.577215664901532861,  0.693147180559945309,
     0.169925001442312363,  0.150374992788438185,
     0.0558683500577583138,  0.0581100239999710876,
     0.0276242584297309125,  0.0305563766579053126,
     0.0164650282258328802,  0.0187846624298170912,
     0.0109328215063027148,  0.0126973126676329530,
     0.00778982818057231804,  0.00914836179621999726,
     0.00583332389328363588,  0.00690104287590348154,
     0.00453267970180679549,  0.00538916811326595459,
     0.00362407767772368790,  0.00432387608605538096]

    A = para[1]
    if A <=0:
        print("Invalid Parameters")
        return
    if nmom >20:
        print("Parameter nmom too large")
        return
    xmom = [para[0]+A*ZMOM[0]]
    if nmom == 1:
        return(xmom)
    xmom.append(A*ZMOM[1])
    if nmom == 2:
        return(xmom)

    for i in range(2,nmom):
        xmom.append(ZMOM[i])
    return(xmom)

#############################################################

def lmrkap(para,nmom):
    EU = 0.577215664901532861
    SMALL = 1e-8
    OFL = 170
    U = para[0]
    A = para[1]
    G = para[2]
    H = para[3]

    if A <= 0 or G <= -1: 
        print("Invalid Parameters")
        return
    if H < 0 and (G*H)<= -1:
        print("Invalid Parameters")
        return
    if nmom > 20:
        print("Parameter nmom too large")
        return

    DLGAM = _spsp.gammaln(1+G)
    ICASE = 1
    if H > 0:
        ICASE = 3
    elif abs(H) < SMALL:
        ICASE = 2
    elif G == 0:
        ICASE = ICASE+3

    if ICASE == 1:
        Beta = []
        for IR in range(1,nmom+1):
            ARG = DLGAM + _spsp.gammaln(-IR/H-G) - _spsp.gammaln(-IR/H)-G*_sp.log(-H)
            if abs(ARG) > OFL:
                print("Calculation of L-Moments Failed")
                return
            Beta.append(_sp.exp(ARG))


    elif ICASE == 2:
        Beta = []
        for IR in range(1,nmom+1):
            Beta.append(_sp.exp(DLGAM-G*_sp.log(IR))*(1-0.5*H*G*(1+G)/IR))

    elif ICASE == 3:
        Beta = []
        for IR in range(1,nmom+1):
            ARG = DLGAM+ _spsp.gammaln(1+IR/H)-_spsp.gammaln(1+G+IR/H)-G*_sp.log(H)
            if abs(ARG) > OFL:
                print("Calculation of L-Moments Failed")
                return
            Beta.append(_sp.exp(ARG))
            
    elif ICASE == 4:
        Beta = []
        for IR in range(1,nmom+1):
            Beta.append(EU+_sp.log(-H)+_spsp.psi(-IR/H))
            
    elif ICASE == 5:
        Beta = []
        for IR in range(1,nmom+1):
            Beta.append(EU+_sp.log(IR))

    elif ICASE == 6:
        Beta = []
        for IR in range(1,nmom+1):
            Beta.append(EU+_sp.log(H)+_spsp.psi(1+IR/H))

    if G == 0:
        xmom = [U+A*Beta[0]]
    else:
        xmom = [U+A*(1-Beta[0])/G]

    if nmom == 1:
        return(xmom)

    ALAM2 = Beta[1]-Beta[0]
    if G == 0:
        xmom.append(A*ALAM2)
    else:
        xmom.append(A*ALAM2/(-G))

    if nmom == 2:
        return(xmom)

    Z0 = 1
    for j in range(3,nmom+1):
        Z0 = Z0*(4.0*j-6)/j
        Z = 3*Z0*(j-1)/(j+1)
        SUMM = Z0*(Beta[j-1]-Beta[0])/ALAM2 - Z
        if j == 3:
            xmom.append(SUMM)
        else:
            for i in range(2,j-1):
                Z = Z*(i+i+1)*(j-i)/((i+i-1)*(j+i))
                SUMM = SUMM - Z*xmom[i]
            xmom.append(SUMM)
    return(xmom)

#############################################################

def lmrnor(para,nmom):

    ZMOM =[0, 0.564189583547756287, 0,   0.122601719540890947,
        0,  0.0436611538950024944, 0,   0.0218431360332508776,
        0,  0.0129635015801507746, 0,   0.00852962124191705402,
        0,  0.00601389015179323333, 0,   0.00445558258647650150,
        0,  0.00342643243578076985, 0,   0.00271267963048139365]

    if para[1] <= 0:
        print("Invalid Parameters")
        return
    if nmom > 20:
        print("Parameter nmom too large")
        return
    xmom = [para[0]]
    if nmom == 1:
        return(xmom)

    xmom.append(para[1]*ZMOM[1])
    if nmom == 2:
        return(xmom)

    for M in range(2,nmom):
        xmom.append(ZMOM[M])

    return(xmom)

#############################################################

def lmrpe3(para,nmom):
    SMALL = 1e-6
    CONST = 1/_sp.sqrt(_sp.pi)
    A0 = 0.32573501
    [A1,A2,A3] = [0.16869150, 0.078327243,-0.0029120539]
    [B1,B2] = [0.46697102,0.24255406]
    C0 = 0.12260172
    [C1,C2,C3] = 0.053730130, 0.043384378, 0.011101277
    [D1,D2] = [0.18324466, 0.20166036]
    [E1,E2,E3] = [2.3807576, 1.5931792, 0.11618371]
    [F1,F2,F3] = [5.1533299, 7.1425260, 1.9745056]
    [G1,G2,G3] = [2.1235833, 4.1670213, 3.1925299]
    [H1,H2,H3] = [9.0551443, 26.649995, 26.193668]

    SD = para[1]
    if SD <= 0:
        print("Invalid Parameters")
        return
    if nmom > 4:
        print("Parameter nmom too large")
        return

    xmom = [para[0]]
    if nmom == 1:
        return(xmom)

    Gamma = para[2]
    if abs(Gamma) < SMALL:
        xmom = [para[0]]
        if nmom == 1:
            return(xmom)

        xmom.append(CONST*Para[1])
        if nmom == 2:
            return(xmom)

        xmom.append(0)
        if nmom == 3:
            return(xmom)

        xmom.append(C0)
        return(xmom)
    else:
        Alpha = 4/(Gamma*Gamma)
        Beta = abs(0.5*SD*Gamma)
        ALAM2 = CONST*_sp.exp(_spsp.gammaln(Alpha+0.5)-_spsp.gammaln(Alpha))
        xmom.append(ALAM2*Beta)
        if nmom == 2:
            return(xmom)

        if Alpha < 1:
            Z = Alpha
            xmom.append((((E3*Z+E2)*Z+E1)*Z+1)/(((F3*Z+F2)*Z+F1)*Z+1))
            if Gamma<0:
                xmom[2] = -xmom[2]
            if nmom == 3:
                return(xmom)

            xmom.append((((G3*Z+G2)*Z+G1)*Z+1)/(((H3*Z+H2)*Z+H1)*Z+1))
            return(xmom)

        else:
            Z = 1.0/Alpha
            xmom.append(_sp.sqrt(Z)*(((A3*Z+A2)*Z+A1)*Z+A0)/((B2*Z+B1)*Z+1))
            if Gamma < 0:
                xmom[2] = -xmom[2]

            if nmom == 3:
                return(xmom)

            xmom.append((((C3*Z+C2)*Z+C1)*Z+C0)/((D2*Z+D1)*Z+1))
            return(xmom)


#############################################################

def lmrwak(para,nmom):
    [XI,A,B,C,D]=para
    fail = 0
    if D >= 1:
        fail = 1
    if (B+D)<= 0 and (B!= 0 or C != 0 or D!=0):
        fail = 1
    if A == 0 and B != 0:
        fail = 1
    if C == 0 and D != 0:
        fail = 1
    if C < 0:
        fail = 1
    if (A+C) < 0:
        fail = 1
    if A == 0 and C == 0:
        fail = 1
    if nmom >= 20:
        fail = 2

    if fail == 1:
        print("Invalid Parameters")
        return
    if fail == 2:
        print("Parameter nmom too large")
        return

    Y=A/(1+B)
    Z=C/(1-D)
    xmom = []
    xmom.append(XI+Y+Z)
    if nmom == 1:
        return
    
    Y=Y/(2+B)
    Z=Z/(2-D)
    ALAM2=Y+Z
    xmom.append(ALAM2)
    if nmom == 2:
        return

    for i in range(2,nmom):
        AM=i+1
        Y=Y*(AM-2-B)/(AM+B)
        Z=Z*(AM-2+D)/(AM-D)
        xmom.append((Y+Z)/ALAM2)

    return(xmom)

#############################################################
def lmrwei(para,nmom):
    if len(para) != 3:
        print("Invalid number of parameters")
        return
    if para[1] <= 0 or para[2] <= 0:
        print("Invalid Parameters")
        return
    
    xmom = lmrgev([0,para[1]/para[2],1/para[2]],nmom)
    xmom[0] = para[0]+para[1] - xmom[0]
    xmom[2] = -xmom[2]
    return(xmom)

