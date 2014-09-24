import scipy as _sp
import scipy.special as _spsp
import scipy.stats as _spst
import math as _math
import sys as _sys
from _otherfunct import is_numeric as _is_numeric

#############################################################
##QUANTILE FUNCTIONS
#############################################################

def quaexp(F,para):
    U = para[0]
    A = para[1]
    if A <= 0:
        print("Parameters Invalid")
        return

    if _is_numeric(F):
        if F <= 0 or F >= 1:
            print("F Value Invalid")
            return
        QUAEXP = U-A*_sp.log(1-F)
        return(QUAEXP)

    else:
        if any([i >= 1 or i<= 0 for i in F]):
            print("F Value Invalid")
            return
        F = _sp.array(F)
        
    QUAEXP = U-A*_sp.log(1-F)
    return(QUAEXP)

#############################################################

def quagam(F,para):
    EPS = 1e-10
    maxit = 30
    QUAGAM = 0
    Alpha = para[0]
    Beta = para[1]
    if Alpha <= 0 or Beta <= 0:
        print("Parameters Invalid")
        return

    if _is_numeric(F):
        if F<=0 or F>= 1:
            print("F Value Invalid")
            return
    else:
        if any([i >= 1 or i<= 0 for i in F]):
            print("F Value Invalid")
            return
        F = _sp.array(F)

    if _is_numeric(F):
        AM1 = Alpha - 1
        if AM1 != 0:
            DLOGG = _spsp.gammaln(Alpha)
            if AM1 <= 0:
                Root = _sp.exp((_sp.log(Alpha*F)+DLOGG)/Alpha)
            else:
                Root = Alpha*(1-1/(9*Alpha) + quastn(F)/_sp.sqrt(9*Alpha))**3

            if Root <= 0.01*Alpha:
                Root = _sp.exp((_sp.log(Alpha*F)+DLOGG)/Alpha)

            for it in range(1,maxit+1):
                FUNC = _spsp.gammainc(Alpha,Root)-F
                RINC = FUNC*_sp.exp(DLOGG+Root-AM1*_sp.log(Root))
                Root = Root-RINC
                if abs(FUNC) <= EPS:
                    QUAGAM = Root*Beta
                    return(QUAGAM)

            if it == maxit:
                print("Failed to converge QUAGAM")
                return()
                
        else:
            QUAGAM = -_sp.log(1-F)*Beta
            return(QUAGAM)

        

    else:
        QUAGAM = []
        for Fiterate in F:
            AM1 = Alpha - 1
            if AM1 != 0:
                
                DLOGG = _spsp.gammaln(Alpha)
                
                if AM1 <= 0:
                    Root = _sp.exp((_sp.log(Alpha*Fiterate)+DLOGG)/Alpha)
                else:
                    Root = Alpha*(1-1/(9*Alpha) + quastn(Fiterate)/_sp.sqrt(9*Alpha))**3

                if Root <= 0.01*Alpha:
                    Root = _sp.exp((_sp.log(Alpha*Fiterate)+DLOGG)/Alpha)

                for it in range(1,maxit+1):
                    FUNC = _spsp.gammainc(Alpha,Root)-Fiterate
                    RINC = FUNC*_sp.exp(DLOGG+Root-AM1*_sp.log(Root))
                    Root = Root-RINC
                    if _is_numeric(FUNC):
                        if abs(FUNC) <= EPS:
                            QUAGAM.append(Root*Beta)
                            break


                if it == maxit:
                    print("Failed to converge QUAGAM")
                    return()
                
            else:
                QUAGAM.append(-_sp.log(1-F)*Beta)
                continue

        if len(QUAGAM) == 1:
            return(QUAGAM[0])
        return(_sp.array(QUAGAM))

#############################################################

def quagev(F,para):
    U = para[0]
    A = para[1]
    G = para[2]
    if A <= 0:
        print("Parameters Invalid")
        return
    if _is_numeric(F):
        if F <= 0 or F >= 1:
            if F == 0 and G < 0:
                QUAGEV = U+A/G
                return(QUAGEV)
            elif F == 1 and G > 0:
                QUAGEV = U+A/G
                return(QUAGEV)
            else:
                print("F Value Invalid")
                return

    else:
        if any([i >= 1 or i <= 0 for i in F]):
            print("F Value Invalid")
            return
        F = _sp.array(F)
        
    
    Y = -_sp.log(-_sp.log(F))
    if G != 0:
        Y = (1-_sp.exp(-G*Y))/G
    QUAGEV = U+A*Y
    return(QUAGEV)

#############################################################

def quaglo(F,para):
    U = para[0]
    A = para[1]
    G = para[2]
    if A <= 0:
        print("Invalid Parameters")
        return
    if _is_numeric(F):
        if F <= 0 or F >= 1:
            if F == 0 and G < 0:
                QUAGLO = U+A/G
                return(QUAGLO)
            elif F == 1 and G > 0:
                QUAGLO = U+A/G
                return(QUAGLO)
            else:
                print("F Value Invalid")
                return
    else:
        if any([i >= 1 or i <= 0 for i in F]):
            print("F Value Invalid")
            return
        F = _sp.array(F)

    Y = _sp.log(F/(1-F))
    if G != 0:
        Y = (1-_sp.exp(-G*Y))/G
    QUAGLO = U+A*Y
    return(QUAGLO)

#############################################################

def quagno(F,para):
    [U,A,G] = para
    if A <= 0:
        print("Invalid Parameters")
        return
    if _is_numeric(F):
        if F <= 0 or F >= 1:
            if F == 0 and G < 0:
                QUAGNO = U+A/G
                return(QUAGNO)
            elif F == 1 and G > 0:
                QUAGNO = U+A/G
                return(QUAGNO)
            else:
                print("F Value Invalid")
                return
        Y = quastn(F)
    else:
        if any([i >= 1 or i <= 0 for i in F]):
            print("F Value Invalid")
            return
        Y = _sp.array(quastn(F))
        
    
    if G != 0:
        Y = (1-_sp.exp(-G*Y))/G
    QUAGNO = U+A*Y
    return(QUAGNO)    
 
#############################################################

def quagpa(F,para):
    U = para[0]
    A = para[1]
    G = para[2]
    if A <= 0:
        print("Invalid parameters")
        return
    if _is_numeric(F):
        if F <= 0 or F >= 1:
            if F == 0:
                QUAGPA = U
                return(QUAGPA)
            elif F == 1 and G > 0:
                QUAGPA = U + A/G
                return(QUAGPA)
            else:
                print("F Value Invalid")
                return

    else:
        if any([i >= 1 or i <= 0 for i in F]):
            print("F Value Invalid")
            return
        F = _sp.array(F)
        
    Y = -_sp.log(1-F)
    if G !=0:
        Y = (1-_sp.exp(-G*Y))/G
    QUAGPA = U+A*Y
    return(QUAGPA)

#############################################################

def quagum(F,para):

    U = para[0]
    A = para[1]
    
    if A <= 0:
        print("Parameters Invalid")
        return
    if _is_numeric(F):
        if F <= 0 or F >= 1:
            print("F Value Invalid")
            return
    else:
        if any([i >= 1 or i <= 0 for i in F]):
            print("F Value Invalid")
            return
        F = _sp.array(F)

    QUAGUM = U-A*_sp.log(-_sp.log(F))
    return(QUAGUM)

#############################################################

def quakap(F,para):
    U = para[0]
    A = para[1]
    G = para[2]
    H = para[3]
    if A <= 0:
        print("Invalid Parameters")
        return

    if _is_numeric(F):
        if F <= 0 or F>= 1:
            if F==0:
                if H<=0 and G < 0:
                    QUAKAP = U+A/G
                if H<= 0 and G>= 0:
                    print("F Value Invalid")
                    return
                if H > 0 and G!= 0:
                    QUAKAP = U+A/G*(1-H**(-G))
                if H > 0 and G == 0:
                    QUAKAP = U+A*_sp.log(H)

                return(QUAKAP)
            
            if F == 1:
                if G <= 0:
                    print("F Value Invalid")
                    return
                else:
                    QUAKAP = U+A/G
                    return(QUAKAP)

    else:
        if any([i >= 1 or i <= 0 for i in F]):
            print("F Value Invalid")
            return
        F = _sp.array(F)
        
 
    Y = -_sp.log(F)
    if H!=0:
        Y = (1-_sp.exp(-H*Y))/H
        
    Y = -_sp.log(Y)
    if G!= 0:
        Y = (1-_sp.exp(-G*Y))/G
    QUAKAP = U+A*Y
    return(QUAKAP)

#############################################################

def quanor(F,para):
    if para[1] <= 0:
        print("Parameters Invalid")
        return
    if _is_numeric(F):
        if F <= 0 or F >= 1:
            print("F Value Invalid")
            return

    else:
        if any([i >= 1 or i <= 0 for i in F]):
            print("F Value Invalid")
            return
        F = _sp.array(F)
        
    QUANOR = para[0]+para[1]*quastn(F)
    return(QUANOR)

#############################################################

def quape3(F,para):
    SMALL = 1e-6

    if para[1]<= 0:
        print("Paremters Invalid")
        return
    Gamma = para[2]
    if _is_numeric(F):
        if F <= 0 or F >= 1:
            if F == 0 and Gamma >0:
                QUAPE3 = para[0]-2*para[1]/Gamma
                return(QUAPE3)
            elif F == 1 and Gamma < 0:
                QUAPE3 = para[0]-2*para[1]/Gamma
                return(QUAPE3)
            else:
                print("F Value Invalid")
                return
    else:
        if any([i >= 1 or i <= 0 for i in F]):
            print("F Value Invalid")
            return
        F = _sp.array(F)


    if abs(Gamma) < SMALL:
        QUAPE3 = para[0] + para[1]*quastn(F)
        return(QUAPE3)

    Alpha = 4/(Gamma*Gamma)
    Beta = abs(0.5*para[1]*Gamma)
    par = [Alpha,Beta]
    if Gamma > 0:
        gammacalc = quagam(F,par)
        if _is_numeric(gammacalc) == False:
            try:
                if len(gammacalc) == 0:
                    return()
            except:
                return()
        QUAPE3 = para[0]-Alpha*Beta+gammacalc
    if Gamma < 0:
        gammacalc = quagam(1-F,par)
        if len(gammacalc) == 0:
            return()
        QUAPE3 = para[0]+Alpha*Beta-gammacalc
    return(QUAPE3)

#############################################################

def quastn(Flist):
    _split1 = 0.425
    _split2 = 5
    const1 = 0.180625
    const2 = 1.6
    [A0,A1,A2,A3,A4,A5,A6,A7,B1,B2,B3,B4,B5,B6,B7] = [0.338713287279636661e1,
     0.133141667891784377e3,  0.197159095030655144e4,
     0.137316937655094611e5,  0.459219539315498715e5,
     0.672657709270087009e5,  0.334305755835881281e5,
     0.250908092873012267e4,  0.423133307016009113e2,
     0.687187007492057908e3,  0.539419602142475111e4,
     0.212137943015865959e5,  0.393078958000927106e5,
     0.287290857357219427e5,  0.522649527885285456e4]

    [C0,C1,C2,C3,C4,C5,C6,C7,D1,D2,D3,D4,D5,D6,D7] = [0.142343711074968358e1,
     0.463033784615654530e1,  0.576949722146069141e1,
     0.364784832476320461e1,  0.127045825245236838e1,
     0.241780725177450612e0,  0.227238449892691846e-1,
     0.774545014278341408e-3,  0.205319162663775882e1,
     0.167638483018380385e1,  0.689767334985100005e0,
     0.148103976427480075e0,  0.151986665636164572e-1,
     0.547593808499534495e-3,  0.105075007164441684e-8]
                                                      
    [E0,E1,E2,E3,E4,E5,E6,E7,F1,F2,F3,F4,F5,F6,F7] = [0.665790464350110378e1,
     0.546378491116411437e1,  0.178482653991729133e1,
     0.296560571828504891e0,  0.265321895265761230e-1,
     0.124266094738807844e-2,  0.271155556874348758e-4,
     0.201033439929228813e-6,  0.599832206555887938e0,
     0.136929880922735805e0,  0.148753612908506149e-1,
     0.786869131145613259e-3,  0.184631831751005468e-4,
     0.142151175831644589e-6,  0.204426310338993979e-14]
    if _is_numeric(Flist):
        Flist = [Flist]
    QUASTNsave = []
    for F in Flist:
        Q = F-0.5
        if abs(Q) > _split1:
            R=F
            if Q >= 0:
                R = 1-F
            if R <= 0:
                print("F Value Invalid")
            R = _sp.sqrt(-_sp.log(R))
            if R > _split2:
                R = R - _split2
                QUASTN=((((((((E7*R+E6)*R+E5)*R+E4)*R+E3)*R+E2)*R+E1)*R+E0)/
                (((((((F7*R+F6)*R+F5)*R+F4)*R+F3)*R+F2)*R+F1)*R+1))
                if Q < 0:
                    QUASTN = -QUASTN
                
            else:
                R=R-const2
                QUASTN=((((((((C7*R+C6)*R+C5)*R+C4)*R+C3)*R+C2)*R+C1)*R+C0)/
                (((((((D7*R+D6)*R+D5)*R+D4)*R+D3)*R+D2)*R+D1)*R+1))
                if Q < 0:
                    QUASTN = -QUASTN


        else:
            R = const1-Q*Q
            QUASTN = Q*((((((((A7*R+A6)*R+A5)*R+A4)*R+A3)*R+A2)*R+A1)*R+A0)/
            (((((((B7*R+B6)*R+B5)*R+B4)*R+B3)*R+B2)*R+B1)*R+1))
        QUASTNsave.append(QUASTN)

    if len(QUASTNsave) == 1:
        return(QUASTNsave[0])
    return(_sp.array(QUASTNsave))
#############################################################

def quawak(F,para):
    ufl = -170
  
    XI = para[0]
    A = para[1]
    B = para[2]
    C = para[3]
    D = para[4]
    
    fail = 0
    if (B+D) <= 0 and (B != 0 or C != 0 or D!= 0):
        fail = 1
    if A == 0 and B != 0:
        fail = 1
    if C == 0 and D != 0:
        fail = 1
    if C <0 or (A+C)< 0:
        fail = 1
    if A == 0 and C == 0:
        fail = 1

    if fail == 1:
        print("Parameters Invalid")
        return

    if _is_numeric(F):
        if F<= 0 or F>= 1:
            if F == 0:
                QUAWAK = XI
            elif F == 1:
                if D > 0:
                    fail = 1
                if D < 0:
                    QUAWAK = XI + A/B - C/D
                if D == 0 and C > 0:
                    fail = 1
                if D == 0 and C == 0 and B == 0:
                    fail = 1
                if D == 0 and C == 0 and B >0:
                    QUAWAK = XI+A/B

                if fail == 1:
                    print("Function Failed")
                else:
                    return(QUAWAK)

    else:
        if any([i >= 1 or i <= 0 for i in F]):
            print("F Value Invalid")
            return
        F = _sp.array(F)

    if _is_numeric(F):
        Z= [-_sp.log(1-F)]
    else:
        Z = -_sp.log(1-F)
    Y1 = Z
    if B == 0:
        Y2 = Z
        if D !=0:
            Y2 = (1-_sp.exp(D*Y2))/(-D)
        QUAWAK = XI+A*Y1+C*Y2
        return(QUAWAK)
    else:
        QUAWAK = []
        for Zit in Z:
            TEMP = -B*Zit

            if TEMP < ufl:
                Y1 = 1/B
            if TEMP >= ufl:
                Y1 = (1-_sp.exp(TEMP))/B

            Y2 = Zit
            if D !=0:
                Y2 = (1-_sp.exp(D*Y2))/(-D)
            QUAWAK.append(XI+A*Y1+C*Y2)

        if len(QUAWAK) == 1:
            return(QUAWAK[0])
        return(QUAWAK)

#############################################################

def quawei(F,para):
    if len(para) != 3:
        print("Invalid number of parameters for QUAWEI")
        return()
    if para[1] <= 0 or para[2] <= 0:
        print("Invalid Parameters for QUAWEI")
        return()
    if _is_numeric(F):
        if F <= 0 or F >= 1:
            print("F Value Invalid for QUAWEI")

    else:
        if any([i >= 1 or i <= 0 for i in F]):
            print("F Value Invalid for QUAWEI")
            return
        F = _sp.array(F)
        
    return(para[0]+para[1]*((-_sp.log(1-F))**(1/para[2])))
