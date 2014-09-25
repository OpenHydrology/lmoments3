import lmoments3 as lmoments


def comparefunc(inp,outp,name,roundnum):
    try:
        check = 0
        for i in range(0,len(inp)):
            if round(inp[i],roundnum) != round(outp[i],roundnum):
                print("ERROR found in "+name)
                print("EXPECTED VALUE: "+str(outp[i]))
                print("ACTUAL VALUE: "+str(inp[i]))
                check = 1
        if check == 0:
            print(name+" Function SUCCESS")
            pass
    except:
        print("ERROR found in "+name)


testdata = [2.0,3.0,4.0,2.4,5.5,1.2,5.4,2.2,7.1,1.3,1.5]
LMU = lmoments.samlmu(testdata)
correctLMU = [3.23636364, 1.14181818, 0.27388535, 0.02335456, -0.04246285]
correctLMU = [3.23636364, 1.14181818, 0.27388535, 0.02335456, -0.04246285]
check = 0
for i in range(0,len(correctLMU)):
    if round(LMU[i],6) != round(correctLMU[i],6):
        print("ERROR found in SAMLMU")
        print("EXPECTED VALUE: "+str(correctLMU[i]))
        print("ACTUAL VALUE: "+str(LMU[i]))
        check = 1
if check == 0:
    print("SAMLMU Function SUCCESS")

print("#######################################")

print("#######################################")
#######################################
##GEV
#######################################
##PELGEV
gevfit = lmoments.pelgev(LMU)
correctgevfit = [2.1792884, 1.3956404, -0.1555609]
comparefunc(gevfit,correctgevfit,"PELGEV",6)

##QUAGEV
gevqua = [lmoments.quagev(0.2,correctgevfit),lmoments.quagev(0.5,correctgevfit),lmoments.quagev(0.8,correctgevfit)]
gevqua2 = lmoments.quagev([0.2,0.5,0.8],correctgevfit)
correctgevqua = [1.539112, 2.705672, 4.537048]
comparefunc(gevqua,correctgevqua,"QUAGEV",6)
comparefunc(gevqua2,correctgevqua,"QUAGEV group",6)

##LMRGEV
gevlmr = lmoments.lmrgev(correctgevfit,4)
correctgevlmr = [3.2363636, 1.1418182, 0.2738854, 0.1998461]
comparefunc(gevlmr,correctgevlmr,"LMRGEV",6)

##CDFGEV
gevcdf = [lmoments.cdfgev(2,correctgevfit),lmoments.cdfgev(5,correctgevfit),lmoments.cdfgev(8,correctgevfit)]
gevcdf2 = lmoments.cdfgev([2,5,8],correctgevfit)
correctgevcdf = [0.3202800, 0.8415637, 0.9606184]
comparefunc(gevcdf,correctgevcdf,"CDFGEV Ind",6)
comparefunc(gevcdf2,correctgevcdf,"CDFGEV Group",6)

##PDFGEV
gevpdf = [lmoments.pdfgev(4,correctgevfit),lmoments.pdfgev(5,correctgevfit),lmoments.pdfgev(6,correctgevfit),lmoments.pdfgev(7,correctgevfit)]
gevpdf2 = lmoments.pdfgev([4,5,6,7],correctgevfit)
correctgevpdf = [0.13388363, 0.07913309, 0.04637529, 0.02757391]
comparefunc(gevpdf,correctgevpdf,"PDFGEV Ind",6)
comparefunc(gevpdf2,correctgevpdf,"PDFGEV group",6)

##LMOMGEV
gevlmom = lmoments.lmomgev(correctgevfit)
correctgevlmom = [3.2363636, 1.1418182, 0.3127273, 0.2281879, 0.1209980]
comparefunc(gevlmom,correctgevlmom,"LMOMGEV",6)

#RANDGEV
try:
    gevrand = lmoments.randgev(10,correctgevfit)
    if len(gevrand) != 10:
        print("RANDGEV FAILED")
except:
    print("RANDGEV FAILED")

print("#######################################")
#######################################
##GLO
#######################################
##PELGLO
glofit = lmoments.pelglo(LMU)
correctglofit = [2.7406580,1.0060517,-0.2738854]
comparefunc(glofit,correctglofit,"PELGLO",6)

##QUAGLO
gloqua = [lmoments.quaglo(0.2,correctglofit),lmoments.quaglo(0.5,correctglofit),lmoments.quaglo(0.8,correctglofit)]
gloqua2 = lmoments.quaglo([0.2,0.5,0.8],correctglofit)
correctgloqua = [1.580189, 2.740658, 4.437061]
comparefunc(gloqua,correctgloqua,"QUAGLO",6)
comparefunc(gloqua2,correctgloqua,"QUAGLO group",6)

##LMRGLO
glolmr = lmoments.lmrglo(correctglofit,4)
correctglolmr = [3.2363636, 1.1418182, 0.2738854, 0.2291777]
comparefunc(glolmr,correctglolmr,"LMRGLO",6)

##CDFGLO
glocdf = [lmoments.cdfglo(2,correctglofit),lmoments.cdfglo(5,correctglofit),lmoments.cdfglo(8,correctglofit)]
glocdf2 = lmoments.cdfglo([2,5,8],correctglofit)
correctglocdf = [0.3052960, 0.8519915, 0.9624759]
comparefunc(glocdf,correctglocdf,"CDFGLO ind",6)
comparefunc(glocdf2,correctglocdf,"CDFGLO group",6)

##PDFGLO
glopdf = [lmoments.pdfglo(4,correctglofit),lmoments.pdfglo(5,correctglofit),lmoments.pdfglo(6,correctglofit),lmoments.pdfglo(7,correctglofit)]
glopdf2 = lmoments.pdfglo([4,5,6,7],correctglofit)
correctglopdf = [0.14033225, 0.07760825, 0.04294245, 0.02463028]
comparefunc(glopdf,correctglopdf,"PDFGLO Ind",6)
comparefunc(glopdf2,correctglopdf,"PDFGLO group",6)

##LMOMGLO
glolmom = lmoments.lmomglo(correctglofit)
correctglolmom = [3.2363636, 1.1418182, 0.3127273, 0.2616792]
comparefunc(glolmom,correctglolmom,"LMOMGLO",6)

#RANDGLO
try:
    glorand = lmoments.randglo(10,correctglofit)
    if len(glorand) != 10:
        print("RANDGLO FAILED")
except:
    print("RANDGLO FAILED")
    
print("#######################################")
#######################################
##GNO
#######################################
##PELGNO
gnofit = lmoments.pelgno(LMU)
correctgnofit = [2.6888917, 1.7664322, -0.5707506]
comparefunc(gnofit,correctgnofit,"PELGNO",6)

##QUAGNO
gnoqua = [lmoments.quagno(0.2,correctgnofit),lmoments.quagno(0.5,correctgnofit),lmoments.quagno(0.8,correctgnofit)]
gnoqua2 = lmoments.quagno([0.2,0.5,0.8],correctgnofit)
correctgnoqua = [1.508372, 2.688892, 4.597378]
comparefunc(gnoqua,correctgnoqua,"QUAGNO",6)
comparefunc(gnoqua2,correctgnoqua,"QUAGNO group",6)

##LMRGNO
gnolmr = lmoments.lmrgno(correctgnofit,4)
correctglolmr = [3.2363636, 1.1418182, 0.2738848, 0.1818274]
comparefunc(gnolmr,correctglolmr,"LMRGNO",6)

##CDFGNO
gnocdf = [lmoments.cdfgno(2,correctgnofit),lmoments.cdfgno(5,correctgnofit),lmoments.cdfgno(8,correctgnofit)]
gnocdf2 = lmoments.cdfgno([2,5,8],correctgnofit)
correctgnocdf = [0.3295539, 0.8357710, 0.9599970]
comparefunc(gnocdf,correctgnocdf,"CDFGNO Ind",6)
comparefunc(gnocdf2,correctgnocdf,"CDFGNO Group",6)

##PDFGNO
gnopdf = [lmoments.pdfgno(4,correctgnofit),lmoments.pdfgno(5,correctgnofit),lmoments.pdfgno(6,correctgnofit),lmoments.pdfgno(7,correctgnofit)]
gnopdf2 = lmoments.pdfgno([4,5,6,7],correctgnofit)
correctgnopdf = [0.13099439, 0.08020770, 0.04842820, 0.02933547]
comparefunc(gnopdf,correctgnopdf,"PDFGNO Ind",6)
comparefunc(gnopdf2,correctgnopdf,"PDFGNO group",6)

##LMOMGNO
gnolmom = lmoments.lmomgno(correctgnofit)
correctgnolmom = [3.2363636, 1.1418182, 0.3127266, 0.2076140, 0.1104614]
comparefunc(gnolmom,correctgnolmom,"LMOMGNO",6)

#RANDGNO
try:
    gnorand = lmoments.randgno(10,correctgnofit)
    if len(gnorand) != 10:
        print("RANDGNO FAILED")
except:
    print("RANDGNO FAILED")
    
print("#######################################")
#######################################
##GPA
#######################################
##PELGPA
gpafit = lmoments.pelgpa(LMU)
correctgpafit = [0.7928727,2.7855796,0.1400000]
comparefunc(gpafit,correctgpafit,"PELGPA",6)

##QUAGPA
gpaqua = [lmoments.quagpa(0.2,correctgpafit),lmoments.quagpa(0.5,correctgpafit),lmoments.quagpa(0.8,correctgpafit)]
gpaqua2 = lmoments.quagpa([0.2,0.5,0.8],correctgpafit)
correctgpaqua = [1.404848, 2.632964, 4.806899]
comparefunc(gpaqua,correctgpaqua,"QUAGPA",6)
comparefunc(gpaqua2,correctgpaqua,"QUAGPA group",6)

##LMRGPA
gpalmr = lmoments.lmrgpa(correctgpafit,4)
correctgpalmr = [3.2363636, 1.1418182, 0.2738854, 0.1230499 ]
comparefunc(gpalmr,correctgpalmr,"LMRGPA",6)

##CDFGPA
gpacdf = [lmoments.cdfgpa(2,correctgpafit),lmoments.cdfgpa(5,correctgpafit),lmoments.cdfgpa(8,correctgpafit)]
gpacdf2 = lmoments.cdfgpa([2,5,8],correctgpafit)
correctgpacdf = [0.3604888, 0.8167330, 0.9597484]
comparefunc(gpacdf,correctgpacdf,"CDFGNO Ind",6)
comparefunc(gpacdf2,correctgpacdf,"CDFGNO Group",6)

##PDFGPA
gpapdf = [lmoments.pdfgpa(4,correctgpafit),lmoments.pdfgpa(5,correctgpafit),lmoments.pdfgpa(6,correctgpafit),lmoments.pdfgpa(7,correctgpafit)]
gpapdf2 = lmoments.pdfgpa([4,5,6,7],correctgpafit)
correctgpapdf = [0.12194724, 0.08343282, 0.05567275, 0.03610412]
comparefunc(gpapdf,correctgpapdf,"PDFGPA Ind",6)
comparefunc(gpapdf2,correctgpapdf,"PDFGPA group",6)

##LMOMGPA
gpalmom = lmoments.lmomgpa(correctgpafit)
correctgpalmom = [3.23636364, 1.14181818, 0.31272727, 0.14050066, 0.07817741]
comparefunc(gpalmom,correctgpalmom,"LMOMGPA",6)

#RANDGPA
try:
    gparand = lmoments.randgpa(10,correctgpafit)
    if len(gparand) != 10:
        print("RANDGPA FAILED")
except:
    print("RANDGPA FAILED")
    
print("#######################################")
#######################################
##GUM
#######################################
##PELGUM
gumfit = lmoments.pelgum(LMU)
correctgumfit = [2.285519, 1.647295]
comparefunc(gumfit,correctgumfit,"PELGUM",6)

##QUAGUM
gumqua = [lmoments.quagum(0.2,correctgumfit),lmoments.quagum(0.5,correctgumfit),lmoments.quagum(0.8,correctgumfit)]
gumqua2 = lmoments.quagum([0.2,0.5,0.8],correctgumfit)
correctgumqua = [1.501596, 2.889274, 4.756363]
comparefunc(gumqua,correctgumqua,"QUAGUM",6)
comparefunc(gumqua2,correctgumqua,"QUAGUM group",6)

##LMRGUM
gumlmr = lmoments.lmrgum(correctgumfit,4)
correctgumlmr = [3.236363, 1.141818, 0.169925, 0.150375]
comparefunc(gumlmr,correctgumlmr,"LMRGUM",6)

##CDFGUM
gumcdf = [lmoments.cdfgum(2,correctgumfit),lmoments.cdfgum(5,correctgumfit),lmoments.cdfgum(8,correctgumfit)]
gumcdf2 = lmoments.cdfgum([2,5,8],correctgumfit)
correctgumcdf = [0.3044484, 0.8249232, 0.9693322]
comparefunc(gumcdf,correctgumcdf,"CDFGUM Ind",6)
comparefunc(gumcdf2,correctgumcdf,"CDFGUM Group",6)

##PDFGUM
gumpdf = [lmoments.pdfgum(4,correctgumfit),lmoments.pdfgum(5,correctgumfit),lmoments.pdfgum(6,correctgumfit),lmoments.pdfgum(7,correctgumfit)]
gumpdf2 = lmoments.pdfgum([4,5,6,7],correctgumfit)
correctgumpdf = [0.15060460, 0.09638151, 0.05733088, 0.03276992]
comparefunc(gumpdf,correctgumpdf,"PDFGUM Ind",6)
comparefunc(gumpdf2,correctgumpdf,"PDFGUM group",6)

##LMOMGUM
gumlmom = lmoments.lmomgum(correctgumfit)
correctgumlmom = [3.2363635, 1.1418182, 0.1940235, 0.1717009, 0.0637914]
comparefunc(gumlmom,correctgumlmom,"LMOMGUM",6)

#RANDGUM
try:
    gumrand = lmoments.randgum(10,correctgumfit)
    if len(gumrand) != 10:
        print("RANDGUM FAILED")
except:
    print("RANDGUM FAILED")
    
print("#######################################")
#######################################
##KAP
#######################################
##PELKAP
kapfit = lmoments.pelkap(LMU)
correctkapfit = [-9.0633543, 17.0127900, 0.9719618, 2.4727933]
comparefunc(kapfit,correctkapfit,"PELKAP",6)
    
##QUAKAP
kapqua = [lmoments.quakap(0.2,correctkapfit),lmoments.quakap(0.5,correctkapfit),lmoments.quakap(0.8,correctkapfit)]
kapqua2 = lmoments.quakap([0.2,0.5,0.8],correctkapfit)
correctkapqua = [1.311688, 2.454434, 5.286237]
comparefunc(kapqua,correctkapqua,"QUAKAP",6)
comparefunc(kapqua2,correctkapqua,"QUAKAP group",6)

##LMRKAP
kaplmr = lmoments.lmrkap(correctkapfit,4)
correctkaplmr = [3.23636364, 1.14181818, 0.27388545, 0.02335466]
comparefunc(kaplmr,correctkaplmr,"LMRKAP",6)

##CDFKAP
kapcdf = [lmoments.cdfkap(2,correctkapfit),lmoments.cdfkap(5,correctkapfit),lmoments.cdfkap(8,correctkapfit)]
kapcdf2 = lmoments.cdfkap([2,5,8],correctkapfit)
correctkapcdf = [0.4185230, 0.7772538, 0.9769973]
comparefunc(kapcdf,correctkapcdf,"CDFKAP Ind",6)
comparefunc(kapcdf2,correctkapcdf,"CDFKAP Group",6)

##PDFKAP
kappdf = [lmoments.pdfkap(4,correctkapfit),lmoments.pdfkap(5,correctkapfit),lmoments.pdfkap(6,correctkapfit),lmoments.pdfkap(7,correctkapfit)]
kappdf2 = lmoments.pdfkap([4,5,6,7],correctkapfit)
correctkappdf = [0.09794121, 0.08128701, 0.07022161, 0.06197593]
comparefunc(kappdf,correctkappdf,"PDFKAP Ind",6)
comparefunc(kappdf2,correctkappdf,"PDFKAP group",6)

##LMOMKAP
kaplmom = lmoments.lmomkap(correctkapfit)
correctkaplmom = [3.236364, 1.141818, 0.09662925, 0.008239735, 0.00005919404]
comparefunc(kaplmom,correctkaplmom,"LMOMKAP",6)

#RANDKAP
try:
    kaprand = lmoments.randkap(10,correctkapfit)
    if len(kaprand) != 10:
        print("RANDKAP FAILED")
except:
    print("RANDKAP FAILED")
    
print("#######################################")
#######################################
##NOR
#######################################
##PELNOR
norfit = lmoments.pelnor(LMU)
correctnorfit = [3.236364, 2.023820]
comparefunc(norfit,correctnorfit,"PELNOR",6)

##QUANOR
norqua = [lmoments.quanor(0.2,correctnorfit),lmoments.quanor(0.5,correctnorfit),lmoments.quanor(0.8,correctnorfit)]
norqua2 = lmoments.quanor([0.2,0.5,0.8],correctnorfit)
correctnorqua = [1.533074, 3.236364, 4.939654]
comparefunc(norqua,correctnorqua,"QUANOR",6)
comparefunc(norqua2,correctnorqua,"QUANOR group",6)

##LMRNOR
norlmr = lmoments.lmrnor(correctnorfit,4)
correctnorlmr = [3.2363636, 1.1418182, 0.0000000, 0.1226017]
comparefunc(norlmr,correctnorlmr,"LMRNOR",6)

##CDFNOR
norcdf = [lmoments.cdfnor(2,correctnorfit),lmoments.cdfnor(5,correctnorfit),lmoments.cdfnor(8,correctnorfit)]
norcdf2 = lmoments.cdfnor([2,5,8],correctnorfit)
correctnorcdf = [0.2706309, 0.8082428, 0.9907083]
comparefunc(norcdf,correctnorcdf,"CDFNOR Ind",6)
comparefunc(norcdf2,correctnorcdf,"CDFNOR Group",6)

##PDFNOR
norpdf = [lmoments.pdfnor(4,correctnorfit),lmoments.pdfnor(5,correctnorfit),lmoments.pdfnor(6,correctnorfit),lmoments.pdfnor(7,correctnorfit)]
norpdf2 = lmoments.pdfnor([4,5,6,7],correctnorfit)
correctnorpdf = [0.18357864, 0.13484509, 0.07759170, 0.03497539]
comparefunc(norpdf,correctnorpdf,"PDFNOR Ind",6)
comparefunc(norpdf2,correctnorpdf,"PDFNOR group",6)

##LMOMNOR
norlmom = lmoments.lmomnor(correctnorfit)
correctnorlmom = [3.2363636, 1.1418182, 0.0000000, 0.1399889, 0.0000000]
comparefunc(norlmom,correctnorlmom,"LMOMNOR",6)

#RANDNOR
try:
    norrand = lmoments.randnor(10,correctnorfit)
    if len(norrand) != 10:
        print("RANDNOR FAILED")
except:
    print("RANDNOR FAILED")
    
print("#######################################")
#######################################
##PE3
#######################################
##PELPE3
pe3fit = lmoments.pelpe3(LMU)
correctpe3fit = [3.236364, 2.199489, 1.646184]
comparefunc(pe3fit,correctpe3fit,"PELPE3",6)

##QUAPE3
pe3qua = [lmoments.quape3(0.2,correctpe3fit),lmoments.quape3(0.5,correctpe3fit),lmoments.quape3(0.8,correctpe3fit)]
pe3qua2 = lmoments.quape3([0.2,0.5,0.8],correctpe3fit)
correctpe3qua = [1.447672, 2.663015, 4.705896]
comparefunc(pe3qua,correctpe3qua,"QUAPE3",6)
comparefunc(pe3qua2,correctpe3qua,"QUAPE3 group",6)

##LMRPE3
pe3lmr = lmoments.lmrpe3(correctpe3fit,4)
correctpe3lmr = [3.2363636, 1.1418182, 0.2738845, 0.1498865]
comparefunc(pe3lmr,correctpe3lmr,"LMRPE3",6)

##CDFPE3
pe3cdf = [lmoments.cdfpe3(2,correctpe3fit),lmoments.cdfpe3(5,correctpe3fit),lmoments.cdfpe3(8,correctpe3fit)]
pe3cdf2 = lmoments.cdfpe3([2,5,8],correctpe3fit)
correctpe3cdf = [0.3462110, 0.8258929, 0.9597978]
comparefunc(pe3cdf,correctpe3cdf,"CDFPE3 Ind",6)
comparefunc(pe3cdf2,correctpe3cdf,"CDFPE3 Group",6)

##PDFPE3
pe3pdf = [lmoments.pdfpe3(4,correctpe3fit),lmoments.pdfpe3(5,correctpe3fit),lmoments.pdfpe3(6,correctpe3fit),lmoments.pdfpe3(7,correctpe3fit)]
pe3pdf2 = lmoments.pdfpe3([4,5,6,7],correctpe3fit)
correctpe3pdf = [0.12681904, 0.08243435, 0.05226950, 0.03260397]
comparefunc(pe3pdf,correctpe3pdf,"PDFPE3 Ind",6)
comparefunc(pe3pdf2,correctpe3pdf,"PDFPE3 group",6)

##LMOMPE3
pe3lmom = lmoments.lmompe3(correctpe3fit)
correctpe3lmom = [3.2363636, 1.1418182, 0.3127263, 0.1711432]
comparefunc(pe3lmom,correctpe3lmom,"LMOMPE3",6)

#RANDPE3
try:
    pe3rand = lmoments.randpe3(10,correctpe3fit)
    if len(pe3rand) != 10:
        print("RANDPE3 FAILED")
except:
    print("RANDPE3 FAILED")
    
print("#######################################")
#######################################
##WAK
#######################################
##PELWAK
wakfit = lmoments.pelwak(LMU)
correctwakfit = [0.7928727, 2.7855796, 0.1400000, 0.0000000, 0.0000000]
comparefunc(wakfit,correctwakfit,"PELWAK",6)

##QUAWAK
wakqua = [lmoments.quawak(0.2,correctwakfit),lmoments.quawak(0.5,correctwakfit),lmoments.quawak(0.8,correctwakfit)]
wakqua2 = lmoments.quawak([0.2,0.5,0.8],correctwakfit)
correctwakqua = [1.404848, 2.632964, 4.806899]
comparefunc(wakqua,correctwakqua,"QUAWAK",6)
comparefunc(wakqua2,correctwakqua,"QUAWAK group",6)

##LMRWAK
waklmr = lmoments.lmrwak(correctwakfit,4)
correctwaklmr = [3.2363636, 1.1418182, 0.2738854, 0.1230499]
comparefunc(waklmr,correctwaklmr,"LMRWAK",6)

##CDFWAK
wakcdf = [lmoments.cdfwak(2,correctwakfit),lmoments.cdfwak(5,correctwakfit),lmoments.cdfwak(8,correctwakfit)]
wakcdf2 = lmoments.cdfwak([2,5,8],correctwakfit)
correctwakcdf = [0.3604888, 0.8167330, 0.9597484]
comparefunc(wakcdf,correctwakcdf,"CDFWAK Ind",6)
comparefunc(wakcdf2,correctwakcdf,"CDFWAK Group",6)

##PDFWAK
wakpdf = [lmoments.pdfwak(4,correctwakfit),lmoments.pdfwak(5,correctwakfit),lmoments.pdfwak(6,correctwakfit),lmoments.pdfwak(7,correctwakfit)]
wakpdf2 = lmoments.pdfwak([4,5,6,7],correctwakfit)
correctwakpdf = [0.12194724, 0.08343282, 0.05567275, 0.03610412]
comparefunc(wakpdf,correctwakpdf,"PDFWAK Ind",6)
comparefunc(wakpdf2,correctwakpdf,"PDFWAK group",6)

##LMOMWAK
waklmom = lmoments.lmomwak(correctwakfit)
correctwaklmom = [3.23636364, 1.14181818, 0.31272727, 0.14050066, 0.07817741]
comparefunc(waklmom,correctwaklmom,"LMOMWAK",6)

##RANDWAK
try:
    wakrand = lmoments.randwak(10,correctwakfit)
    if len(wakrand) != 10:
        print("RANDWAK FAILED")
except:
    print("RANDWAK FAILED")
    
print("#######################################")    
#######################################
##WEI
#######################################
##PELWEI
weifit = lmoments.pelwei(LMU)
correctweifit = [0.6740393, 2.7087887, 1.1750218]
comparefunc(weifit,correctweifit,"PELWEI",6)

##QUAWEI
weiqua = [lmoments.quawei(0.2,correctweifit),lmoments.quawei(0.5,correctweifit),lmoments.quawei(0.8,correctweifit)]
weiqua2 = lmoments.quawei([0.2,0.5,0.8],correctweifit)
correctweiqua = [1.429808, 2.656981, 4.735337]
comparefunc(weiqua,correctweiqua,"QUAWEI",6)
comparefunc(weiqua2,correctweiqua,"QUAWEI group",6)

##LMRWEI
weilmr = lmoments.lmrwei(correctweifit,4)
correctweilmr = [3.2363636, 1.1418182, 0.2738853, 0.1413359]
comparefunc(weilmr,correctweilmr,"LMRWEI",6)

##CDFWEI
weicdf = [lmoments.cdfwei(2,correctweifit),lmoments.cdfwei(5,correctweifit),lmoments.cdfwei(8,correctweifit)]
weicdf2 = lmoments.cdfwei([2,5,8],correctweifit)
correctweicdf = [0.3507727, 0.8233116, 0.9600031 ]
comparefunc(weicdf,correctweicdf,"CDFWEI Ind",6)
comparefunc(weicdf2,correctweicdf,"CDFWEI Group",6)

##PDFWEI
weipdf = [lmoments.pdfwei(4,correctweifit),lmoments.pdfwei(5,correctweifit),lmoments.pdfwei(6,correctweifit),lmoments.pdfwei(7,correctweifit)]
weipdf2 = lmoments.pdfwei([4,5,6,7],correctweifit)
correctweipdf = [0.07149587, 0.04550752, 0.02836919, 0.01738135]
comparefunc(weipdf,correctweipdf,"PDFWEI Ind",6)
comparefunc(weipdf2,correctweipdf,"PDFWEI group",6)

##LMOMWEI
weilmom = lmoments.lmomwei(correctweifit)
correctweilmom = [1.88828496, 1.14181818, 0.31272723, 0.16137985, 0.09159867]
comparefunc(weilmom,correctweilmom,"LMOMWEI",6)

##RANDWEI
try:
    weirand = lmoments.randwei(10,correctweifit)
    if len(weirand) != 10:
        print("RANDWAK FAILED")
except:
    print("RANDWAK FAILED")
    
#######################################
##NLogL
#######################################
data = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 2.4, 3.5, 1.4, 6.5, 1.2, 6.8, 5.4, 3.4]
gamfit = lmoments.pelgam(lmoments.samlmu(data,5))
test1 = lmoments.NlogL(data,"GAM",gamfit)
test2 = lmoments.NlogL(data,"GAM")
print("NLogL Test Succeeded")


