#----------------------------------------------------------------------
#----------------------------------------------------------------------
# fjk2Dmanifolds : Compute manifolds of fixed points of
#                            the stroboscopic map of the 2D 
#                            reduction of the JFK model
#
# Jen Creaser Apr 2021
#----------------------------------------------------------------------
#----------------------------------------------------------------------

import itertools
import math
import time
#--------------------------------------------------------------------------
# Set the variable names and the parameter names:
load(
     unames = { 1: 'A1',  2: 'C1', 3: 'A2',  4: 'C2'} ,
     parnames = {1: 'B_min', 2: 'I', 3: 'mu', 4: 'taux', 5:'k', 7:'delta',
     9:'gapA',10:'gapC',11:'T', 13:'A',14:'C',15:'MAs',16:'MCs',17:'MAe',18:'MCe',
     20:'L1',21:'L2',22:'v1A',23:'v1C',24:'v2A',25:'v2C',26:'N',27:'L2Norm'}
    )

#=========================================================
# Notes:
    # Man_cont_loads has more than 17 repeats of the orbits in it.
#=========================================================

# set colours for text
headcol = '\x1b[1;33;40m'  # bold yellow
linecol = '\x1b[0;36;40m'  # blue
sublinecol = '\x1b[0;34;40m'    # dark blue
endcol = '\x1b[0;32;40m'    # green

#=========================================================
# Functions
#=========================================================
def manifolds_I(Ivals):

    print(headcol + '\nComputing manifolds for different I (lux) values' + '\x1b[0m')
    print(linecol + '\nContinue PO in I' + '\x1b[0m')
    rn(c='fjk2Dm.1',UZSTOP={'I':max(Ivals)+1},sv='Istart')
    rn(c='fjk2Dm.1',s='Istart',IRS='UZ1',DS='-',UZR={'I':Ivals},UZSTOP={'L2Norm':0.15},sv='Istart')
    
    rl('Istart') 
    compute_orbits('Istart','I',2)    

#=========================================================
def manifolds_T(Ival,Tvals):    
    
    print(headcol + '\nComputing manifolds for different tau_c values' + '\x1b[0m')
    
    parName = 'I' + str("%.2f" % Ival).replace(".", "pt") + '_taux'
    print parName
    
    print(linecol + '\nContinue PO in I' + '\x1b[0m')
    startsol = 'Tstart_I' + str("%.2f" % Ival).replace(".", "pt") 
    dl(startsol)
    print startsol

    rn(c='fjk2Dm.1',sv='start',UZSTOP={'I':[Ival+1]}) # Inc then Dec for multiple I values
    rn(c='fjk2Dm.1',s='start',IRS='UZ1',NPR=0,DS='-',sv='startI',UZR={'I':[Ival]},UZSTOP = {'I':Ival-1})

    data=sl('startI')
    UZdat=data('UZ')
    print UZdat
    lenUZ=len(UZdat)
   
    print(linecol + '\nContinue PO in taux - intrinsic frequency' + '\x1b[0m') 
    for x in range(1,lenUZ):
        stuz = 'UZ' + str("%d" %x)
        print stuz 
        rn(c='fjk2Dm.T',s='startI',IRS=stuz,NMX=500,UZSTOP={27:0.5,4:24.0},sv='startT')
        rn(c='fjk2Dm.T',s='startT',IRS='UZ1',ILP=1,NMX=500,DS='-',UZR={'taux':Tvals},UZSTOP={27:0.5,4:24.0},ap=startsol)

    rl(startsol) 
    compute_orbits(startsol,parName,4)

#=========================================================
def manifolds_N(Ival,Nvals):    
    
    print(headcol + '\nComputing manifolds for different N values' + '\x1b[0m')
    
    parName = 'I' + str("%.2f" % Ival).replace(".", "-") + '_N'
    print parName
    
    startsol = 'Nstart_I' + str("%.2f" % Ival).replace(".", "pt") 
    dl(startsol)
    print startsol

    print(linecol + '\nContinue PO in I' + '\x1b[0m')
    
    rn(c='fjk2Dm.1',sv='start',UZSTOP={'I':[Ival+10]}) # Inc then dec in I for multiple I values
    rn(c='fjk2Dm.1',s='start',IRS='UZ1',NPR=0,DS='-',sv='startI',UZR={'I':[Ival]},UZSTOP = {'L2Norm':0.15})

    data=sl('startI')
    UZdat=data('UZ')
    print UZdat
    lenUZ=len(UZdat)
   
    print(linecol + '\nContinue PO in N - day length' + '\x1b[0m') 
    for x in range(1,lenUZ):
        stuz = 'UZ' + str("%d" %x)
        print stuz 
        rn(c='fjk2Dm.N',s='startI',IRS=stuz,UZSTOP={'N':Nvals},ap=startsol)
        ## For multiple Nvals values use this code:
        #rn(c='fjk2Dm.N',s='startI',IRS=stuz,UZSTOP={'N':max(Nvals)+0.01},NPR=100,sv='startN')
        #rn(c='fjk2Dm.N',s='startN',IRS='UZ1',ILP=1,DS='-',UZR={'N':Nvals},NPR=100,UZSTOP={'N':min(Nvals)-0.01},ap=startsol)
    
    rl(startsol) 
    compute_orbits(startsol,parName,26)

#=========================================================   
def compute_orbits(startsol,parName,pval):
    data=sl(startsol)
    UZdat=data('UZ')
    print UZdat
    lenUZ=len(UZdat)
    s=data()
    
    print(linecol + '\nCompute manifolds for ' + parName + '\x1b[0m')

    for i in range(1,lenUZ+1):
        pt='UZ%d'%(i)
        sol=s(pt)
        spar = sol.PAR(pval)
        spartxt = str("%.2f" % spar).replace(".", "-") 
        orbsave = 'Manifold_orbs17_%s_UZ%d_%s'%(parName,i,spartxt)
        print orbsave
        
        print(sublinecol + '\n' + pt + parName + spartxt + '\x1b[0m')    
        #V1
        rn(c='fjk2Dm.orb', s = startsol,IRS=pt, PAR={30: 1.0}, NPR=100,NMX = 100000,UZSTOP={'delta':[-2e-1,2e-1],'MAs':[-2.5,2.5]}, sv= orbsave)
        rn(c='fjk2Dm.orb', s = startsol,IRS=pt,PAR={30: 1.0}, NPR=100,NMX = 100000,DS = '-',UZSTOP={'delta':[-2e-1,2e-1],'MAs':[-2.5,2.5]}, ap = orbsave)
        print(sublinecol + '\nVector 1 done' + '\x1b[0m')
        time.sleep(1)
 
        #V2
        rn(c='fjk2Dm.orb', s = startsol,IRS=pt,PAR={30: 2.0}, NPR=100,NMX = 100000,UZSTOP={'delta':[-2e-1,2e-1],'MAs':[-2.5,2.5]}, ap = orbsave)
        rn(c='fjk2Dm.orb', s = startsol,IRS=pt,PAR={30: 2.0}, NPR=100,NMX = 100000,DS = '-',UZSTOP={'delta':[-2e-1,2e-1],'MAs':[-2.5,2.5]}, ap = orbsave)
        print(sublinecol + '\nVector 2 done' + '\x1b[0m')
        rl(orbsave)
        time.sleep(1) # to give auto time to catch up so it saves ok
        print(endcol + '\nSaved as ' + orbsave + '\x1b[0m')

        wait()


    cl()

#----------------------------------------------------------------
#----------------------------------------------------------------
# MAIN PROGRAM
#----------------------------------------------------------------
#----------------------------------------------------------------

ld('fjk2Dmanifolds')
#-------------------------------------
# Compute manifolds for I
#-------------------------------------
Ivals =[50.0,1000.0]
manifolds_I(Ivals)

#-------------------------------------
# Compute manifolds for taux
#-------------------------------------
#Ival = 50.0
#Tvals = [23.6,24.6]
#manifolds_T(Ival,Tvals) # Note the unstable solution is complex for these

#Ival = 150.0
#Tvals = [24.6]
#manifolds_T(Ival,Tvals)

#Ival = 1000
#Tvals = [23.2,24.2,25.2]
#manifolds_T(Ival,Tvals) 

#----------------------------------
# Compute manifolds for N
#-------------------------------------
#Ival = 50.0
#Nvals = 17.0
#Nvals = [4.0,10.0,12.0,16.0,20.0]
#manifolds_N(Ival,Nvals)

cl()


