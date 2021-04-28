#----------------------------------------------------------------------
#----------------------------------------------------------------------
# fjk2Dbifurcations : Compute bifurcations of the periodic 
#                               orbits in the 2D reduction of the FJK 
#                               model
#
# Jen Creaser Apr 2021
#----------------------------------------------------------------------
# to run code type the following into the command line
#                 auto_d fjk2Dbifurcations.py
#----------------------------------------------------------------------

import itertools
import math
import numpy as np

#--------------------------------------------------------------------------
# Set the variable names and the parameter names:
load(
     unames = { 1: 'A1',  2: 'C1', 3: 'A2',  4: 'C2'} ,
     parnames = {1: 'B_min', 2: 'lux', 3: 'mu', 4: 'taux', 5:'k', 7:'delta',
     9:'gapA',10:'gapC',11:'T', 13:'A',14:'C',15:'MAs',16:'MCs',17:'MAe',18:'MCe',
     20:'L1r',21:'L2r',26:'L1i',27:'L2i',28:'NS',30:'N',31:'L2Norm'}
    )

# set colours for text
headcol = '\x1b[1;33;40m'  # bold yellow
linecol = '\x1b[0;36;40m'  # blue
sublinecol = '\x1b[0;34;40m'    # dark blue
endcol = '\x1b[0;32;40m'    # green


#=========================================================
#=========================================================
# Light Intensity - I
#=========================================================
#=========================================================
# Continue just the po
def bif_luxpo():

    print(headcol + '\nContinue PO in lux' + '\x1b[0m')
    rn(c='fjk2Db.luxpo',sv='bif_lux_po')
    rn(c='fjk2Db.luxpo',DS='-',ap='bif_lux_po')
    
    rl('bif_lux_po')

#=========================================================
# Continue the po with the variational problem
def bif_lux():

    print(headcol + '\nContinue PO in lux with variational problem' + '\x1b[0m')
    rn(c='fjk2Db.lux',sv='bif_lux_start')
    rn(c='fjk2Db.lux',DS='-',s='bif_lux_start',IRS=8,sv='bif_lux_start')
    rl('bif_lux_start')
    

#=========================================================
#=========================================================
# Intrinstic period - taux
#=========================================================
#=========================================================
def bif_taux_l50():
    print(headcol + '\nContinue PO in taux for lux = 50' + '\x1b[0m')
    tnam = 'bif_taux_l50'
    snam = 'bif_lux_start' # computed above
    dl(tnam)   
     
    for x in [5,7]: # choose starting solution numbers for I=50
        rn(c='fjk2Db.taux',s=snam,IRS = x,UZSTOP={20:[-1.0,1.0],21:[1.0-1.0]},UZR={},sv='tstart')
        rn(c='fjk2Db.taux',s='tstart',IRS = 'UZ1',DS='-',UZSTOP={20:[-1.0,1.0],21:[1.0-1.0]},UZR={},ap=tnam)
        rl(tnam)
        wait()
        
    rn(c='fjk2Db.taux',s=snam,IRS = 1,ap=tnam)
    rn(c='fjk2Db.taux',s=snam,IRS = 1 ,DS='-',ap=tnam)
    rl(tnam)
    wait() 
    
#=========================================================
def bif_taux_l150():        
    print(headcol + '\nContinue PO in taux for lux=150' + '\x1b[0m')
    snam = 'bif_lux_start'
    tnam2 = 'bif_lux_taux_ext150'
    dl(tnam2)
    rn(c='fjk2Db.taux',s=snam,IRS = 9,DS='-',ap=tnam2)
    rn(c='fjk2Db.taux',s=snam,IRS = 9,ap=tnam2)    
    rl(tnam2)
    wait()
    
#=========================================================
def bif_taux_l1000():        
    print(headcol + '\nContinue PO in taux for lux=1000' + '\x1b[0m')
    snam = 'bif_lux_start'
    tnam2 = 'bif_taux_l1000'
    dl(tnam2)
    rn(c='fjk2Db.taux',s=snam,IRS = 10,DS='-',ap=tnam2)
    rn(c='fjk2Db.taux',s=snam,IRS = 10,ap=tnam2)    
    rl(tnam2)
    wait()
    
#=========================================================
#=========================================================
# Day Length - N
#=========================================================
#=========================================================
def bif_N_l50():
    
    print('\x1b[1;33;40m' + '\nContinue PO in N for lux=50' + '\x1b[0m')
    snam = 'bif_lux_start'
    nnam = 'bif_N_l50'
    dl(nnam)    

    rn(c='fjk2Db.N',s=snam,IRS=1,ap=nnam)	
    rn(c='fjk2Db.N',s=snam,IRS=1,DS='-',ap=nnam)
    wait()
    
    for x in [5,7]:                                                                                                                
        rn(c='fjk2Db.N',s=snam,IRS = x,UZSTOP={20:[-1.0,1.0],21:[1.0-1.0]},UZR={},sv='nstart') 
        rn(c='fjk2Db.N',s='nstart',IRS = 'UZ1',DS='-',UZSTOP={20:[-1.0,1.0],21:[1.0-1.0]},UZR={},ap=nnam)  
        wait()

    rl(nnam)
    wait()
    
#=========================================================
def bif_N_l150():
    
    print('\x1b[1;33;40m' + '\nContinue PO in N for lux=150' + '\x1b[0m')
    snam = 'bif_lux_start'
    nnam = 'bif_N_l150'
    dl(nnam)    
            
    rn(c='fjk2Db.N',s=snam,IRS = 9,ap=nnam)
    rn(c='fjk2Db.N',s=snam,IRS = 9,DS='-',ap=nnam)
    rl(nnam)
    pl(nnam)
    wait()
    
#=========================================================
def bif_N_l1000():
    
    print('\x1b[1;33;40m' + '\nContinue PO in N for lux=1000' + '\x1b[0m')
    snam = 'bif_lux_start'
    nnam = 'bif_N_l1000'
    dl(nnam)    
            
    rn(c='fjk2Db.N',s=snam,IRS = 10,ap=nnam)
    rn(c='fjk2Db.N',s=snam,IRS = 10,DS='-',ap=nnam)
    rl(nnam)
    pl(nnam)
    wait()
    

#=========================================================
#=========================================================
# Lux and taux
#=========================================================
#=========================================================
def bif_lt():    
    uzpar = 2           # 2 is Lux - this is just to identify a UZ point for the lux-N computation below
    uzvals = 200.0   

    print(headcol + '\nLux-taux bif diag - following folds'  + '\x1b[0m')
    snam = 'bif_taux_l150' # solns 1,5 lambda1real=1;  2,6 lambda2real=1 are folds
    fnam = 'bif_lux_taux_ext'
    dl(fnam)

    for x in [1,5]:
        pt='UZ%d'%(x)
        print(linecol + '\n Following Fold 1 ' + pt + '\x1b[0m')   
        rn(c='fjk2Db.lt',s=snam, IRS = pt, NPR=0,   ICP =  [2,4,13,14,28,26,21,27,31], UZSTOP={4:[3.0,24.01,23.99],2:[1E-3,1000.0,1001.0]},ap=fnam,UZR={uzpar:[uzvals]})
        rn(c='fjk2Db.lt',s=snam,DS='-',IRS=pt,NPR=0,ICP=  [2,4,13,14,28,26,21,27,31], UZSTOP={4:[3.0,24.01,23.99],2:[1E-3,1000.0,1001.0]}, ap=fnam,UZR={uzpar:[uzvals]})

    for x in [2]: # don't need UZ6 as well, 2 does it all.
        pt='UZ%d'%(x)
        print(linecol + '\n Following Fold 2 ' + pt + '\x1b[0m')   
        rn(c='fjk2Db.lt',s=snam, IRS = pt,NPR=0,  ICP =  [2,4,13,14,28,20,26,27,31],UZR={uzpar:[uzvals]}, UZSTOP={4:[3.0],2:[1E-3,1000.0,1001.0]},ap=fnam)
        rn(c='fjk2Db.lt',s=snam,DS='-',IRS=pt,NPR=0,ICP= [2,4,13,14,28,20,26,27,31],UZR={uzpar:[uzvals]}, UZSTOP={4:[3.0],2:[1E-3,1000.0,1001.0]}, ap=fnam)
     
    print(headcol + '\nLux-taux bif diag - following NS'  + '\x1b[0m')
    snam = 'bif_taux_l1000' # point 5 is NS
     
    for x in [1,5]:
        pt='UZ%d'%(x)
        print(linecol + '\n Following NS ' + pt + '\x1b[0m')   
        rn(c='fjk2Db.lt',s=snam, IRS = pt,  ICP =  [2,4,13,14,20,26,21,27,31],DSMIN=1e-3, DSMAX=1e0,UZR={uzpar:[uzvals]},  ap=fnam)
        rn(c='fjk2Db.lt',s=snam, DS='-',IRS = pt,  ICP =  [2,4,13,14,20,26,21,27,31],DSMIN=1e-3, DSMAX=1e0,UZR={uzpar:[uzvals]},ap=fnam)
     

#=========================================================
#=========================================================
# Lux and N for taux fixed
#=========================================================
#=========================================================
def bif_IN():   

     fnam = 'bif_lux_N_taux24pt2'
     dl(fnam) 
     
     print(headcol + '\nFollow folds in lux and N' + '\x1b[0m')
     startsol = 'bif_lux_taux_ext' # UZ 1 for f1 and 5 for f2
    
     for x in [3,4]: # first fold
         pt='UZ%d'%(x)
         print(linecol + '\n Following Fold 1 ' + pt + '\x1b[0m')    
         rn(c='fjk2Db.lt',s=startsol, IRS = pt, DSMIN = 1E-3,DSMAX = 1e0,ICP = [2,30,13,14,28,26,21,27,31],UZSTOP={30:[0.0,24.0],2:[1E-3,2000.0]},ap=fnam)
         rn(c='fjk2Db.lt',s=startsol, DS='-',DSMIN = 1E-3,DSMAX = 1e0,IRS = pt, ICP =  [2,30,13,14,28,26,21,27,31],UZSTOP={30:[0.0,24.0],2:[1E-3,2000.0]},ap=fnam)

     for x in []: # second fold 
        pt='UZ%d'%(x)
        print(linecol + '\n Following Fold 2 ' + pt + '\x1b[0m')    
        rn(c='fjk2Db.lt',s=startsol, IRS = pt, DSMIN = 1E-3,DSMAX = 1e0,ICP =  [2,30,13,14,28,20,26,27,31],UZSTOP={30:[0.0,24.0],2:[1E-3,2000.0]},ap=fnam)
        rn(c='fjk2Db.lt',s=startsol, DS='-',DSMIN = 1E-3,DSMAX = 1e0,IRS = pt, ICP =  [2,30,13,14,28,20,26,27,31],UZSTOP={30:[0.0,24.0],2:[1E-3,2000.0]},ap=fnam)


     for x in []: # NS
        pt='UZ%d'%(x)
        print(linecol + '\n Following NS in N and lux ' + pt + '\x1b[0m')    
        rn(c='fjk2Db.lt',s=startsol, IRS = pt, ICP =  [2,30,13,14,20,26,21,27,31],DSMIN = 1E-3,DSMAX = 1e0,UZSTOP={30:[0.0,24.0],2:[1e-3,2000.0]},UZR={},ap=fnam)
        rn(c='fjk2Db.lt',s=startsol,DS='-',IRS=pt,ICP=[2,30,13,14,20,26,21,27,31],DSMIN = 1E-3,DSMAX = 1e0,UZSTOP={30:[0.0,24.0],2:[1e-3,2000.0]},UZR={},ap=fnam)

     wait()
     # There are no ns points for tau=24.2
     
#=========================================================
#=========================================================
# taux and N for lux fixed
#=========================================================
#=========================================================
def bif_tN():   

     startsol = 'bif_lux_taux_ext' # UZ 1 for f1 and 5 for f2
     fnam = 'bif_taux_N_lux200'
     dl(fnam) 
     
     print(headcol + '\nFollow folds in taux and N' + '\x1b[0m')
    
     for x in [2,4,5]: # first fold
         pt='UZ%d'%(x)
         print(linecol + '\n Following Fold 1 ' + pt + '\x1b[0m')    
         rn(c='fjk2Db.lt',s=startsol, IRS = pt, DSMIN = 1E-4,DSMAX = 1e-2,ICP = [4,30,13,14,28,26,21,27,31],UZSTOP={30:[0.0,24.0],4:[3.0,30.0]},UZR={},ap=fnam)
         rn(c='fjk2Db.lt',s=startsol,IRS=pt,DS='-',DSMIN=1E-4,DSMAX= 1e-2,ICP = [4,30,13,14,28,26,21,27,31],UZSTOP={30:[0.0,24.0],4:[3.0,30.0]},UZR={},ap=fnam)

     for x in [6]: # second fold 
        pt='UZ%d'%(x)
        print(linecol + '\n Following Fold 2 ' + pt + '\x1b[0m')    
        rn(c='fjk2Db.lt',s=startsol, IRS = pt, DSMIN = 1E-4,DSMAX = 1e-2,ICP = [4,30,13,14,28,20,26,27,31],UZSTOP={30:[0.0,24.0],4:[3.0,30.0]},UZR={},ap=fnam)
        rn(c='fjk2Db.lt',s=startsol,IRS=pt,DS='-',DSMIN=1E-4,DSMAX= 1e-2,ICP = [4,30,13,14,28,20,26,27,31],UZSTOP={30:[0.0,24.0],4:[3.0,30.0]},UZR={},ap=fnam)

     for x in [8]: # NS
        pt='UZ%d'%(x)
        print(linecol + '\n Following NS in N and lux ' + pt + '\x1b[0m')    
        rn(c='fjk2Db.lt',s=startsol, IRS = pt, ICP =  [4,30,13,14,20,26,21,27,31],DSMIN = 1E-4,DSMAX = 1e-2,UZSTOP={30:[0.0,24.0],4:[3.0,30.0]},UZR={},ap=fnam)
        rn(c='fjk2Db.lt',s=startsol,DS='-',IRS=pt,ICP=[4,30,13,14,20,26,21,27,31],DSMIN = 1E-4,DSMAX = 1e-2,UZSTOP={30:[0.0,24.0],4:[3.0,30.0]},UZR={},ap=fnam)


#=========================================================
def bif_tN_ext():   

     startsol = 'bif_lux_taux_ext150' 
     fnam2 = 'bif_taux3-30_N_lux150'
     dl(fnam2) 
     
     print(linecol + '\n Finding starting points' + '\x1b[0m')    
     data=sl(startsol)
     UZdat=data('UZ') 
     print UZdat
     lenUZ=len(UZdat)
     s=data()

     print(headcol + '\nFollow folds in taux and N' + '\x1b[0m')
     
     for x in range(1, lenUZ+1):  # all of them, stops at MX
        pt='UZ%d'%(x)
        sol=s(pt)
        #  20:'L1r',21:'L2r',26:'L1i',27:'L2i',
        L1rval = sol.PAR(20) # extract eigenvalues
        L1ival = sol.PAR(26) # extract eigenvalues
        L2rval = sol.PAR(21) # extract eigenvalues
        L2ival = sol.PAR(27) # extract eigenvalues
	L1sq = round(math.sqrt(L1rval*L1rval + L1ival*L1ival),2)
	L2sq = round(math.sqrt(L2rval*L2rval + L2ival*L2ival),2)
    print(headcol + '\nChecking...' + '\x1b[0m')
	print L1sq, L2sq
	#wait()

        fnam = '%s_%s'%(fnam2,pt) # Note that this must be cleared manually

	if L1sq==1:
	    print(linecol + '\n Following Fold 1 ' + pt + '\x1b[0m')    
            rn(c='fjk2Db.lt',s=startsol, IRS = pt, DSMIN = 1E-4,DSMAX = 1e-2,ICP = [4,30,13,14,28,26,21,27,31],UZSTOP={30:[-0.1,0.0,24.0,24.1],4:[3.0,30.0]},UZR={},ap=fnam2)
            rn(c='fjk2Db.lt',s=startsol,IRS=pt,DS='-',DSMIN=1E-4,DSMAX= 1e-2,ICP = [4,30,13,14,28,26,21,27,31],UZSTOP={30:[-0.1,0.0,24.0,24.1],4:[3.0,30.0]},UZR={},ap=fnam2)
	elif L2sq==1:
           print(linecol + '\n Following Fold 2 ' + pt + '\x1b[0m')    
           rn(c='fjk2Db.lt',s=startsol, IRS = pt, NMX=10000, DSMIN = 1E-4,DSMAX = 1e-2,ICP = [4,30,13,14,28,20,26,27,31],UZSTOP={30:[-0.1,0.0,24.0,24.1],4:[3.0,30.0]},UZR={},ap=fnam2)
           rn(c='fjk2Db.lt',s=startsol,IRS=pt, NMX=10000, DS='-',DSMIN=1E-4,DSMAX= 1e-2,ICP = [4,30,13,14,28,20,26,27,31],UZSTOP={30:[-0.1,0.0,24.0,24.1],4:[3.0,30.0]},UZR={},ap=fnam2)

 
#=========================================================
def bif_tN_ns():

    print(headcol + '\nLux-taux bif diag - following NS'  + '\x1b[0m')
    startsol = 'bif_taux_l1000' # point 5 is NS
    fnam = 'bif_taux_N_lux1000'

    for x in [1,5]: # NS
        pt='UZ%d'%(x)
        print(linecol + '\n Following NS in N and taux ' + pt + '\x1b[0m')    
        rn(c='fjk2Db.lt',s=startsol, NMX=1000,IRS = pt, ICP =  [4,30,13,14,20,26,21,27,31],DSMIN = 1E-4,DSMAX = 1e-2,UZSTOP={30:[0.0,24.0],4:[3.0,30.0]},UZR={},ap=fnam)
        rn(c='fjk2Db.lt',s=startsol,DS='-',NMX=1000,IRS=pt,ICP=[4,30,13,14,20,26,21,27,31],DSMIN = 1E-4,DSMAX = 1e-2,UZSTOP={30:[0.0,24.0],4:[3.0,30.0]},UZR={},ap=fnam)



#=========================================================
#=========================================================
# 3D bif
#=========================================================
#=========================================================
def bif_3D_NSstart(uzvals):

    print(headcol + '\nIts going to be a 3D bifurcation diagram dude!'  + '\x1b[0m')
    
    INt_start = 'NSstart'
    dl(INt_start)
    # Continue in I and taux to find the starting points in N
    
    print(headcol + '\nFollowing NS in I and taux to get start points'  + '\x1b[0m')
    snam = 'bif_taux_l1000' # point 5 is NS
     
    for x in [1,5]:
        pt='UZ%d'%(x)
        print(linecol + '\n' + pt + '\x1b[0m')   
        rn(c='fjk2Db.lt',s=snam, IRS = pt, ICP =  [2,4,13,14,20,26,21,27,31],DSMIN=1e-3, DSMAX=1e0,UZR={2:uzvals},UZSTOP={},ap=INt_start)
    rl(INt_start)
    
#=========================================================
def bif_3D_NS(uzvals):    

    INt_start = 'NSstart'
    # Now identify all the uz points
    print(linecol + '\n Finding starting points' + '\x1b[0m')    
    data=sl(INt_start)
    UZdat=data('UZ') 
    print UZdat
    lenUZ=len(UZdat)
    s=data()
     
    for x in range(1, lenUZ+1):  # all of them, stops at MX
        pt='UZ%d'%(x)
        sol=s(pt)
        spar = sol.PAR(2) # extract I value
        spartxt = str("%.2f" % spar).replace(".", "pt") # for non integer I
        orbsave = 'bif_INt_NS0_I%s_%s'%(spartxt,pt) # Note that this must be cleared manually

        print(endcol + '\n Following ' + orbsave + '\x1b[0m')    
        rn(c='fjk2Db.lt',s=INt_start,NPR=0,ILP=1, NMX=1000,IRS = pt, ICP =  [4,30,13,14,20,26,21,27,31],DSMIN = 1E-3,DSMAX = 1e-1,UZSTOP={30:[-0.5,0.0,24.0,24.05],4:[21.5,0.0,30.0]},UZR={},sv=orbsave)
        rn(c='fjk2Db.lt',s=INt_start,DS='-',NPR=0,ILP=1,NMX=1000,IRS=pt,ICP=[4,30,13,14,20,26,21,27,31],DSMIN = 1E-3,DSMAX = 1e-1,UZSTOP={30:[-0.5,0.0,24.0,24.05],4:[21.5,0.0,30.0]},UZR={},ap=orbsave)


#=========================================================
def bif_3D_f1start(uzvals):

    dl('F1startStar')

    # starboard
    print(linecol + '\n Fold 1: find starboard start points from I bif ' + '\x1b[0m')    
    startsol = 'bif_lux_start'# UZ 6 is f1
    pt='UZ6'
    rn(c='fjk2Db.lt',s=startsol, IRS = pt, DSMIN = 1E-3,DSMAX = 1e0,ICP =  [2,30,13,14,28,26,21,27,31],UZSTOP={30:23.0},UZR = {},sv='1')
    rn(c='fjk2Db.lt',s='1', DS='-',DSMIN = 1E-3,DSMAX = 1e0,IRS = 'UZ1', ICP =  [2,30,13,14,28,26,21,27,31],UZSTOP={30:3.0 },UZR = {30:[4.0,22.0]},sv='2')

    for x in [1,2]:
        rn(c='fjk2Db.lt',s='2',DS='-', IRS = 'UZ%d'%(x), DSMIN = 1E-3,DSMAX = 1e0,ICP =  [2,4,13,14,28,26,21,27,31],UZSTOP={2:1000},UZR = {},sv='3')
        rn(c='fjk2Db.lt',s='3',DSMIN = 1E-3,DSMAX = 1e0,IRS = 'UZ1', ICP =  [2,4,13,14,28,26,21,27,31],UZSTOP={4:24.01,2:1e-3 },UZR = {2:uzvals},ap='F1startStar')
    wait()

    rl('F1startStar')

    dl('F1startPort')
    # port
    print(linecol + '\n Fold 1: find port start points from taux bif ' + '\x1b[0m')    
    snam = 'bif_taux_l150' # points UZ5 is taux<24.2
    rn(c='fjk2Db.lt',s=snam, IRS = 'UZ5',DS='-',ICP =  [4,30,13,14,28,26,21,27,31], UZSTOP={4:[23.8]},sv='start')
    
    # continue in I and N
    rn(c='fjk2Db.lt',s='start', IRS = 'UZ1', DSMIN = 1E-3,DSMAX = 1e0,ICP =  [2,30,13,14,28,26,21,27,31],UZSTOP={2:2000.0},UZR = {2:uzvals},ap= 'F1startPort')
    rn(c='fjk2Db.lt',s='start', IRS='UZ1',DS='-',DSMIN = 1E-3,DSMAX = 1e0,ICP =  [2,30,13,14,28,26,21,27,31],UZSTOP={2: 2000.0},UZR = {2:uzvals},ap='F1startPort')

    rl('F1startPort')
    
#=========================================================
def bif_3D_f11(uzvals):

    print(headcol + '\nIts going to be a 3D bifurcation diagram dude!'  + '\x1b[0m')
    
    INt_start = 'F1startStar'
         
    # Now identify all the uz points
    print(linecol + '\n Finding starting points' + '\x1b[0m')    
    data=sl(INt_start)
    UZdat=data('UZ') 
    print UZdat
    lenUZ=len(UZdat)
    s=data()
     
    for x in range(1, lenUZ+1):  # all of them, stops at MX
        pt='UZ%d'%(x)
        sol=s(pt)
        spar = sol.PAR(2) # extract I value
        spartxt = str("%.2f" % spar).replace(".", "pt") # for non integer I
        orbsave = 'bif_INt_F11_I%s_%s'%(spartxt,pt) # Note that this must be cleared manually

        print(endcol + '\n Following Starboard ' + orbsave + ' ' + pt + '\x1b[0m')    
        rn(c='fjk2Db.lt',s=INt_start,IRS = pt, NPR=0,ICP =  [4,30,13,14,28,26,21,27,31],DSMIN = 1E-4,DSMAX = 1e-2,NMX=4000,UZSTOP={30:[-0.5,0.0,12.0,24.05],4:[23.9,30.0]},UZR={},ILP=1,sv=orbsave)
        rn(c='fjk2Db.lt',s=INt_start,DS='-',IRS=pt,NPR=0,ICP=[4,30,13,14,28,26,21,27,31],DSMIN =1E-4,DSMAX = 1e-2,NMX=4000,UZSTOP={30:[-0.5,0.0,12.0,24.05],4:[23.9,30.0]},UZR={},ILP=1,ap=orbsave)

#=========================================================
def bif_3D_f12(uzvals):

    print(headcol + '\nIts going to be a 3D bifurcation diagram dude!'  + '\x1b[0m')
    
    INt_start = 'F1startPort'
         
    # Now identify all the uz points
    print(linecol + '\n Finding starting points' + '\x1b[0m')    
    data=sl(INt_start)
    UZdat=data('UZ') 
    print UZdat
    lenUZ=len(UZdat)
    s=data()
     
    for x in range(1, lenUZ+1):  # all of them, stops at MX
        pt='UZ%d'%(x)
        sol=s(pt)
        spar = sol.PAR(2) # extract I value
        spartxt = str("%.2f" % spar).replace(".", "pt") # for non integer I
        orbsave = 'bif_INt_F12_I%s_%s'%(spartxt,pt) # Note that this must be cleared manually

        print(endcol + '\n Following Port ' + orbsave + ' ' + pt + '\x1b[0m')    
        rn(c='fjk2Db.lt',s=INt_start,IRS = pt,NPR=0, ICP =  [4,30,13,14,28,26,21,27,31],DSMIN = 1E-4,DSMAX = 1e-2,NMX=4000,UZSTOP={30:[-0.5,0.0,12.0,24.05],4:[24.5]},UZR={},ILP=1,sv=orbsave)
        rn(c='fjk2Db.lt',s=INt_start,DS='-',IRS=pt,NPR=0,ICP=[4,30,13,14,28,26,21,27,31],DSMIN= 1E-4,DSMAX = 1e-2,NMX=4000,UZSTOP={30:[-0.5,0.0,12.0,24.05],4:[24.5]},UZR={},ILP=1,ap=orbsave)

#=========================================================
def bif_3D_f2start(uzvals):

    #Get f2 start from IN bif for tau=24.2  
    print(linecol + '\n Following Fold 2 from I bif ' + '\x1b[0m')    
    startsol = 'bif_lux_start'# UZ 3 is f2
    
    dl('F2start')

    pt='UZ3'
    rn(c='fjk2Db.lt',s=startsol, IRS = pt, DSMIN = 1E-3,DSMAX = 1e0,ICP =  [2,30,13,14,28,20,26,27,31],UZSTOP={2:2000.0},UZR = {2:uzvals},sv = 'F2start')
    rn(c='fjk2Db.lt',s=startsol, DS='-',DSMIN = 1E-3,DSMAX = 1e0,IRS = pt, ICP =  [2,30,13,14,28,20,26,27,31],UZSTOP={2: 2000.0},UZR = {2:uzvals},ap='F2start')

    rl('F2start')
    
#=========================================================
def bif_3D_f2(uzvals):
    print(headcol + '\nIts going to be a 3D bifurcation diagram dude!'  + '\x1b[0m')
    
    INt_start = 'F2start'
            
    # Now identify all the uz points
    print(linecol + '\n Finding starting points' + '\x1b[0m')    
    data=sl(INt_start)
    UZdat=data('UZ') 
    print UZdat
    lenUZ=len(UZdat)
    s=data()
     
    for x in range(90, lenUZ+1):  # all of them, stops at MX
        pt='UZ%d'%(x)
        sol=s(pt)
        spar = sol.PAR(2) # extract I value
        spartxt = str("%.2f" % spar).replace(".", "pt") # for non integer I
        orbsave = 'bif_INt_F20_I%s_%s'%(spartxt,pt) # Note that this must be cleared manually

        print(endcol + '\n Following ' + orbsave + ' ' + pt + '\x1b[0m')           
        rn(c='fjk2Db.lt',s=INt_start,IRS = pt,ILP=1, NPR=0,ICP =  [4,30,13,14,28,20,26,27,31],DSMIN = 1E-3,DSMAX = 1e0,NMX=10000,UZSTOP={30:[-0.5,0.0,12.0,24.0],4:[0.0,30.0]},UZR={},sv=orbsave)
        rn(c='fjk2Db.lt',s=INt_start,DS='-',IRS=pt,ILP=1,NPR=0,ICP=[4,30,13,14,28,20,26,27,31],DSMIN = 1E-3,DSMAX = 1e0,NMX=10000,UZSTOP={30:[-0.5,0.0,12.0,24.0],4:[0.0,30.0]},UZR={},ap=orbsave)



     
#=========================================================
#=========================================================
# MAIN PROGRAM
#=========================================================
#=========================================================


ld('fjk2Dbifurcations')


# uncomment the bits you want to compute
# some computation depend on ones that preceed it


#------------------------------
## Bifurcations in lux
#----------------------------
#bif_luxpo()       # just p.o.
bif_lux()           # with variational problem
#wait()

#------------------------------
## Bifurcations in taux 
#------------------------------
#bif_taux_l50()   
#bif_taux_l150() 
#bif_taux_l1000() 
#wait()

#------------------------------
## Bifurcations in N 
#--------------------------
#bif_N_l50() 
#bif_N_l150()
#bif_N_l1000()
#wait()

#------------------------------
# 2 Parameter bifs
#------------------------------
# lux and taux
#bif_lt()

# lux and N
#bif_IN()

# taux and N
#bif_tN()
#bif_tN_ext()
#bif_tN_ns()
#wait()

#------------------------------
# 3D bif
#------------------------------
#cl()
#uzvals = np.arange(0.0, 610.0, 10.0).tolist()
#print uzvals

#bif_3D_NSstart(uzvals)
#bif_3D_NS(uzvals)
#wait()

#bif_3D_f1start(uzvals)
#bif_3D_f11(uzvals)
#bif_3D_f12(uzvals)
#wait()

#bif_3D_f2start(uzvals)
#bif_3D_f2(uzvals)
#wait()

cl()



