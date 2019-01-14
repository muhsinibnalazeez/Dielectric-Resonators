# -*- coding: utf-8 -*-

"""
@author: MuhsinIbnAlAzeez
"""

import sys,os
sys.path.insert(1,os.path.join(sys.path[0],'..'))



from Dielectric_Resonator import *
import matplotlib.pyplot as plt
import numpy as np

def func(a,H,e2):
	x=34/(a*np.sqrt(e2))
	y=(a/H)+3.45
	z=x*y
	return z

values=[(0.007,0.010,38),(0.0084,0.015,38),(0.007,0.015,38),(0.004,0.010,38),
        (0.005,0.015,66),(0.007,0.010,66)]
fExp=[4.660,3.315,3.450,5.372,2.850,3.506]
Fexp=[]
Fexp1=[]
F1=[]
F2a=[]
F2H=[]
F3=[]
Feff=[]
FFF=[4.707,3.298,3.455,5.409,2.873,3.528]
FF=[]
FF1=[]
D_H=[]
fd=[]
with open('experiment.txt','w') as file:
    for i,x in enumerate(values):
        H=x[0]
        a=x[1]/2
        e2=x[2]
        file.write(str(H)+' & '+str(x[1])+' & '+str(e2)+' & '+str(fExp[i]))
        fexp=fExp[i]
        Ff=FFF[i]
        D=Resonator(a,H,e2,1)
        Da,De=D.accurate()
        fa=round(Da.f/10**9,4)
        f1=round(D.approx_1().f/10**9,4)
        f2a=round(D.approx_2a().f/10**9,4)
        f2H=round(D.approx_2H().f/10**9,4)
        f3=round(((2*Da.f)-De.f)/10**9,4)
        feff=round(De.f/10**9,4)
        fd.append((func(a,H,e2)*10**6*a)**2)
        Fexp.append((fexp*10**9*2*a)**2)
        Fexp1.append((fexp*10**9*2*a)**2*e2)
        F1.append((f1*10**9*2*a)**2)
        F2a.append((f2a*10**9*2*a)**2)
        F2H.append((f2H*10**9*2*a)**2)
        F3.append((f3*10**9*2*a)**2)
        Feff.append((feff*10**9*2*a)**2)
        FF.append((Ff*10**9*2*a)**2)
        FF1.append((Ff*10**9*2*a)**2*e2)
        D_H.append(2*a/H)
        perr=(feff-fexp)*100/fexp
        perr1=(fa-fexp)*100/fexp
        file.write(' & '+str(feff)+' & '+str(round(perr,3))+' & '+str(round(perr1,3))+'\\\\'+'\n')

D_H2=D_H[:-1]
D_H=D_H[:-2]
F1=F1[:-2]
F2a=F2a[:-2]
F2H=F2H[:-2]
F3=F3[:-2]
Feff=Feff[:-2]
Fexp2=Fexp[:-2]
fd=fd[:-2]
FF2=FF[:-2]

plt.figure(0)
plt.plot(D_H,F1,D_H,F2a,D_H,F2H,D_H,F3,D_H,Feff,D_H,Fexp2,D_H,FF2)
plt.xlabel(r'$(D/H)$')
plt.ylabel(r'$(FD)^2((Hz\cdot m)^{2})$')
plt.show()

plt.figure(2)
plt.plot(D_H2,Fexp1[:-1],'--',label=r'$f_{exp}$')
plt.plot(D_H2,FF1[:-1],'v',label=r'$fg$')

plt.figure(1)
plt.plot(D_H,fd)
plt.show()
        
        
