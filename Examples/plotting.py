# -*- coding: utf-8 -*-
"""
@author: MuhsinIbnAlAzeez
"""

import sys,os
sys.path.insert(1,os.path.join(sys.path[0],'..'))

from Dielectric_Resonator import *
import matplotlib.pyplot as plt


F1=[]
F2H=[]
F2a=[]
F3=[]
Feff=[]
Fexp=[]
D_H=[]

values=[(0.007,0.010,38),(0.0084,0.015,38),(0.007,0.015,38),(0.004,0.010,38),
        (0.005,0.015,66)]#,(0.007,0.010,66)]
fExp=[4.660,3.315,3.450,5.372,2.850,3.506]

with open('experiment.txt','w') as file:
    for i,x in enumerate(values):
        H=x[0]
        a=x[1]/2
        e2=x[2]
        file.write(str(H)+' & '+str(x[1])+' & '+str(e2)+' & '+str(fExp[i]))
        fexp=fExp[i]
        epsilon=0.0001 if e2==38 else 0.0005
        eps1=H/10
        eps2=2*a/100
        D=Resonator(a,H,e2,1)
        Da,De=D.accurate()
        fa=round(Da.f/10**9,4)
        f1=round(D.approx_1().f/10**9,4)
        f2a=round(D.approx_2a().f/10**9,4)
        f2H=round(D.approx_2H().f/10**9,4)
        f3=round(((2*Da.f)-De.f)/10**9,4)
        feff=round(De.f/10**9,4)
        Fexp.append((fexp*10**9*2*a)**2*e2)
        F1.append((f1*10**9*2*a)**2*e2)
        F2a.append((f2a*10**9*2*a)**2*e2)
        F2H.append((f2H*10**9*2*a)**2*e2)
        F3.append((f3*10**9*2*a)**2*e2)
        Feff.append((feff*10**9*2*a)**2*e2)
        D_H.append(2*a/H)
        file.write(' & '+str(feff)+'\\\\'+'\n')
        
        
        
plt.figure(0)
plt.plot(D_H,F1,D_H,F2a,D_H,F2H,D_H,F3,D_H,Feff,D_H,Fexp)
plt.xlabel(r'$(D/H)$')
plt.ylabel(r'$(FD)^2\epsilon_{2}((Hz\cdot m)^{2})$')
plt.show()

plt.figure(1)
plt.plot(D_H,Fexp,'--',D_H,Feff,'v')
plt.xlabel(r'$(D/H)$')
plt.ylabel(r'$(FD)^2\epsilon_{2}((Hz\cdot m)^{2})$')
plt.show()














