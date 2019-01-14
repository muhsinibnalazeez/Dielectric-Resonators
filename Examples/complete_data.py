# -*- coding: utf-8 -*-
"""
@author: MuhsinIbnAlAzeez
"""

import matplotlib.pyplot as plt
import numpy as np

A=[0.01/2,0.015/2,0.015/2,0.01/2,0.015/2]#,0.01/2]
F1=[6.8739,4.9007,5.2634,8.4934,4.7561,5.2158]
F2H=[6.4496,4.4300,4.5603,7.0249,3.6391,4.8646]
F2a=[5.6887,4.2417,4.7103,7.7916,4.4730,4.3223]
F3=[5.0962,3.6241,3.8357,6.0510,3.2647,3.8756]
Feff=[4.7077,3.2981,3.4547,5.4091,2.8731,3.5285]
Fexp=[4.660,3.315,3.450,5.372,2.850,3.506]
D_H=[1.4286,1.7857,2.1428,2.5000,3.0000]#,1.4286]

#F1=np.array(F1[:-1])
#F2H=np.array(F2H[:-1])
#F2a=np.array(F2a[:-1])
#F3=np.array(F3[:-1])
#Feff=np.array(Feff[:-1])
#Fexp=np.array(Fexp[:-1])
#D_H=np.array(d_h[:-1])

f1=[]
f2H=[]
f2a=[]
f3=[]
feff=[]
fexp=[]

for (i,a) in enumerate(A):
    if i!=4:
        e2=38
    else:
        e2=66
    ff1=(F1[i]*10**9*2*a)**2*e2
    ff2H=(F2H[i]*10**9*2*a)**2*e2
    ff2a=(F2a[i]*10**9*2*a)**2*e2
    ff3=(F3[i]*10**9*2*a)**2*e2
    ffeff=(Feff[i]*10**9*2*a)**2*e2
    ffexp=(Fexp[i]*10**9*2*a)**2*e2
    f1.append(ff1)
    f2H.append(ff2H)
    f2a.append(ff2a)
    f3.append(ff3)
    feff.append(ffeff)
    fexp.append(ffexp)


plt.figure(0)
plt.plot(D_H,f1,label=r'$f_{1}$')
plt.plot(D_H,f2a,label=r'$f_{2a}$')
plt.plot(D_H,f2H,label=r'$f_{2H}$')
plt.plot(D_H,f3,label=r'$f_{3}$')
plt.plot(D_H,feff,label=r'$f_{eff}$')
plt.plot(D_H,fexp,'^',label=r'$f_{exp}$')
plt.xlabel(r'$(D/H)$')
plt.ylabel(r'$(FD)^2\epsilon_{2}((Hz\cdot m)^{2})$')
plt.legend(loc='best')
plt.show()

plt.figure(1)
plt.plot(D_H,fexp,'^',label='Experimental values')
plt.plot(D_H,feff,'-',label='computed values')
plt.xlabel(r'$(D/H)$')
plt.ylabel(r'$(FD)^2\epsilon_{2}((Hz\cdot m)^{2})$')
plt.legend(loc='best')
plt.show()
    
    
    
    
    
    
    
    
