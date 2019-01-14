# -*- coding: utf-8 -*-

"""
@author: MuhsinIbnAlAzeez
"""

import scipy.special as spl
import math
epsilon=0.0005
iterate=10**5
eps1=0.0007
eps2=0.0001
class Resonator(object):
    c=2.9979e8
    
    def __init__(self,a,H,e2,p=1,f=None,mode="TE"):
        self.mode=mode
        if mode=='TE':
            self.Xmn=spl.jnp_zeros(0,1)[0]
        elif mode=='TM':
            self.Xmn=spl.jn_zeros(0,1)[0]
        else:
            raise ValueError("mode must be either 'TE' or 'TM'" )
        self.a=a
        self.H=H
        self.e2=e2
        self.p=p
        self.f=f
        self.err=0
    def __str__(self):
        y='%e' %self.f if self.f!=None else None
        x='''Dielectric Resonator model.
+=======================================================+
|                     parameters                        |
+==========================+============================+
|   Radius                 |  a= {}m{}|
|   Height                 |  H= {}m{}|
|   relative permitivity   |  e2= {}{}|
|   {}01{} mode frequency   |  f= {}{}|
+==========================+============================+
    '''.format(self.a,' '*(22-len(str(self.a))),self.H,' '*(22-len(str(self.H))),self.e2,' '*(22-len(str(self.e2))),self.mode,self.p,y,' '*(23-len(str(y))))
        return x
    
    
    def approx_2H(self):
        a=self.a
        H=self.H
        kc=self.Xmn/a
        f=kc*self.c/(2*math.pi*math.sqrt(self.e2))+1000
        #f_max=kc*self.c/(2*math.pi)-1000
        def X(f):
            k0=2*math.pi*f/self.c
            B=math.sqrt((k0*k0*self.e2)-(kc*kc))
            A=math.sqrt((kc*kc)-(k0*k0))
            x=(2/B)*math.atan(A/B)
            return x-H
        
#        if X(f_min)>X(f_max):
#            f_min,f_max=f_max,f_min
        
#        f=(f_min+f_max)/2
        
        while abs(X(f))>=epsilon:
            f+=iterate
#            if X(f)>0:
#                f_max=f
#            else:
#                f_min=f
#            f=(f_min+f_max)/2
        print('the Frequency =','%e' %f)
        print('x: ',abs(X(f))+H)
        self.err+=abs(X(f))
        self.err/=2
        if self.mode=='TE':
            return Resonator(a,H,self.e2,self.p,f)
        else:
            return Resonator(a,H,self.e2,self.p,f,mode='TM')
    
    
    def approx_2a(self):
        a=self.a
        H=self.H
        b=self.p*math.pi/H
        f_2a=b*self.c/(2*math.pi*math.sqrt(self.e2))+1000
        #f_max=f_min
        #f_max=(b*c/(2*math.pi))-1000
        def X(f1a):
            k0=2*math.pi*f1a/self.c
            kco=math.sqrt(b*b-k0*k0)
            kci=math.sqrt(k0*k0*self.e2-b*b)
            pi=(kci*a)
            po=(kco*a)
            j0=spl.jn(0,pi)
            dj0=-spl.jn(1,pi)
            K0=spl.kn(0,po)
            dK0=-spl.kn(1,po)
            if self.mode=='TE':
                u=(pi*j0/dj0)
                v=-(po*K0/dK0)
            else:
                u=(pi*j0/dj0)
                v=-(self.e2*po*K0/dK0)
            return u,v
        u=X(f_2a)[0]
        v=X(f_2a)[1]
        while u<=v:
            f_2a+=iterate
            u=X(f_2a)[0]
            v=X(f_2a)[1]
        
#        while X(f_max)<0:
#            f_max+=1e8
#        
#        f1a=(f_min+f_max)/2
#        
#        
#        while abs(X(f1a))>=epsilon:
#            if X(f1a)>0:
#                f_max=f1a
#            else:
#                f_min=f1a
#            f1a=(f_min+f_max)/2
        
        print('the Frequency =','%e' %f_2a)
        print('x: ',abs(X(f_2a)[0]-X(f_2a)[1]))
        self.err+=abs(X(f_2a)[0]-X(f_2a)[1])
        self.err/=2
        if self.mode=='TE':
            return Resonator(a,H,self.e2,self.p,f_2a)
        else:
            return Resonator(a,H,self.e2,self.p,f_2a,mode='TM')
    
    
    def approx_1(self,label='f'):
        '''
        Approximation 1
        flags:
            a : find a from given f and H
            H : find H from given f and a
            f : find f from given a and H
        '''
        a=self.a
        H=self.H
        f=self.f
        def f_1(a,H):
            return (self.c/(2*math.pi*math.sqrt(self.e2)))*math.sqrt((self.Xmn/a)**2+(self.p*math.pi/H)**2)
        if label=='H':
            f_1H=f+1000
            while f_1H>=f:
                f_1H=f_1(a,H)
                H+=eps1
            print('height','%.4f' %H)
            print('error: ',abs(f-f_1(a,H))/f)
            self.err+=abs(f-f_1(a,H))/f
            self.err/=2
            if self.mode=='TE':
                return Resonator(a,H,self.e2,self.p,f)
            else:
                return Resonator(a,H,self.e2,self.p,f,mode='TM')
        elif label=='a':
            if self.mode=='TE':
                f_1H=f+1000
                while f_1H>=f:
                    f_1H=f_1(a,H)
                    a+=eps2
            else:
                f_1H=f-1000
                while f_1H<=f:
                    f_1H=f_1(a,H)
                    a-=eps2
                
            print('radius','%.5f' %a)
            print('error: ',abs(f-f_1(a,H))/f)
            self.err+=abs(f-f_1(a,H))/f
            self.err/=2
            if self.mode=='TE':
                return Resonator(a,H,self.e2,self.p,f)
            else:
                return Resonator(a,H,self.e2,self.p,f,mode='TM')
        elif label=='f':
            f=f_1(a,H)
            print('frequency:','%e' %f)
            if self.mode=='TE':
                return Resonator(a,H,self.e2,self.p,f)
            else:
                return Resonator(a,H,self.e2,self.p,f,mode='TM')
        else:
            raise ValueError("label must be either 'a','H' or 'f'")
    
    
    def approx_3(self):
        a=self.a
        H=self.H
        b=self.p*math.pi/H
        f_3=b*self.c/(2*math.pi*math.sqrt(self.e2))+1000
        #f_max=f_min
        #f_max=(b*c/(2*math.pi))-1000
        
        def X(f1a):
            k0=2*math.pi*f1a/self.c
            kco=math.sqrt(b*b-k0*k0)
            kci=math.sqrt(k0*k0*self.e2-b*b)
            pi=(kci*a)
            po=(kco*a)
            j0=spl.jn(0,pi)
            dj0=-spl.jn(1,pi)
            K0=spl.kn(0,po)
            dK0=-spl.kn(1,po)
            if self.mode=='TE':
                u=(pi*j0/dj0)
                v=-(po*K0/dK0)
            else:
                u=(pi*j0/dj0)
                v=-(self.e2*po*K0/dK0)
            return u,v,kci
        
#        while X(f_max)[0]<0:
#            f_max+=1e7
#        
#        f1a=(f_min+f_max)/2
        
        
        while X(f_3)[0]<=X(f_3)[1]:
            f_3+=iterate
#            if X(f1a)[0]>0:
#                f_max=f1a
#            else:
#                f_min=f1a
#            f1a=(f_min+f_max)/2
        
        
        
        print('the Frequency =',f_3)
        print('x: ',abs(X(f_3)[0]-X(f_3)[1]))
        self.err+=abs(X(f_3)[0]-X(f_3)[1])
        self.err/=2
        kci=X(f_3)[-1]
        print(kci)
        #For the next iteration
        F=math.sqrt(kci*kci*self.c*self.c/(self.e2*4*math.pi**2))+1000
        #f_max=kci*self.c/(2*math.pi)-1000
        def X(f):
            '''
            helper function
            '''
            k1=2*math.pi*f/self.c
            B=math.sqrt(k1*k1*self.e2-kci*kci)
            u=(2*math.pi*f/self.c)**2*(self.e2-1)
            v=(B*B)*(1+(math.tan(H*B/2))**2)
            return u,v
        u=X(F)[0]
        v=X(F)[1]
        while u>=v:
            F+=iterate
            k1=2*math.pi*F/self.c
            B=math.sqrt(k1*k1*self.e2-kci*kci)
            u=(2*math.pi*F/self.c)**2*(self.e2-1)
            v=(B*B)*(1+(math.tan(H*B/2))**2)
        
        #printing the result
        print('frequency fi:','%e' %F)
        print('Error: ',abs(X(F)[0]-X(F)[1]))
        self.err+=abs(X(F)[0]-X(F)[1])
        self.err/=2
        if self.mode=='TE':
            return Resonator(a,H,self.e2,self.p,F)
        else:
            return Resonator(a,H,self.e2,self.p,F,mode='TM')
    
    def accurate(self):
        '''Help to accurately evaluate the resonant frequency
of Dielectric resonator.

Returns: tuple (final DR, effective DR)'''
        def dprint(x):
            print('='*20)
            print(x)
            print('='*20)
        
        
        #executing
        dprint('(a,H)+2H->f1H')
        D1=self.approx_2H()
        dprint('(f1H,a)+1->H1')
        D2=D1.approx_1('H')
        dprint('(a,H1)+2a->f1a')
        D3=D2.approx_2a()
        dprint('(f1a,H1)+1->ae')
        D4=D3.approx_1('a')
        ae=D4.a
        dprint('(a,H)+2a->f2a')
        D5=self.approx_2a()
        dprint('(f2a,H)+1->a1')
        D6=D5.approx_1('a')
        dprint('(a1,H)+2H->f2H')
        D7=D6.approx_2H()
        dprint('(f2H,a1)+1->He')
        D8=D7.approx_1('H')
        He=D8.H
        dprint('(ae,He)+1->fe')
        if self.mode=='TE':
            De=Resonator(ae,He,self.e2,self.p)
        else:
            De=Resonator(ae,He,self.e2,self.p,mode='TM')
        De=De.approx_1()
        dprint('(a,H)+3->fi')
        Di=self.approx_3()
        dprint('F=(fe+fi)/2')
        F=(De.f+Di.f)/2
        print('#'*20)
        print('the frequency F =','%e' %F)
        print('error:',self.err)
        if self.mode=='TE':
            return Resonator(self.a,self.H,self.e2,self.p,F),De
        else:
            return Resonator(self.a,self.H,self.e2,self.p,F,mode='TM'),De
