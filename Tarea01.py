import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from astropy import constants as cons

lambda1= np.loadtxt('sun_AM0.dat', usecols = [0])
flux1= np.loadtxt('sun_AM0.dat', usecols = [1])
lambda2=lambda1*10#[A^-1]
flux2=flux1*100 #[ergs*s^-1*cm^-2*A^-1]

semilogx(lambda2,flux2)
xlabel('$Longitud\; de\; onda\; (\lambda) \;[A^{-1}]$')
ylabel('$Flujo\; [ergs\cdot s^{-1} \cdot cm^{-2} \cdot A^{-1}]$')
title('Espectro solar')
grid(True)
savefig("EspectroSolar.png")
show()

#---------------------------Parte 2-------------------#

def intmed(x,y):
    s=0
    for i in range(len(x)-1):
        s=s+(x[i+1]-x[i])*y[i]
    return s

Lum1=intmed(lambda2,flux2)
print Lum1


#----------------------------parte 3---------------------#

def BT(x):
    T=5777
    h=cons.h.cgs.value
    c=cons.c.cgs.value
    k=cons.k_B.cgs.value
    c1=2*pi*h*c**2
    c2=(h*c)/k*T
    return (c1/x**5)/(e**(c2/x)-1)

def intmed2(x):
    s=0
    for i in range(len(x)-1):
        s=s+(x[i+1]-x[i])*BT(x[i])
    return s

T=5777
h=cons.h.cgs.value
c=cons.c.cgs.value
k=cons.k_B.cgs.value
c3=(2*pi*h)/c**2
c4=(k*T)/h




E1=intmed2(lambda2)
E2=(pi**4/15)*c3*(c4**4)
print E1
print E2
