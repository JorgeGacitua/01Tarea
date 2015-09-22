import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from astropy import constants as cons
from scipy import integrate

lambda1= np.loadtxt('sun_AM0.dat', usecols = [0])
flux1= np.loadtxt('sun_AM0.dat', usecols = [1])
lambda2=lambda1*10#[A^-1]
flux2=flux1*100 #[ergs*s^-1*cm^-2*A^-1]

semilogx(lambda2,flux2)
xlabel('$Longitud\; de\; onda\; (\lambda) \;[\AA^{-1}]$')
ylabel('$Flujo\; [ergs\cdot s^{-1} \cdot cm^{-2} \cdot \AA^{-1}]$')
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
print "Luminosidad1=" + str(Lum1)


#----------------------------parte 3---------------------#
T=5777
h=cons.h.cgs.value
c=cons.c.cgs.value
k=cons.k_B.cgs.value

def BT(x):
    c1=2*pi*h*c**2
    c2=(h*c)/k*T
    return (c1/x**5)/(e**(c2/x)-1)

def intmed2(x):
    s=0
    for i in range(len(x)-1):
        s=s+(x[i+1]-x[i])*BT(x[i])
    return s

def F1(x):
    return (x**3)/(np.exp(x)-1)
def intmed3(x):
    s=0
    for i in range(len(x)-1):
        s=s+(x[i+1]-x[i])*F1(x[i])
    return s

E1=intmed2(lambda2)
E2=((2*pi*h/c**2)*(k*T/h)**4)*((pi**2)/4)
E3=intmed3(lambda2)

print "Energia1="+str(E1)
print "Energia2="+str(E2)
print "Energia3="+str(E3)

#---------------------------Parte 4------------------------#

Lum2=np.trapz(lambda2,x=flux2)
print "Luminosidad2=" + str(Lum2)
F=lambda x:(x**3)/(np.exp(x)-1)
c5=((2*pi*h/c**2)*(k*T/h)**4)
ee1=integrate.quad(F,0,np.inf)
E4=c5*ee1[0]
print "Energia4=" + str(E4)
