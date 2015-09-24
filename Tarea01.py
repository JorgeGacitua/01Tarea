import time
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

def intrap(x,y):
    s=0
    for i in range(len(x)-1):
        s=s+(x[i+1]-x[i])*y[i]
    return s

Lum1=intrap(lambda2,flux2)
print "Luminosidad1=" + str(Lum1)


#----------------------------parte 3---------------------#
T=5777
h=cons.h.cgs.value
c=cons.c.cgs.value
k=cons.k_B.cgs.value
au=cons.au.value*100 # en centimetros

def F1(x):
    f=(x**3)/(np.exp(x)-1)
    return f
def F2(x):
    f=((np.tan(x)**3)*(1+np.tan(x)**2))/(np.exp(np.tan(x))-1)
    return f

def rango(a,b,step):
    '''
    Crea un vector con valores entre 'a' y 'b' con un espaciasdo 'step'
    '''
    while a<=b:
            yield a
            a+=step

def intmed(f,a,b):
    m=(a+b)/2
    s=(b-a)*f(m)
    return s



def simpson(f,a,b):
    n=100000
    s=f(a)+f(b)
    h=(b-a)/n
    for i in rango(1,n-1,2):
        s=s+4*f(a+i*h)
    for j in rango(2,n-2,2):
        s=s+2*f(a+j*h)
    return s*h/3

C=((2*pi*h/c**2)*(k*T/h)**4)
t0=time.time()
E1=C*(intmed(F2,0.0,0.25)+simpson(F2,0.25,(np.pi/2)-0.25)+intmed(F2,(np.pi/2)-0.25,(np.pi/2)))
tf=time.time()
E2=C*((np.pi**4)/15) #resultado de calcular la integral analíticamente


print "Energia1="+str(E1)

print "Tiempo calculando simpson="+str(tf-t0)+" segundos"
print "Energia2="+str(E2)

radio=au*np.sqrt(E1/(4*np.pi*Lum1))
print "Radio="+str(radio)
#---------------------------Parte 4------------------------#
t0=time.time()
Lum2=np.trapz(lambda2,x=flux2)
tf=time.time()
print "Luminosidad2=" + str(Lum2)
print "Tiempo calculando trapecio="+str(tf-t0)+" segundos"


F=lambda x:(x**3)/(np.exp(x)-1)
c5=((2*pi*h/c**2)*(k*T/h)**4)
t0=time.time()
ee1=integrate.quad(F,0,np.inf)
E4=c5*ee1[0]
tf=time.time()
print "Energia4=" + str(E4)
print "Tiempo calculando quad="+str(tf-t0)+" segundos"
