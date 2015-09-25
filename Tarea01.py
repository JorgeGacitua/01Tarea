import time
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
from astropy import constants as cons
from scipy import integrate

lambda1= np.loadtxt('sun_AM0.dat', usecols = [0])
flux1= np.loadtxt('sun_AM0.dat', usecols = [1])
lambda2=lambda1*10#[A]
flux2=flux1*100 #[ergs*s^-1*cm^-2*A^-1]

semilogx(lambda2,flux2)
xlabel('$Longitud\; de\; onda\; (\lambda) \;[\AA$')
ylabel('$Flujo\; [ergs\cdot s^{-1} \cdot cm^{-2} \cdot \AA^{-1}]$')
title('Espectro solar')
grid(True)
savefig("EspectroSolar.png")
show()

#---------------------------Parte 2-------------------#
au=cons.au.value
def intrap(x,y):
    s=0
    for i in range(len(x)-1):
        s=s+(x[i+1]-x[i])*y[i]
    return s

Lum1=4*np.pi*(au**2)*intrap(lambda1,flux1)
print "Luminosidad1=" + str(Lum1)


#----------------------------parte 3---------------------#
T=5777
h=cons.h.value
c=cons.c.value
k=cons.k_B.value


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
    '''
    Calcula la integral de la funcion 'f' entre los
    valores 'a y 'b usando el metodo del valor medio
    '''
    m=(a+b)/2
    s=(b-a)*f(m)
    return s



def simpson(f,a,b):
    '''
    Calcula la integral de la funcion 'f' entre los
    valores 'a y 'b usando el metodo de Simpson
    '''
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

radio=np.sqrt(Lum1/(4*np.pi*E1))
print "Radio="+str(radio)+" metros"
#---------------------------Parte 4------------------------#

t0=time.time()
Lum2=4*np.pi*(au**2)*np.trapz(flux1,x=lambda1)
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
