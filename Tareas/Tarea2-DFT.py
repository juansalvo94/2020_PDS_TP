# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 20:03:30 2020

@author: Juan Diego
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.signal as sig
import scipy as sci
from timeit import default_timer as timer
import time

def mi_funcion_sen(Vmax,Vdc,Frec,Phase,Muestras,SamplingFrec):
    """ La funcion entrega una senoidal con los parametros que figuran al costado
        Devuelve los vectores temporales y de tension.
        tt,xx
        """
    Ts = 1/SamplingFrec
    tt=np.arange(0.0,Muestras*Ts,Ts)
    x = Vmax * np.sin( 2*np.pi*Frec*tt + Phase )+Vdc
       
    return tt,x
       
def My_DFT(signal,N):
	output = np.arange(0.0,N,dtype="complex_")
	for k in range(N):  # For each output element
		s = complex(0)
		for t in range(N):  # For each input element
			angle = 2j * np.pi * t * k / N
			s += signal[t] * np.exp(-angle)
		output[k] = s
	return output
    
""" Tarea 2.b
muestras = []
for x in range(4,7):
    muestras.append(2**x)
    
Tiempos = []
Tiempos_txt =[]
for N in muestras:  


    fs = 1000 # Hz
    a0 = 1 # Volts
    p0 = 0 # radianes
    f0 = 25.5  # Hz 
    ts = 1/fs

    tt,xx = mi_funcion_sen (a0,0,f0,p0,N,fs)  

    Tini = timer()
    dft =My_DFT( xx,N)
    dft_x = np.arange(0,N)
    dft_abs = np.abs(dft)*2/N
    dft_phs = np.angle(dft,deg = True)
    Tfin = timer()
    Tiempo_dft = Tfin-Tini
    Tiempos.append(Tiempo_dft)
    Tiempos_txt.append(str(Tiempo_dft))
  
    
txt = str(Tiempos[0])
print (txt)
"""
""" fft
Tini = timer()

t = np.arange(N)
fft = sci.fft.fft(xx)
fft_abs = np.abs(fft)*2/N #Se utilzia 2/N para respetar parseval y que se mantenga la PSD.
#Se podria pasar a DB para mejorar los detalles en valores minimos.
fft_abs_db = 20*np.log10(fft_abs[0:N//2])
fft_phs = np.angle(fft,deg = True)

Tfin = timer()
Tiempo_fft = Tfin-Tini

freq = np.arange(N/2)
#freq = np.fft.fftfreq(t.shape[-1])*N

plt.figure(1)
#plt.plot(freq, fft)
plt.plot(freq, fft_abs_db, "x")
plt.title("abs(scipy.fft.fft)") 
#Ver para plotear sin lineas
plt.grid(which='both', axis='both')
plt.show()
"""
""" tarea 3.a
bins = [0,0.01,0.25,0.5]

fCentral = []
pAdyacente = []
restoDeFrecuencias= []

for mod in bins:
    N = 1000
    fs = 1000 # Hz
    a0 = 1 # Volts
    p0 = 0 # radianes
    f0 = 250+mod  # Hz 
    ts = 1/fs
    
    tt,xx = mi_funcion_sen (a0,0,f0,p0,N,fs)
    
    t = np.arange(N)
    fft = sci.fft.fft(xx)
    fft_abs = np.abs(fft)*2/N #Se utilzia 2/N para respetar parseval y que se mantenga la PSD.
    #Se podria pasar a DB para mejorar los detalles en valores minimos.
    freq = np.arange(N/2)
    
    fCentral.append(str(fft_abs[250]))
    pAdyacente.append(str(fft_abs[251]))
    rdf = 0.0
    for f in freq:
        if f != 250:
            rdf += (fft_abs[int(f)])**2
    restoDeFrecuencias.append(str(rdf))
    
print (fCentral)
print (pAdyacente)
print (restoDeFrecuencias)
 """   
#Ej 3.b
N = 1000
fs = 1000 # Hz
a0 = 1 # Volts
p0 = 0 # radianes
f0 = 250+0.5  # Hz 
ts = 1/fs
tt,xx0 = mi_funcion_sen (a0,0,f0,p0,N,fs)
xx1 = xx0
xx2 = xx0
xx3 = xx0
Mj = [0,100,N,10*N]
aux = np.arange(N//10)
aux.fill(0.0)
xx1 = np.hstack((xx1,aux))

aux = np.arange(N)
aux.fill(0.0)
xx2 = np.hstack((xx2,aux))

aux = np.arange(10*N)
aux.fill(0.0)
xx3 = np.hstack((xx3,aux))

fft0 = sci.fft.fft(xx0)
fft_abs0 = np.abs(fft0)*2/N
fft_abs0 = fft_abs0[0:N//2]
fft1 = sci.fft.fft(xx1)
fft_abs1 = np.abs(fft1)*2/N
fft_abs1 = fft_abs1[0:N//2]
fft2 = sci.fft.fft(xx2)
fft_abs2 = np.abs(fft2)*2/N
fft_abs2 = fft_abs2[0:N//2]
fft3 = sci.fft.fft(xx3)
fft_abs3 = np.abs(fft3)*2/N
fft_abs3 = fft_abs3[0:N//2]
 #Se utilzia 2/N para respetar parseval y que se mantenga la PSD.
#Se podria pasar a DB para mejorar los detalles en valores minimos.
freq = np.arange(N/2)
        

plt.figure(1)
plt.plot(freq,fft_abs0)
plt.plot(freq,fft_abs1)
plt.plot(freq,fft_abs2)
plt.plot(freq,fft_abs3)
plt.title("fft's")
plt.grid(which='both', axis='both')





"""
plt.figure(2)
plt.plot(dft_x, dft_abs,dft_x,dft_phs)
#plt.plot(dft_x, dft_abs)
plt.title("abs(My_DFT)")
plt.grid(which='both', axis='both')

plt.figure (3)
plt.plot(freq, np.real(dft) - np.real(fft))
plt.title("diferencia parte real (My_DFT- fft)")
plt.grid(which='both', axis='both')  

plt.figure (4)
plt.plot(freq, np.imag(dft) - np.imag(fft))
plt.title("diferencia parte imag (My_DFT- fft)")
plt.grid(which='both', axis='both')  

plt.figure (5)
plt.plot(tt, xx)
plt.title("Señal en T")
plt.grid(which='both', axis='both') 

"""
""" Comentarios de entrega N°1:
    La DFT manual me entrega los valores en los extremos del eje X.
    Al aumentar la frecuencia, los valores que deberian alejarse del centro, 
    se acercan al mismo.
    Logicamente, ocurre lo contrario (y correcto) en el algoritmo de fft que 
    provee numpy.
    Le dedique unos minutos a esto para cranearlo, pero no hubo caso.
    
    No se si me falto trasponer la matriz de resultados o que, pero en caso de
    tener que trasponerla por el momento no sabría como.
    """