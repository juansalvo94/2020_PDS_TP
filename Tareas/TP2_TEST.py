# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 20:42:38 2020

@author: Juan Diego
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.signal as sig
import scipy as sci
import scipy.signal.windows as win


def generador_senoidal (fs, f0, N, a0=1, p0=0):
    """ 
    
    brief:  Generador de señales senoidal, con argumentos
    
    fs:     frecuencia de muestreo de la señal [Hz]
    N:      cantidad de muestras de la señal a generar
    f0:     frecuencia de la senoidal [Hz]
    a0:     amplitud pico de la señal [V]
    p0:     fase de la señal sinusoidal [rad]
    
    como resultado la señal devuelve:
    
    signal: senoidal evaluada en cada instante 
    tt:     base de tiempo de la señal
    """    

    # comienzo de la función  
    Ts = 1/fs
    tt=np.arange(0.0,N*Ts,Ts)        
    signal = a0 * np.sin( 2*np.pi*f0*tt + p0 )  
    # fin de la función  
    return tt, signal
def generador_cuadrada (fs, f0, N, a0=1, p0=0,dt = 0.5):
    """
    Parameters
    ----------
    fs : float
        frecuencia de sampleo.
    f0 : float
        frecuencia dela señal.
    N : int
        cantidad de muestras.
    a0 : float, optional
        amplitud. The default is 1.
    p0 : float, optional
        fase. The default is 0.
    dt: float
        duty de la señal. The default is 0.5
        
    Returns
    -------
    signal: senoidal evaluada en cada instante 
    tt:     base de tiempo de la señal

    """
    Ts = 1/fs
    tt=np.arange(0.0,N*Ts,Ts)
    signal = a0 * sig.square( 2*np.pi*f0*tt + p0,dt )
    
    return tt,signal    
def generador_dienteDeSierra (fs, f0, N, a0=1, p0=0):
    """
    Parameters
    ----------
    fs : float
        frecuencia de sampleo.
    f0 : float
        frecuencia dela señal.
    N : int
        cantidad de muestras.
    a0 : float, optional
        amplitud. The default is 1.
    p0 : float, optional
        fase. The default is 0.

    Returns
    -------
    tt : TYPE
        base de tiempo de la señal.
    signal : TYPE
        senoidal evaluada en cada instante .

    """
    Ts = 1/fs
    tt=np.arange(0.0,N*Ts,Ts)
    signal = a0* sig.sawtooth(2*np.pi*f0*tt + p0)
    
    
    return tt,signal
def mi_analizador (signal, muestras, plotSi = False, dbSi = False):
    """
    Funcion para analizar el espectro de una señal

    Parameters
    ----------
    signal : Array of float64
        Señal que quiero analizar
    muestras : int
        Cantidad de muestras que tiene la señal
    plotSi: bool
        Define si plotea la funcion de analizador  o no
    dbSi: bool
        Define si la informacion se entrega en veces o en dB

    Returns
    -------
    fft_abs : Array of float64
        vector modulo de la fft normalizado para cumplir parseval
    fft_phs : Array of float64
        Vector de fase de la fft
    eje : Array 
        Eje de frecuencias.

    """
    tam = signal.size
    fft = sci.fft.fft(signal)
    fft_abs = np.abs(fft)*2/muestras #Se utilzia 2/N para respetar parseval y que se mantenga la PSD.
    fft_abs = fft_abs[0:tam//2]
    fft_phs = np.angle(fft,deg = True)
    fft_phs = fft_phs[0:tam//2]
    
    eje=np.arange(0,muestras/2,(muestras)/tam)
    
    if dbSi:
        fft_abs = 20*np.log10(fft_abs[0:muestras//2])

    if plotSi:
        plt.plot(eje, fft_abs)
        plt.title("abs(scipy.fft.fft)")
        plt.grid(which='both', axis='both')
        plt.show()
        
    return fft_abs,fft_phs,eje

#Punto 1.a

N  = 1000 # muestras
Nz = 19*N
fs = 1000 # Hz
a0 = 1 # Volts
p0 = 0 # radianes
f0 = 1  # Hz ->500

winBartlett =   1.001002003*win.bartlett(N) 
winHann =       1.001001001*win.hann(N)
winBlackman =   1.191667858*win.blackman(N)
#winFlatTop =    2.321661568*win.flattop(N)
winFlatTop =    win.flattop(N)
fftFlatTop = sci.fft.fft(winFlatTop)
print(2*np.abs(fftFlatTop[0])/N)
winFlatTop = win.flattop(N)/(2*np.abs(fftFlatTop[0])/N)

tt=np.arange(0,N)

plt.figure(1,(21,9))
plt.plot(tt, winBartlett,label='Bartlet')
plt.plot(tt, winHann, label = 'Hann')
plt.plot(tt, winBlackman, label = 'Blackman')
plt.plot(tt, winFlatTop, label = 'Flat-Top')
plt.title('Ventanas: w(k)' )
plt.xlabel('Muestras')
plt.ylabel('Amplitud')
plt.grid(which='both', axis='both')
plt.legend()

#Se aplica Zero-Padding para mejorar la resolucion de la fft
zeros = np.zeros(Nz)
winBartlettZero = np.concatenate((winBartlett,zeros))
winHannZero = np.concatenate((winHann,zeros))
winBlackmanZero = np.concatenate((winBlackman,zeros))
winFlatTopZero = np.concatenate((winFlatTop,zeros))
#buscar la forma de  nivelar todas las fft de las ventanas

fftBartlett = mi_analizador(winBartlettZero,N,False,False)
fftHann = mi_analizador(winHannZero,N,False,False)
fftBlackman = mi_analizador(winBlackmanZero,N,False,False)
fftFlatTop = mi_analizador(winFlatTopZero,N,False,False)

plt.figure(2,(21,9))
plt.plot(fftBartlett[2], fftBartlett[0],label='Bartlet')
plt.plot(fftHann[2], fftHann[0], label = 'Hann')
plt.plot(fftBlackman[2], fftBlackman[0], label = 'Blackman')
plt.plot(fftFlatTop[2], fftFlatTop[0], label = 'Flat-Top')
plt.title('FFT de las ventanas' )
plt.xlabel('Muestras')
plt.ylabel('Amplitud')
plt.grid(which='both', axis='both')
plt.legend()


#Punto 2
"""
N  = 1000 # muestras
fs = 1000 # Hz
p0 = 0 # radianes
f1 = 250  # Hz ->500
f2 = 260
a1 = 1
a2 = 0.01
tt,xx1 = generador_senoidal(fs,f1,N,a1,p0)
tt,xx2 = generador_senoidal(fs,f2,N,a2,p0)
xx0 = xx1 + xx2

#Aplico Zero Padding
zeros = np.zeros(19*N)
#xx1 = np.concatenate((xx1,zeros))
#xx2 = np.concatenate((xx2,zeros))
#xx0 = np.concatenate((xx0,zeros))

ffta1 = np.abs(sci.fft.fft(xx1)*2/N)
ffta1 = ffta1[0:xx1.size//2]
ffta1 = 20*np.log10(ffta1)
ffta2 = np.abs(sci.fft.fft(xx2)*2/N)
ffta2 = ffta2[0:xx2.size//2]
ffta2 = 20*np.log10(ffta2)
ffta0 = np.abs(sci.fft.fft(xx0)*2/N)
ffta0 = ffta0[0:xx0.size//2]
ffta0 = 20*np.log10(ffta0)
fft_eje = np.arange(0,N/2,(N/2)/ffta1.size)
    
plt.figure(2)
plt.plot(fft_eje,ffta1,label='a1')
plt.plot(fft_eje,ffta2,label='a2')
plt.title('Señales separadas' )
plt.xlabel('Muestras')
plt.ylabel('Amplitud [dB]')
plt.grid(which='both', axis='both')
plt.legend()

plt.figure(3)
plt.plot(fft_eje,ffta0,label='a1+a2')
plt.title('Señales Juntas' )
plt.xlabel('Muestras')
plt.ylabel('Amplitud [dB]')
plt.grid(which='both', axis='both')
plt.legend()
"""