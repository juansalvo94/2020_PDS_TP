# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 20:12:58 2020

@author: Juan Diego
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.signal as sig
import scipy as sci
from timeit import default_timer as timer
import time

"""A esta altura del curso, ya disponemos de dos funcionalidades que nos serán
 muy útiles a lo largo del curso, un generador de señales (senoidales, ruido blanco) 
 y un analizador de espectro. Si bien el analizador está en un estado embrional, 
 aprovecharemos esta tarea para comenzar a desarrollar esta herramienta que será 
 fundamental para lo que resta de curso.
 
Revisá los videos acerca de la DFT y realizá una función que te permita visualizar 
el espectro de la señal que le ingreses. Debería tener una interfaz parecida a esta:
    
mi_analizador( xx, ...., (otros parámetros que vos consideres útiles)  )

Como resultado de la función me gustaría poder ver el módulo (o la densidad de 
potencia espectral) y su fase. Todas estas señales deberán visualizarse entre 0 y fs/2 [Hz].
A esta altura deberías contar ya con algunas funciones desarrolladas, por lo que 
espero como entrega, un enlace a tu repositorio en github y eventualmente un testbench
 donde ensayes tu función o preferentemente un notebook con alguna experimentación 
 que hayas realizado y tus observaciones.
"""

def mi_funcion_sen(Vmax,Vdc,Frec,Phase,Muestras,SamplingFrec):
    """ La funcion entrega una senoidal con los parametros que figuran al costado
        Devuelve los vectores temporales y de tension.
        tt,xx
        """
    Ts = 1/SamplingFrec
    tt=np.arange(0.0,Muestras*Ts,Ts)
    x = Vmax * np.sin( 2*np.pi*Frec*tt + Phase )+Vdc
       
    return tt,x
       

def mi_analizador (signal, muestras, plotSi, dbSi):
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
    fft = sci.fft.fft(signal)
    fft_abs = np.abs(fft)*2/muestras #Se utilzia 2/N para respetar parseval y que se mantenga la PSD.
    fft_abs = fft_abs[0:muestras//2]
    fft_phs = np.angle(fft,deg = True)
    fft_phs = fft_phs[0:muestras//2]
    eje = np.arange(muestras/2)
    
    if dbSi:
        fft_abs = 20*np.log10(fft_abs[0:muestras//2])

    if plotSi:
        plt.plot(eje, fft_abs)
        plt.title("abs(scipy.fft.fft)")
        plt.grid(which='both', axis='both')
        plt.show()
        
    return fft_abs,fft_phs,eje

N  = 10000 # muestras
fs = 1000 # Hz
a0 = 1 # Volts
p0 = 0 # radianes
f0 = 10  # Hz 
ts = 1/fs

tt,xx = mi_funcion_sen (a0,0,f0,p0,N,fs)

dens, phs, ejex = mi_analizador(xx,N, True, False)
"""
plt.figure(1)
plt.stem(ejex, dens,use_line_collection = True)
plt.title("abs(scipy.fft.fft)")
plt.grid(which='both', axis='both')
plt.show()
"""
