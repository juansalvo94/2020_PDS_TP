# -*- coding: utf-8 -*-
"""
Created on Tue Aug 25 20:03:30 2020

@author: Juan Diego
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.signal as sig

def mi_funcion_sen(Vmax,Vdc,Frec,Phase,Muestras,SamplingFrec):
    """ La funcion entrega una senoidal con los parametros que figuran al costado
        Devuelve los vectores temporales y de tension.
        tt,xx
        """
    Ts = 1/SamplingFrec
    tt=np.arange(0.0,Muestras*Ts,Ts)
    x = Vmax * np.sin( 2*np.pi*Frec*tt + Phase )+Vdc
       
    return tt,x

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

N  = 1000 # muestras
fs = 1000 # Hz
a0 = 1 # Volts
p0 = 0 # radianes
f0 = 1  # Hz ->500Hz
f1 = 5
    
tt,xx = generador_senoidal (fs, f0, N, a0, p0)
tt,xx2 = generador_dienteDeSierra (fs, f1, N, a0, p0)


xx_m = np.mean(xx)
xx_var = np.var(xx)
print("Media de la distribucion:0 -- Estimacion de la media: {:g}".format(xx_m))
print("Varianza de la distribucion: -- Estimacion de la varianza: {:g}".format(xx_var))

#xx_ac = sig.correlate(xx,xx)
#xx2_ac = sig.correlate(xx2,xx2)

plt.figure(1)
line_hdls = plt.plot(tt, xx)
line_hdls = plt.plot(tt, xx2)

plt.figure(1)
plt.title('Señal: ' + 'Senoidal' )
plt.xlabel('tiempo [segundos]')
plt.ylabel('Amplitud [V]')
plt.grid(which='both', axis='both')
    
# presentar una leyenda para cada tipo de señal
axes_hdl = plt.gca()
    
# este tipo de sintaxis es *MUY* de Python
axes_hdl.legend(line_hdls, sig_type['descripcion'], loc='upper right'  )
"""    
plt.figure(2)
plt.plot(xx_ac)
plt.plot(xx2_ac)
plt.title("Secuencia de autocorrelacion")
plt.ylabel("Autocorrelacion")
plt.show()
"""
