# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 17:48:52 2019

@author: alan_
"""
#### En este script se desarrollaron las funciones de la actividad 1
### y se implementaron con señales sinteticas (superposicion de dos senos)

import numpy as np
from matplotlib import pyplot as plt
from scipy import signal as sig
from scipy import fftpack as fft
from scipy.signal import hann
from scipy.signal import chebwin

def PSD(matrix):
    N1=np.size(matrix,1)
    N2=np.size(matrix,0)
    #rxx=sig.convolve2d(matrix,matrix,mode='same')
    rxx=matrix
    rxx=rxx*(chebwin(N1,at=100).reshape(1,N1)*np.ones((N2,1)))
    sxx=np.fft.fft2(rxx)
    mag_sxx=(np.abs(sxx[:,0:N1/2])*ts)**2
    mag_sxx=10*np.log10(np.mean(mag_sxx,0))
    F=np.linspace(0,fs/2,len(mag_sxx))
    return F,mag_sxx

def inv_fitting(t,y):
    A=[]
    A.append(t)
    A.append(np.ones(len(t)))
    A=np.array(A).T
    inv=pseudo_inv(A)
    slcion=np.dot(y.reshape(1,len(y)),inv)
    return slcion

def pseudo_inv(A):
    U, Sigma, V = np.linalg.svd(A, full_matrices=False,compute_uv=True)
    Sigma_pseu=1./Sigma
    inv=np.matrix(U)*np.diag(Sigma_pseu)*np.matrix(V)
    return inv


def barlett_par (v):
    x1=np.arange(0,v/2)
    x2=np.arange(v/2,v)
    triangular=np.append(2*x1/float(v),2-(2*x2/float(v)))
    return triangular

def MAverage(x,N):
    L=len(x)-N
    y=np.zeros(L)
    for i in range(L):
        y[i]=x[i]-((1./(N+1))*np.sum(x[i:i+N]))
    return y

def Mwindow(x,N,overlap):
    L=len(x)
    ind=0
    m=[]
    y=[]
    for i in range(L):
       if ((i+1)-ind) % N == 0:
            y.append(x[ind]-np.mean(x[ind:i+1]))
            m.append(ind)
            i=int(i-round(N*overlap))
            ind=i+1
    return y,np.array(m)*ts

def FFT (x,v1,v2,v3):
    v=[v1,v2,v3]
    lab=['hann','barlett','chebwin']
    f=np.linspace(0,fs/2,N/2)
    plt.figure(figsize=(12,5))
    for i in range(len(v)):
        Fw=20*np.log10(np.abs(fft.fft(p*v[i])))
        plt.plot(f[0:N/10],Fw[0:N/10],linewidth=4,label=lab[i])
        legend = plt.legend(loc='upper right', shadow=True, fontsize='medium')
        legend.get_frame().set_facecolor('pink') 
        plt.title('FFT con 3 ventanas')
        plt.grid()        

def RF(den,num,fm):
    w, h = sig.freqz(den,num)
    h[h==0] = 1E-5

    H    = 20*np.log10( np.abs(h) )
    W    = np.angle  (h)
    W    = np.unwrap (W)
    W    = np.degrees(W)
    w    = np.linspace(0,fs/2,H.shape[0] )

    return w, W, H

#DECLARACION DE CONSTANTES
fs=320
ts=1./fs
D=1
F=[]
F.append(12)
F.append(17.5)
#SUPERPOSICION DE SEÑALES
F=np.array(F)[np.newaxis]
t=np.arange(0,D,1./fs)[:,np.newaxis]
p=2*np.pi*t*F
p=np.sin(p)
p=p.sum(axis=1)
N=len(p)
#FFT VENTANAS
FFT(p,hann(N),barlett_par(N),chebwin(N,at=100))

#%%#    estimacion de psd con periodograma de welch

f=plt.figure(figsize=(8,5))
N=len(p)
v=np.asarray([2,7,12,27])           # numero de ventanas a usar 
plt.figure(figsize=(10,8))
for k in range (len(v)):
    x=0
    if N % v[k]!=0:
        Nzeros=v[k]-(N % v[k])
        x=np.append(p,np.zeros(Nzeros))      # agregar ceros para mejorar la estimación, y para que el reshape se pueda hacer
    else:
        x=p
    Nv=len(x)/v[k]
    matrix=x.reshape(v[k],Nv)
    F,sxx=PSD(matrix)
    plt.subplot(2,2,k+1)
    plt.plot(F,sxx,label='{0} ventanas, estimada'.format(v[k]))
    plt.psd(x,Fs=fs,NFFT=Nv,noverlap=Nv/2,window=chebwin(Nv,at=100))        # psd real segun pyplot para comparar nuestra estimación
    legend = plt.legend(loc='upper right', shadow=True, fontsize='small')
    legend.get_frame().set_facecolor('pink') 
    plt.xlabel('Hz')
    plt.ylabel('dBs')
    plt.grid()

#%%#    ventanas de media movil
    
filt,tw=Mwindow(p,20,.5)
t=np.arange(0,ts*N,ts)
plt.figure()
plt.plot(t,p,label='original')
plt.plot(tw,filt,label='MAverage')
legend = plt.legend(loc='upper right', shadow=True, fontsize='small')
legend.get_frame().set_facecolor('pink') 
plt.ylabel('Amplitud')
plt.xlabel('tiempo (s)')

#%%     quitando tendencia con pseudo inversa

slcion=inv_fitting(t,p)
plt.figure()
plt.plot(t,p)
plt.plot(t,(t*slcion[0,0])+slcion[0,1])


#referencia: OPPENHEIM PROCESAMIENTO SEÑALES DISCRETAS
#            G. Strang, Linear Algebra and Its Applications, 2nd Ed., Orlando, FL, Academic Press, Inc., 1980, pp. 139-142.