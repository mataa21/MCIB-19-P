# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 05:44:13 2019

@author: alan_
"""

### Se implementaran las funciones desarrolladas en la parte 1 en las registros.

import numpy as np
from matplotlib import pyplot as plt
from scipy import signal as sig
from scipy import fftpack as fft
from scipy.signal import chebwin

def PSD(matrix):
    N1=np.size(matrix,1)
    N2=np.size(matrix,0)
    #rxx=sig.convolve2d(matrix,matrix,mode='same')
    rxx=matrix
    rxx=rxx*(chebwin(N1,at=100).reshape(1,N1)*np.ones((N2,1)))
    sxx=np.fft.fft(rxx,axis=1)
    mag_sxx=(np.abs(sxx[:,0:N1//2])*ts)**2
    mag_sxx=10*np.log10(np.mean(mag_sxx,0))
    F=np.linspace(0,fs//2,len(mag_sxx))
    return F,mag_sxx

def Periodograma(p,v):
    N=len(p)
    if N % v!=0:
        Nzeros=v-(N % v)
        x=np.append(p,np.zeros(Nzeros))      # agregar ceros para mejorar la estimación, y para que el reshape se pueda hacer
    else:
        x=p
    Nv=len(x)//v
    matrix=x.reshape(v,Nv)
    F,sxx=PSD(matrix)
    plt.plot(F[0:len(F)//4],sxx[0:len(F)//4],label=lavel[j])
    legend = plt.legend(loc='upper right', shadow=True, fontsize='small')
    legend.get_frame().set_facecolor('pink') 
    plt.xlabel('Hz')
    plt.ylabel('dBs')
    plt.title(senal[i])
    plt.grid()

def RF(den,num,fm):
    w, h = sig.freqz(den,num)
    h[h==0] = 1E-5

    H    = 20*np.log10( np.abs(h) )
    W    = np.angle  (h)
    W    = np.unwrap (W)
    W    = np.degrees(W)
    w    = np.linspace(0,fs//2,H.shape[0] )

    return w, W, H

def pseudo_inv(A):
    U, Sigma, V = np.linalg.svd(A, full_matrices=False,compute_uv=True)
    Sigma_pseu=1/Sigma
    inv=np.matrix(U)*np.diag(Sigma_pseu)*np.matrix(V)
    return inv

def inv_fitting(t,y):
    A=[]
    A.append(t)
    A.append(np.ones(len(t)))
    A=np.array(A).T
    inv=pseudo_inv(A)
    slcion=np.dot(y.reshape(1,len(y)),inv)
    return slcion


def detrend(y,v):
    N=len(y)
    L=N//v
    ind=0
    ydet=[]
    for i in range(N):
        if (i+1) % L == 0:
            t=np.arange(0,L)
            fit=inv_fitting(t,y[ind:i+1])
            ydet.append(y[ind:i+1]-(t*fit[0,0]+fit[0,1]))
            ind=i+1
    ydet = np.ndarray.flatten(np.array(ydet))
    return ydet
            
def barlett_par (v):
    x1=np.arange(0,v//2)
    x2=np.arange(v//2,v)
    triangular=np.append(2*x1/v,2-(2*x2/v))
    return triangular

def MAverage(x,N):
    L=len(x)-N
    y=np.zeros(L)
    for i in range(L):
        y[i]=x[i]-((1/(N+1))*np.sum(x[i:i+N]))
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

def FFT(s):
    F = fft.fft(s)
    F = 20*np.log10(np.abs(F))
    N= F.shape[0]
    W = chebwin(N,at=100)
    Fw = fft.fft(s*W)
    Fw = 20*np.log10(np.abs(Fw))
    f = np.linspace(0,fs//2,N//2)
    plt.plot(f[0:N//8],F[0:N//8],label=lavel[j])
    legend = plt.legend(loc='upper right', shadow=True, fontsize='small')
    legend.get_frame().set_facecolor('pink') 
    plt.title(senal[i])
    plt.grid()   

## cargar registros y definir constantes
    
senal=['oximetro','respirograma','ecg','emg','pulso']
name=["first.txt","second tab.txt","third.txt","four tab.txt","five.txt"]
lavel = ['Primera Etapa','Segunda etapa', 'Tercera etapa','Cuarta etapa','Quinta etapa']
via="/Users/alan_/OneDrive/Documentos/Metodos computacionales/MCIB-19-P/data/"
fs=1000
ts=1/fs
data=[]
for i in range(len(name)):
    data.append(np.loadtxt(via+name[i]))
#%% Estimación de la densidad espectral de potencia de los registros
    
for i in range(len(senal)):
    plt.figure(figsize=(8,5))
    for j in range(len(lavel)):
        Periodograma(data[j][:,i+1],80)
        
#%%   Filtrado de los registros
data2=data.copy()


#Diseño de los Filtros
fc1 = 2*np.array([.3,48.])/fs  #Frecuencia de Corte ya normalizada del Pasa Banda  entre .5 y 58Hz
b1,a1 = sig.butter(3,fc1, btype='bandpass') #Diseño Pasa Banda
fc2 = 2*1.5/fs  #Frecuencia de Corte ya normalizada para 30 Hz
b2,a2 = sig.butter(5,fc2, btype='low') #Diseño pasa bajas

f0 = 60.0  # Frequency to be removed from signal (Hz)
Q = 30.0  # Quality factor
# Design notch filter
b, a = sig.iirnotch(w0=f0*2/fs,Q=Q)

# Frequency response
freq, h = sig.freqz(b, a)
# Plot
fig, ax = plt.subplots(2, 1, figsize=(8, 6))
ax[0].plot((fs/2)*freq/max(freq), 20*np.log10(abs(h)), color='blue')
ax[0].set_title("Frequency Response")
ax[0].set_ylabel("Amplitude (dB)", color='blue')
ax[0].grid()
ax[1].plot((fs/2)*freq/max(freq), np.unwrap(np.angle(h))*180/np.pi, color='green')
ax[1].set_ylabel("Angle (degrees)", color='green')
ax[1].set_xlabel("frecuencia (Hz)")
ax[1].grid()
plt.show()


#Evaluación Respuesta en frecuencia del filtro
w1,W1,H1 = RF(b1,a1,fs)
w2,W2,H2 = RF(b2,a2,fs)
#Graficar la Magnitud Pasa banda
plt.figure(figsize=(10,4))
plt.semilogx(w1[0:500],H1[0:500],'r',linewidth=2)
plt.title("Filtro pasa banda")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Magnitud (dB)")
plt.grid(True)
plt.show()
#Graficar la Fase Pasa Banda
plt.figure(figsize=(10,4))
plt.semilogx(w1,W1,'r',linewidth=2)
plt.title("Filtro pasa banda")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Fase (grados)")
plt.grid(True)
plt.show()
#Graficar la Magnitud Pasa bajas
plt.figure(figsize=(10,4))
plt.semilogx(w2[0:500],H2[0:500],'r',linewidth=2)
plt.title("Filtro pasa bajas")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Magnitud (dB)")
plt.grid(True)
plt.show()
#Graficar la Fase Pasa Bajas
plt.figure(figsize=(10,4))
plt.semilogx(w2,W2,'r',linewidth=2)
plt.title("Filtro pasa bajas")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Fase (grados)")
plt.grid(True)
plt.show()
plt.show

        

#%%  
acnd=[]
for i in range(len (lavel)):
    sign=[]
 #Filtrado de las Señales
    ox = sig.filtfilt(b2,a2,data2[i][:,1])  #Señal oxímetro registro 1
    resp = sig.filtfilt(b2,a2,data2[i][:,2])  #Señal respirograma registro 1
    ecg = sig.filtfilt(b1,a1,data2[i][:,3])   #Señal ECG registro 1
    emg = sig.filtfilt(b1,a1,data2[i][:,4])   #Señal EMG registro 1
    sign.append( sig.filtfilt(b,a,ox))  #Señal oxímetro registro 1
    sign.append( sig.filtfilt(b,a,resp))  #Señal respirograma registro 1
    sign.append( sig.filtfilt(b,a,ecg))   #Señal ECG registro 1
    sign.append( sig.filtfilt(b,a,emg))   #Señal EMG registro 1
    sign.append(sig.filtfilt(b,a,data2[i][:,5]))
    acnd.append(sign)
# eliminando tendencia lineal
v=[6,6,201,30,23]
for i in range(len(lavel)):
    plt.figure()
    for j in range(len(senal)):
        acnd[i][j]=detrend(acnd[i][j][:],v[j])
    
#densidad espectral de potencia despues de acondicionar
        
for i in range(len(senal)):
    plt.figure(figsize=(8,5))
    for j in range(len(lavel)):
        Periodograma(acnd[j][i][:],80)
        
#visualizacion de registro 3 despues de acondicionar
for i in range(len(senal)):
    plt.figure()
    for j in range(1):
        plt.plot(data2[j][63000:85000,0],acnd[j+3][i][63000:85000],label=lavel[j+3])
        legend = plt.legend(loc='lower left', shadow=True, fontsize='small')
        legend.get_frame().set_facecolor('pink') 
        plt.title(senal[i])
        plt.ylabel('Amplitud')
        plt.xlabel('segundos')
        plt.grid()   
        
np.savez_compressed(via+'acondicionada.npz',oxigenacion=acnd[3][0][:],respiracion=acnd[3][1][:],ecg=acnd[3][2][:],emg=acnd[3][3][:],pulso=acnd[3][4][:],sr=fs)
np.savez_compressed(via+'data.npz',x=data)
