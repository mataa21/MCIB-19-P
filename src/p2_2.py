# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 17:48:52 2019

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
    mag_sxx=(np.abs(sxx[:,0:N1/2])*ts)**2
    mag_sxx=10*np.log10(np.mean(mag_sxx,0))
    F=np.linspace(0,fs/2,len(mag_sxx))
    return F,mag_sxx

def Periodograma(p,v):
    N=len(p)
    if N % v!=0:
        Nzeros=v-(N % v)
        x=np.append(p,np.zeros(Nzeros))      # agregar ceros para mejorar la estimación, y para que el reshape se pueda hacer
    else:
        x=p
    Nv=len(x)/v
    matrix=x.reshape(v,Nv)
    F,sxx=PSD(matrix)
    plt.plot(F[0:len(F)/4],sxx[0:len(F)/4],label=lavel[j])
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
    w    = np.linspace(0,fs/2,H.shape[0] )

    return w, W, H

def pseudo_inv(A):
    U, Sigma, V = np.linalg.svd(A, full_matrices=False,compute_uv=True)
    Sigma_pseu=1./Sigma
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
    L=N/v
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

def FFT(s):
    F = fft.fft(s)
    F = 20*np.log10(np.abs(F))
    N= F.shape[0]
    W = chebwin(N,at=100)
    Fw = fft.fft(s*W)
    Fw = 20*np.log10(np.abs(Fw))
    f = np.linspace(0,fs/2,N/2)
    plt.plot(f[0:N/8],F[0:N/8],label=lavel[j])
    legend = plt.legend(loc='upper right', shadow=True, fontsize='small')
    legend.get_frame().set_facecolor('pink') 
    plt.title(senal[i])
    plt.grid()   

senal=['oximetro','respirograma','ecg','emg','pulso']
name=["first.txt","second tab.txt","third.txt","four tab.txt","five.txt"]
lavel = ['Primera Etapa','Segunda etapa', 'Tercera etapa','Cuarta etapa','Quinta etapa']
via="/Users/alan_/OneDrive/Documentos/Metodos computacionales/MCIB-19-P/data/"
fs=1000
ts=1./fs
data=[]
for i in range(len(name)):
    data.append(np.loadtxt(via+name[i]))
#%% Estimación de la densidad espectral de potencia de los registros
    
for i in range(len(senal)):
    plt.figure(figsize=(8,5))
    for j in range(len(lavel)):
        Periodograma(data[j][:,i+1],80)
        
#%%   Filtrado de los registros
        
 #Declaración de las señales a filtrar
oxi_1 = data[1][:,1]   #Oxímetro primer registro
oxi_2 = data[3][:,1]   #Oxímetro segundo registro
respi_1 = data[1][:,2]   #Resporigrama primer registro
respi_2 = data[3][:,2]   #Respirograma segundo registro
ecg_1 = data[1][:,3]   #ECG primer registro
ecg_2 = data[3][:,3]   #ECG segundo registro
emg_1 = data[1][:,4]   #EMG primer registro
emg_2 = data[3][:,4]   #EMG segundo registro
#Diseño de los Filtros
fc1 = 2*np.array([.5,58.])/fs  #Frecuencia de Corte ya normalizada del Pasa Banda  entre .5 y 58Hz
b1,a1 = sig.butter(3,fc1, btype='band') #Diseño Pasa Banda
fc2 = 2.*30/fs  #Frecuencia de Corte ya normalizada para 30 Hz
b2,a2 = sig.butter(3,fc2, btype='low') #Diseño pasa bajas
 #Filtrado de las Señales
oxi_1_fil = sig.filtfilt(b2,a2,oxi_1)  #Señal oxímetro registro 1
oxi_2_fil = sig.filtfilt(b2,a2,oxi_2)   #Señal oxímetro registro 2
respi_1_fil = sig.filtfilt(b2,a2,respi_1)  #Señal respirograma registro 1
respi_2_fil = sig.filtfilt(b2,a2,respi_2)  #Señal respirograma registro 2
ecg_1_fil = sig.filtfilt(b1,a1,ecg_1)   #Señal ECG registro 1
ecg_2_fil = sig.filtfilt(b1,a1,ecg_2)  #Señal ECG registro 2
emg_1_fil = sig.filtfilt(b1,a1,emg_1)   #Señal EMG registro 1
emg_2_fil = sig.filtfilt(b1,a1,emg_2)  #Señal EMG registro 2
#Evaluación Respuesta en frecuencia del filtro
w1,W1,H1 = RF(b1,a1,fs)
w2,W2,H2 = RF(b2,a2,fs)
#Graficar la Magnitud Pasa banda
plt.figure(figsize=(10,4))
plt.plot(w1[0:500],H1[0:500],'r',linewidth=2)
plt.title("Filtro pasa banda")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Magnitud (dB)")
plt.grid(True)
plt.show()
#Graficar la Fase Pasa Banda
plt.figure(figsize=(10,4))
plt.plot(w1,W1,'r',linewidth=2)
plt.title("Filtro pasa banda")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Fase (grados)")
plt.grid(True)
plt.show()
#Graficar la Magnitud Pasa bajas
plt.figure(figsize=(10,4))
plt.plot(w2[0:500],H2[0:500],'r',linewidth=2)
plt.title("Filtro pasa bajas")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Magnitud (dB)")
plt.grid(True)
plt.show()
#Graficar la Fase Pasa Bajas
plt.figure(figsize=(10,4))
plt.plot(w2,W2,'r',linewidth=2)
plt.title("Filtro pasa bajas")
plt.xlabel("Frecuencia (Hz)")
plt.ylabel("Fase (grados)")
plt.grid(True)
plt.show()
 #Señales Filtradas
#Oxímetro
plt.figure(figsize=(10,4))
plt.plot(oxi_1_fil[2000:4000], label = ' Oximetro primer')
plt.legend(loc = 'upper right')
plt.plot(oxi_2_fil[2000:4000], label = ' Oximetro tercer ')
plt.legend(loc = 'upper right')
plt.title('Oximetro 1 y 3 Filtradas')

plt.show
#Respirograma
plt.figure(figsize=(10,4))
plt.plot(respi_1_fil[2000:4000], label = 'Senal Respirograma primer registro')
plt.legend(loc = 'upper right')
plt.plot(respi_2_fil[2000:4000], label = 'Senal Respirograma tercer registro')
plt.legend(loc = 'upper right')
plt.title('Senales de Respirograma del Registro 1 y 3 Filtradas')

plt.show
#ECG
plt.figure(figsize=(10,4))
plt.plot(ecg_1_fil[2000:4000], label = 'Senal ECG primer registro')
plt.legend(loc = 'upper right')
plt.plot(ecg_2_fil[2000:4000], label = 'Senal ECG tercer registro')
plt.legend(loc = 'upper right')
plt.title('Senales de ECG del Registro 1 y 3 Filtradas')

plt.show
#Oxímetro
plt.figure(figsize=(10,4))
plt.plot(emg_1_fil[2000:4000], label = 'Senal EMG primer registro')
plt.legend(loc = 'upper right')
plt.plot(emg_2_fil[2000:4000], label = 'Senal EMG tercer registro')
plt.legend(loc = 'upper right')
plt.title('Senales de EMG del Registro 1 y 3 Filtradas')

plt.show

#%%        
# Grafica de los registros temporales para percibir tendencias

for i in range(len(senal)):
    plt.figure()
    for j in range(len(lavel)):
        plt.plot(data[j][0:4000,0],data[j][0:4000,i+1],label=lavel[j])
        legend = plt.legend(loc='lower left', shadow=True, fontsize='small')
        legend.get_frame().set_facecolor('pink') 
        plt.title(senal[i])
        plt.ylabel('Amplitud')
        plt.xlabel('segundos')
        plt.grid()   
        

#%%            
# eliminando tendencia lineal
        
v=[3,10,3,3,8]
for i in range(len(senal)):
    plt.figure()
    for j in range(len(lavel)):
        ydetrend=detrend(data[j][0:4000,i+1],v[i])
        plt.plot(np.arange(0,len(ydetrend)*ts,ts),ydetrend,label=lavel[j])
        legend = plt.legend(loc='lower left', shadow=True, fontsize='small')
        legend.get_frame().set_facecolor('pink') 
        plt.title(senal[i])
        plt.ylabel('Amplitud')
        plt.xlabel('segundos')
        plt.grid()   
        