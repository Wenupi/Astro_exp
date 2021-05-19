# Solo para omitir warnings, en general NO deberían usarlo
import warnings
warnings.filterwarnings("ignore")

#import pyfits #modulo para leer archivos fits
from astropy.io import fits
import matplotlib.pyplot as plt #modulo para graficar
import numpy as np #este modulo es para trabajar con matrices como en matlab
import scipy as sp
from scipy import ndimage
import pandas as pd
from astropy.stats import sigma_clip

cubo = fits.open("southgal_fixbadc.fits")  #abrir objeto cubo de datos
data = cubo[0].data  #extraer matriz de datos
header= cubo[0].header  #extraer el header del archivo fits


def values(h,j):  # funcion para obtener los valores de los ejes l, b, v
    N=h['NAXIS'+str(j)];
    val=np.zeros(N)
    for i in range(0,N):
        val[i] = (i+1-float(h['CRPIX'+str(j)]))*float(h['CDELT'+str(j)]) +\
                  float(h['CRVAL'+str(j)])
    return val


# Estos seran los tres arreglos con los valores reales de los tres ejes
# del cubo
velocidad=values(header,1)
longitud=values(header,2)
latitud=values(header,3)

lon_i = np.where(longitud== 325.0)[0][0]
lat_i = np.where(latitud== -0.25)[0][0]

i_l=-1
i_b=-1
i_k=-1

#print(lon_i,lat_i)
"""
La figura 1 es determinando el ruido "a mano", abajo se hace automatico
=======================================================================
"""
fig1, ax1 = plt.subplots(3, 1, figsize=(3.5, 9))
ax1[0].plot(velocidad, data[lat_i][lon_i][:])
ax1[0].set_xlabel('Velocidad')
ax1[0].set_ylabel('Temperatura', fontsize=18)

T = data[i_b][i_l][:] #i_b = 14 i_l =200
ruido = np.where(T<0.5 ) #unidades de K  # otra forma es T>0.5  
ax1[1].plot(velocidad, T, '.',color='r', label='Signal')
ax1[1].plot(velocidad[ruido], T[ruido],'.', color='b', label='Noise')
ax1[1].legend()
ax1[1].set_xlabel('Velocidad', fontsize=18)
ax1[1].set_ylabel('Temperatura', fontsize=18)
ax1[1].set_title('Separando a mano', fontsize=18)

T = data[i_b][i_l][:]
ruido = np.where(T<0.5 ) #otra forma es T>0.5
ax1[2].plot(velocidad, T,color='r', label='Signal')
ax1[2].plot(velocidad[ruido], T[ruido], color='b', label='Noise')
ax1[2].legend()
ax1[2].set_xlabel('Velocidad', fontsize=18)
ax1[2].set_ylabel('Temperatura', fontsize=18)
ax1[2].set_title('Separando a mano', fontsize=18)
fig1.tight_layout()    
fig1.savefig("espectro_basico")
"""
============================================================
"""
#r = sigma_clip(T, sigma=3)
r = sigma_clip(T, sigma_lower=3, sigma_upper=3)
rms = np.sqrt(np.mean(r**2))
rmask = r.mask
v_tan = velocidad[rmask][0]
print ('v_tan con mascara =', v_tan)

# Se crea una funcion que para una longitud(l) fija, se recorre latitud(b) y se calcula el rms de las
# velocidades
# Esta misma funcion recorre el cubo de las velocidades asociadas a l y b, hasta que se llega a una
# velocidad que es 5 veces mayor que el rms, esta ultima se guarda un arreglo


def fmin(l,latitud,vs):
    #recorre latitud
    for q in range(33):
        T1 = data[q][l][:]
        rms = np.sqrt(np.mean(T1**2))   #calcula rms
        #recorre velocidad
        for w in range(306):
            if data[q][l][w]>=5*rms:  #buscamos que no sea ruido
                vs[q]=velocidad[w]  #guardamos la primera v donde T mayor a 5rms
                break


vmin = np.zeros(385)
bvmin = np.zeros(385)
R = np.zeros(385)
R0 = 8.5 # kPc
vsol = 220

# maximorum
# Se recorren las longitudes y se busca la velocidad más negativa (mayor en modulo), se guarda esta
# y su latitud asociada
# Se obtiene un arreglo de R con la ecuacion R =| R0 · cos(l π/180 ) |


for i in range(385):
    vs = np.zeros(33)
    fmin(i,latitud,vs)
    v1 = vs[0]
    b1 = latitud[0]
    for j in range(32):
        if vs[j+1]<v1:
            v1 = vs[j+1]
            b1 = latitud[j+1]
    vmin[i] = v1
    bvmin[i] = b1
    R[i] = np.abs(R0*sp.sin(longitud[i]*sp.pi/180.))  # R0 sin(l)

    
# Se obtiene la Vtan con Vtan = −Vmin − Vsol · sin(lπ/180 ), donde Vmin es la velocidad mayor en
# modulo para l, y Vsol es la velocidad de rotacion del sol.    
#velocidad de rotacion

vR = np.zeros(385)
for i in range(385):
    vR[i] = vmin[i]*(np.abs(sp.sin(longitud[i]*sp.pi/180.))/\
                     sp.sin(longitud[i]*sp.pi/180.)) +\
                     np.abs(vsol*sp.sin(longitud[i]*sp.pi/180.))

#curva de rotacion
fig2, ax2 = plt.subplots(figsize=(3.5, 2.5))
ax2.plot(R,vR, 'maroon')
ax2.plot(R,vR, 'r.')
ax2.grid()
ax2.set_xlabel("R", fontsize=10)
ax2.set_ylabel("Vtan", fontsize=10)
ax2.set_title("Velocidad de rotacion en funcion de R", fontsize=10)
fig2.tight_layout()    
fig2.savefig("primera_curva_rotacion")
