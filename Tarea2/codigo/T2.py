# %%
import warnings
warnings.filterwarnings("ignore")
#import pyfits #modulo para leer archivos fits
from astropy.io import fits  # leer archivos fits
import matplotlib.pyplot as plt  # graficar
import numpy as np  # para manejar los datos
import scipy as sp
from scipy import ndimage
import pandas as pd  # tablas y manejo de datos
from astropy.stats import sigma_clip
from scipy.optimize import curve_fit  # fitear curvas
from matplotlib.font_manager import FontProperties  # fuente de los gráficos
import math

G = 4.302e-6  # [km3kg-1s-2] creo
cubo = fits.open("southgal_fixbadc.fits")  # se abre el cubo de datos
data = cubo[0].data  # extracción matriz de datos
header= cubo[0].header  # extracción del header


def values(h,j):
    """
    Funcion para obtener los valores de los ejes l, b, v
    """
    N = h['NAXIS'+str(j)];
    val = np.zeros(N)
    for i in range(0,N):
        val[i] = (i + 1 - float(h['CRPIX' + str(j)]))*float(h['CDELT' + str(j)]) +\
                  float(h['CRVAL'+str(j)])
    return val


# Estos serán los tres arreglos con los valores reales de los tres ejes
# del cubo
velocidad=values(header,1)
longitud=values(header,2)
latitud=values(header,3)
"""
# Índices de una longitud y latitud específica
lon_i = np.where(longitud== 325.0)[0][0]
lat_i = np.where(latitud== -0.25)[0][0]

i_l = -1
i_b = -1
i_k = -1
"""
#print(lon_i,lat_i)
"""
La figura 1 es determinando el ruido "a mano", abajo se hace automatico
=======================================================================
"""
#fig1, ax1 = plt.subplots(3, 1, figsize=(3.5, 9))
#ax1[0].plot(velocidad, data[lat_i][lon_i][:])
"""
ax1[0].set_xlabel('Velocidad')
ax1[0].set_ylabel('Temperatura', fontsize=18)

#T = data[i_b][i_l][:] #i_b = 14 i_l =200
T = data[14][200][:]
ruido = np.where(T<0.5 ) #unidades de K  # otra forma es T>0.5  
ax1[1].plot(velocidad, T, '.',color='r', label='Signal')
ax1[1].plot(velocidad[ruido], T[ruido],'.', color='b', label='Noise')
ax1[1].legend()
ax1[1].set_xlabel('Velocidad', fontsize=18)
ax1[1].set_ylabel('Temperatura', fontsize=18)
ax1[1].set_title('Separando a mano', fontsize=18)

#T = data[i_b][i_l][:]
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
"""
============================================================
"""
"""
T = data[14][200][:]
#r = sigma_clip(T, sigma=3)
r = sigma_clip(T, sigma_lower=3, sigma_upper=3)
rms = np.sqrt(np.mean(r**2))
rmask = r.mask
v_tan = velocidad[rmask][0]
print ('v_tan con mascara =', v_tan)
"""

# Se crea una funcion que para una longitud(l) fija, se recorre latitud(b) y se calcula el rms de las
# velocidades
# Esta misma funcion recorre el cubo de las velocidades asociadas a l y b, hasta que se llega a una
# velocidad que es 5 veces mayor que el rms, esta ultima se guarda un arreglo


def fmin(l, latitud, vs):
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
R0 = 8.5  # kPc
vsol = 220  # km/s
omegasol = vsol/R0

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

# %%
"""
======================P1======================
"""
#curva de rotacion
fig2, ax2 = plt.subplots(2, 1, figsize=(5.5, 5), sharex=True)
fig2.subplots_adjust(hspace=0)
ax2[0].plot(R, vR, 'maroon')
ax2[0].plot(R, vR, 'r.')
ax2[1].plot(R, vR/R+omegasol, 'maroon')
ax2[1].plot(R, vR/R+omegasol, 'r.')
ax2[1].grid()
ax2[0].set_ylim(50, 255)
ax2[0].set_yticks(np.arange(60, 255, 40))
#ax2[1].set_yticks(np.arange(60, 130, 40))
ax2[0].grid()
ax2[1].set_xlabel("$R$ [kpc]", fontsize=10)
ax2[0].set_ylabel("$V_{tan}$ [km/s]", fontsize=10)
ax2[1].set_ylabel(r"$\omega_{tan}$ [rad/s]", fontsize=10)
#ax2.set_title("Velocidad de rotacion en funcion de R", fontsize=10)
fig2.show()
#fig2.tight_layout()    
fig2.savefig("img/curva_rotacion.pdf")
# %%
"""
======================P2======================
"""
def altura_z(l, b_max):
    return R0*np.cos(l*sp.pi/180.)*np.tan(b_max*sp.pi/180.)

l_list = np.arange(0, 385)
fig4, ax4 = plt.subplots(figsize=(5.5, 2.5))
ax4.scatter(R, altura_z(l_list, bvmin), color='blue', marker='.')
ax4.plot(R, altura_z(l_list, bvmin), color='green', linewidth=0.5)
ax4.set_xlabel('R [kpc]')
ax4.set_ylabel('z [kpc]')
ax4.grid()
fig4.show()
fig4.tight_layout()
fig4.savefig("img/corrugacion.pdf")


# %%
"""
======================P3======================
"""
G = 4.302e-6
# Funciones de la velocidad tangencial para cada modelo
def distribucion_masa_1(r, M_0):
    return np.sqrt(G*M_0/r)


def distribucion_masa_2(r, S):
    return np.sqrt(G*(np.pi*r**2*S)/r)


def distribucion_masa_3(r, rho):
    return np.sqrt(G*(4/3*np.pi*r**3*rho)/r)


def distribucion_masa_4(r, S, M_0):
    return np.sqrt(G*(np.pi*r**2*S + M_0)/r)


def distribucion_masa_5(r, rho, M_0):
    return np.sqrt(G*(4/3*np.pi*r**3*rho + M_0)/r)


popt_1, pcov_1 = curve_fit(distribucion_masa_1, R, vR)
popt_2, pcov_2 = curve_fit(distribucion_masa_2, R, vR)
popt_3, pcov_3 = curve_fit(distribucion_masa_3, R, vR)
popt_4, pcov_4 = curve_fit(distribucion_masa_4, R, vR)
popt_5, pcov_5 = curve_fit(distribucion_masa_5, R, vR)

font = FontProperties()
font.set_family('serif')  # fuente para los plots

fig3, ax3 = plt.subplots(5, 1, figsize=(5.2, 9), sharex=True)
fig3.subplots_adjust(hspace=0)
ax3[0].plot(R, distribucion_masa_1(R, *popt_1))
ax3[0].plot(R, vR, 'maroon')
ax3[0].plot(R, vR, 'r.')
ax3[0].set_ylim(50, 255)
ax3[0].set_yticks(np.arange(60, 255, 40))
ax3[0].grid()
ax3[0].text(0.8, 0.15, 'Masa central', fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[0].transAxes,
         bbox=dict(facecolor='white', alpha=1))
"""
ax3[0].text(0.9, 0.15, '$M_0$='+str(format(popt_1[0],'.3E')), fontsize=8.7,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[0].transAxes,
         bbox=dict(facecolor='white', alpha=1))
"""
ax3[1].plot(R, distribucion_masa_2(R, *popt_2))
ax3[1].plot(R, vR, 'maroon')
ax3[1].plot(R, vR, 'r.')
ax3[1].set_ylim(50, 255)
ax3[1].set_yticks(np.arange(60, 255, 40))
ax3[1].grid()
ax3[1].text(0.8, 0.15, 'Disco uniforme', fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[1].transAxes,
         bbox=dict(facecolor='white', alpha=1))
"""
ax3[1].text(0.8, 0.3, 'S='+str(popt_2[0]), fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[1].transAxes,
         bbox=dict(facecolor='white', alpha=1))
"""
ax3[2].plot(R, distribucion_masa_3(R, *popt_3))
ax3[2].plot(R, vR, 'maroon')
ax3[2].plot(R, vR, 'r.')
ax3[2].set_ylim(50, 255)
ax3[2].set_yticks(np.arange(60, 255, 40))
ax3[2].grid()
ax3[2].set_ylabel('$V_{tan}$ [km/s]')
ax3[2].text(0.8, 0.15, 'Esfera uniforme', fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[2].transAxes,
         bbox=dict(facecolor='white', alpha=1))
"""
ax3[2].text(0.8, 0.3, r'$\rho=$'+str(popt_3[0]), fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[2].transAxes,
         bbox=dict(facecolor='white', alpha=1))
"""
ax3[3].plot(R, distribucion_masa_4(R, *popt_4))
ax3[3].plot(R, vR, 'maroon')
ax3[3].plot(R, vR, 'r.')
ax3[3].set_ylim(50, 255)
ax3[3].set_yticks(np.arange(60, 255, 40))
ax3[3].grid()
ax3[3].text(0.76, 0.15, 'Disco + masa central', fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[3].transAxes,
         bbox=dict(facecolor='white', alpha=1))
"""
ax3[3].text(0.76, 0.3, 'S='+str(popt_4[0])+'$M_0=$'+str(popt_4[1]), fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[3].transAxes,
         bbox=dict(facecolor='white', alpha=1))
"""
ax3[4].plot(R, distribucion_masa_5(R, *popt_5))
ax3[4].plot(R, vR, 'maroon')
ax3[4].plot(R, vR, 'r.')
ax3[4].set_ylim(50, 255)
ax3[4].set_yticks(np.arange(60, 255, 40))
ax3[4].grid()
ax3[4].set_xlabel('$R$ [kpc]')
ax3[4].text(0.76, 0.15, 'Esfera + masa central', fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[4].transAxes,
         bbox=dict(facecolor='white', alpha=1))
"""
ax3[4].text(0.76, 0.3, r'$\rho =$'+str(popt_5[0])+'$M_0=$'+str(popt_5[1]), fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[4].transAxes,
         bbox=dict(facecolor='white', alpha=1))
"""
fig3.show()
#fig3.tight_layout()
fig3.savefig("img/fiteo_rotacion.pdf")
# %%