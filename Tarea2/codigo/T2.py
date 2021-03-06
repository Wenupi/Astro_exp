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
from astropy import units as u
from matplotlib.font_manager import FontProperties  # fuente de los gráficos
import math

G = 4.302e-6  # Constante gravitacion de Newton
cubo = fits.open("southgal_fixbadc.fits")  # se abre el cubo de datos
data = cubo[0].data  # extracción matriz de datos
header = cubo[0].header  # extracción del header


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
velocidad = values(header,1)
longitud = values(header,2)
latitud = values(header,3)


# Se crea una funcion que para una longitud(l) fija, se recorre latitud(b) y se calcula el rms de las
# velocidades
# Esta misma funcion recorre el cubo de las velocidades asociadas a l y b, hasta que se llega a una
# velocidad que es 5 veces mayor que el rms, esta ultima se guarda un arreglo


def fmin(l, latitud, vs):
    #recorre latitud
    for q in range(33):
        T1 = data[q][l][:]
        noise = sigma_clip(T1, sigma=5)  # solo se toma la señal mayor a 5 sigma
        rms = np.sqrt(np.mean(noise**2))   #calcula rms del ruido
        #recorre velocidad
        for w in range(306):
            if data[q][l][w]>=5*rms:  #buscamos que no sea ruido
                vs[q]=velocidad[w]  #guardamos la primera v donde T mayor a 5rms
                break


vmin = np.zeros(385)
bvmin = np.zeros(385)
R = np.zeros(385)
#R0 = 8.5*10**3*3.086e+13  # [kPc]->[km]
R0 = 8.5*10**3*u.parsec  # [pc]
vsol = 220*u.kilometer/u.second  # [km/s]
omegasol = vsol/(R0.to(u.kilometer))

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
    R[i] = np.abs((R0.to(u.kilometer)).value*sp.sin(longitud[i]*sp.pi/180.))  # R0 sin(l) [km]

    
# Se obtiene la Vtan con Vtan = −Vmin − Vsol · sin(lπ/180 ), donde Vmin es la velocidad mayor en
# modulo para l, y Vsol es la velocidad de rotacion del sol.    
#velocidad de rotacion

vR = np.zeros(385)
for i in range(385):
    vR[i] = vmin[i]*(np.abs(sp.sin(longitud[i]*sp.pi/180.))/\
                     sp.sin(longitud[i]*sp.pi/180.)) +\
                     np.abs((vsol.value)*sp.sin(longitud[i]*sp.pi/180.))


# %%
"""
======================P1======================
"""
t = 3.2408*10**(-14)*10**(-3)  # kilometro -> kilo parsec

#curva de rotacion
fig2, ax2 = plt.subplots(2, 1, figsize=(5.5, 5), sharex=True)
fig2.subplots_adjust(hspace=0)
ax2[0].plot(R*t, vR, 'maroon')
ax2[0].plot(R*t, vR, 'r.')
ax2[1].plot(R*t, -vmin/R+omegasol.value, 'maroon')
ax2[1].plot(R*t, -vmin/R+omegasol.value, 'r.')
ax2[1].grid()
ax2[0].set_ylim(50, 255)
ax2[0].set_yticks(np.arange(60, 255, 40))
ax2[0].grid()
ax2[1].set_xlabel("$R$ [kpc]", fontsize=10)
ax2[0].set_ylabel("$V_{tan}$ [km/s]", fontsize=10)
ax2[1].set_ylabel(r"$\omega_{tan}$ [rad/s]", fontsize=10)
fig2.show()
#fig2.tight_layout()    
fig2.savefig("img/curva_rotacion.pdf")
# %%
"""
======================P2======================
"""
def altura_z(l, b_max):
    return R0.value*np.cos(l*sp.pi/180.)*np.tan(b_max*sp.pi/180.)  #*1e-3*3.24078e-14

l_list = np.arange(0, 385)
fig4, ax4 = plt.subplots(figsize=(5.3, 2.5))
ax4.scatter(R*t, altura_z(longitud, bvmin)*10**(-3), color='blue', marker='.')
ax4.plot(R*t, altura_z(longitud, bvmin)*10**(-3), color='green', linewidth=0.5)
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

t_2 = 3.086*10**(13)  # parsec -> kilometro
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


popt_1, pcov_1 = curve_fit(distribucion_masa_1, R*t, vR, bounds=(0, [np.inf]))
popt_2, pcov_2 = curve_fit(distribucion_masa_2, R*t, vR, bounds=(0, [np.inf]))
popt_3, pcov_3 = curve_fit(distribucion_masa_3, R*t, vR, bounds=(0, [np.inf]))
popt_4, pcov_4 = curve_fit(distribucion_masa_4, R*t, vR, bounds=(0, [np.inf, np.inf]))
popt_5, pcov_5 = curve_fit(distribucion_masa_5, R*t, vR, bounds=(0, [np.inf, np.inf]))

font = FontProperties()
font.set_family('serif')  # fuente para los plots

fig3, ax3 = plt.subplots(5, 1, figsize=(5.2, 9), sharex=True)
fig3.subplots_adjust(hspace=0)
ax3[0].plot(R*t, distribucion_masa_1(R*t, *popt_1))
ax3[0].plot(R*t, vR, 'maroon')
ax3[0].plot(R*t, vR, 'r.')
ax3[0].set_ylim(50, 255)
ax3[0].set_yticks(np.arange(60, 255, 40))
ax3[0].grid()
ax3[0].text(0.8, 0.15, 'Masa central', fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[0].transAxes,
         bbox=dict(facecolor='white', alpha=1))


ax3[1].plot(R*t, distribucion_masa_2(R*t, *popt_2))
ax3[1].plot(R*t, vR, 'maroon')
ax3[1].plot(R*t, vR, 'r.')
ax3[1].set_ylim(50, 255)
ax3[1].set_yticks(np.arange(60, 255, 40))
ax3[1].grid()
ax3[1].text(0.8, 0.15, 'Disco uniforme', fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[1].transAxes,
         bbox=dict(facecolor='white', alpha=1))


ax3[2].plot(R*t, distribucion_masa_3(R*t, *popt_3))
ax3[2].plot(R*t, vR, 'maroon')
ax3[2].plot(R*t, vR, 'r.')
ax3[2].set_ylim(50, 255)
ax3[2].set_yticks(np.arange(60, 255, 40))
ax3[2].grid()
ax3[2].set_ylabel('$V_{tan}$ [km/s]')
ax3[2].text(0.8, 0.15, 'Esfera uniforme', fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[2].transAxes,
         bbox=dict(facecolor='white', alpha=1))


ax3[3].plot(R*t, distribucion_masa_4(R*t, *popt_4))
ax3[3].plot(R*t, vR, 'maroon')
ax3[3].plot(R*t, vR, 'r.')
ax3[3].set_ylim(50, 255)
ax3[3].set_yticks(np.arange(60, 255, 40))
ax3[3].grid()
ax3[3].text(0.76, 0.15, 'Disco + masa central', fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[3].transAxes,
         bbox=dict(facecolor='white', alpha=1))


ax3[4].plot(R*t, distribucion_masa_5(R*t, *popt_5))
ax3[4].plot(R*t, vR, 'maroon')
ax3[4].plot(R*t, vR, 'r.')
ax3[4].set_ylim(50, 255)
ax3[4].set_yticks(np.arange(60, 255, 40))
ax3[4].grid()
ax3[4].set_xlabel('$R$ [kpc]')
ax3[4].text(0.76, 0.15, 'Esfera + masa central', fontsize=11,
         fontproperties=font, horizontalalignment='center',
         verticalalignment='center', transform=ax3[4].transAxes,
         bbox=dict(facecolor='white', alpha=1))

fig3.show()
#fig3.tight_layout()
fig3.savefig("img/fiteo_rotacion.pdf")

error_1 = np.mean((distribucion_masa_1(R*t, *popt_1)-vR)**2)
error_2 = np.mean((distribucion_masa_2(R*t, *popt_2)-vR)**2)
error_3 = np.mean((distribucion_masa_3(R*t, *popt_3)-vR)**2)
error_4 = np.mean((distribucion_masa_4(R*t, *popt_4)-vR)**2)
error_5 = np.mean((distribucion_masa_5(R*t, *popt_5)-vR)**2)

# %%