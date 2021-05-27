import Cosmology as cosmo
import read_jc_obs as jc

from matplotlib import pyplot as plt  # para usar legend
from pylab import figure, plot, errorbar, xlabel, ylabel, xlim, title
from numpy import zeros, size, histogram, log10, pi, sqrt, exp, append, arange, loadtxt

########################## FUNCIÓN DE CALZETTI ################################

Av = 1
Rv = 4.05
Zsun = 0.0134       # Metalicidad solar, Asplund et al. 2009
s = 1.6             # Exponente s, si lambda>2000A
a_disc = 1.68       # En la fórmula de <NH>
mp = 1.67e-27       # Masa del protón en kg
costheta = 0.30          # scattering forward oriented
sectheta = 1./costheta
albedo = 0.56
Mpc_to_cm = 3.086e+24    # pasar de Mpc a cm
    #se usa porque rhalf_mass_dis viene dado en Mpc
Msun_to_kg = 1.989e+30   # masa del Sol en kg, 
    #se usa porque Mcold_disc viene dado en unidades de masas solares

cosmo.set_cosmology(0.307, 0.048, 0.693, 0.823, 0.96)   # definir cosmología de Planck


def calzettiIR(waveA):  # IR: 6300A < λ < 22000A
    wave = waveA/10000.     # la pasa a micrometros
    x = 1./wave             # x es el inverso
    k = 2.659*(-1.857 + 1.040*x) + Rv
    return k*Av/Rv

def calzettiOPT(waveA):  # Optical/NIR: 1200A < λ < 6300A
    wave = waveA/10000.   # la pasa a micrometros
    x = 1./wave           
    k = 2.659*(-2.156 + 1.509*x - 0.198*x**2 + 0.011*x**3) + Rv
    return k*Av/Rv

def att_Cal_1(waveA, Mcold_disc, rhalf_mass_disc, Z_disc):
    
    if waveA < 6300:            # Dependiendo de λ, calcula Alambda/A_V
        Al_Av = calzettiOPT(waveA)
    else:
        Al_Av = calzettiIR(waveA)
        
    
    mean_col_dens_disc_log = log10(Mcold_disc)+log10(Msun_to_kg)-log10(1.4*mp*pi) \
        -2.*log10(a_disc*rhalf_mass_disc*Mpc_to_cm)
    
    # calcula el log10 de tau_lambda
    tau_disc = log10(Al_Av) + s*log10(Z_disc/Zsun) + mean_col_dens_disc_log - log10(2.1e+21)
    tau_disc = 10.**tau_disc    # deshacer el logaritmo
    
    al_disc = (sqrt(1.-albedo))*tau_disc     # calcula el a_lambda
    
    # calcula la atenuación
    attenuation_disc = -2.5*log10((1.-exp(-al_disc*sectheta))/(al_disc*sectheta))
    
    return attenuation_disc

def att_Cal_2(logL, attenuation_disc, line_factor):
    
    # calcula el log10 de la L atenuada
    logL_attenuated = logL - 0.4*attenuation_disc*(line_factor)
    
    return logL_attenuated

###############################################################################

# Se dibujan 4 plots distintos, cada uno con un z diferente

z = [0.99, 0.69, 0.51, 0.17]

# PLOT HISTOGRAMA
lbox = 500     # lado de la caja en Mpc/h
vol = lbox**3  # volumen en Mpc^3/h^3
    
lmin = 39.
lmax = 43.
dl = 0.1
lbins = arange(lmin,lmax,dl)
lhist = lbins
vbins = append(lbins,lmax) 
    
types = ['unatt','att_continuo', 'att_Saito']
lines = ['-','-.','--']

fig = 0 # parámetro de control de las 4 figuras

for redshift in z:
    
    # DATOS EXPERIMENTALES DEL OII
    obs_path = '/Users/pablofernandez/Documents/5º curso/TFG2/Datos/Comparat_2016/'

    ox, oy, el, eh = jc.read_jc_lf(obs_path,redshift,
                               infile='O2_3728-data-summary-Planck15.txt')

    Merror = zeros((2,size(el)))  # guarda las barras de error
    Merror[0,:]=el
    Merror[1,:]=eh
    
    # DATOS MILLENIUM: 4 ARCHIVOS DIFERENTES
    if redshift == 0.99:
        data = loadtxt('LOII_41.txt', float, delimiter=',')
    if redshift == 0.69:
        data = loadtxt('LOII_45.txt', float, delimiter=',')
    if redshift == 0.51:
        data = loadtxt('LOII_48.txt', float, delimiter=',')
    if redshift == 0.17:
        data = loadtxt('LOII_56.txt', float, delimiter=',')

    L_OII_noext = log10(data[:, 0])+40   # para que esté en log10(L[erg/s h^-2])
    L_OII_ext = log10(data[:, 1])+40
    Mcold_disc = data[:, 2]
    metals_coldgas = data[:, 3]
    rdisk = data[:, 4]

    Z_disc = metals_coldgas/Mcold_disc
    Mcold_disc = Mcold_disc*1e10       # (en masas solares/h)
    rhalf_mass_disc = rdisk            # (en Mpc/h)
    
    # PLOT HISTOGRAMA
    
    fig = fig + 1
    figure(fig)

    v_datos = zeros((size(types),size(L_OII_noext)))
    # sin atenuar
    v_datos[0,:] = L_OII_noext
        
        
    for j in range(size(L_OII_noext)):
        
        attenuation_disc = att_Cal_1(3737, Mcold_disc[j], rhalf_mass_disc[j], Z_disc[j])
        
        # atenuación en el continuo
        v_datos[1,j] = att_Cal_2(v_datos[0,j], attenuation_disc, 1)
        # atenuación con factor de Saito
        v_datos[2,j] = att_Cal_2(v_datos[0,j], attenuation_disc, 5/(redshift+2.2))
            
    for i in range(size(types)):
        H, bins_edges = histogram(v_datos[i],bins=vbins)
        y = log10(H/(vol*dl))
             
        plot(lhist,y,linestyle=lines[i],label='('+types[i]+')')    
        
    leg = plt.legend(loc=0)
    
    # PLOT DATOS EXPERIMENTALES
    
    errorbar(ox, oy, Merror, fmt='o', color='grey', capsize=2)
    
    title('z = %.2f' %redshift) 
    xlabel("log$_{10}$(L[OII][erg s$^{-1}$ h$^{-2}$]])",size=15)
    ylabel("log$_{10}$($\Phi$[Mpc$^{-3}$ dex$^{-1}$ h$^{-3}$])",size=15)