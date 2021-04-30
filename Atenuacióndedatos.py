import read_jc_obs as jc
import Cosmology as cosmo

from matplotlib import pyplot as plt  # para usar legend
from pylab import figure, plot, xlabel, ylabel, errorbar, xlim, title
from numpy import zeros, size, histogram, log10, pi, sqrt, exp, append, arange

from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table

from Att_Cardelli import att_Car

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

def att_Cal_2(logL, attenuation_disc):
    
    # calcula el log10 de la L atenuada
    logL_attenuated = logL - 0.4*attenuation_disc
    
    return logL_attenuated

###############################################################################

############# DATOS EXPERIMENTALES DEL OII ###############

obs_path = '/Users/pablofernandez/Documents/5º curso/TFG2/Datos/Comparat_2016/'
redshift = 0.94

ox_OII, oy_OII, el_OII, eh_OII = jc.read_jc_lf(obs_path,redshift,
                               infile='O2_3728-data-summary-Planck15.txt')
# 'ox': luminosidad y 'oy': funcion luminosidad (ambos en escala logaritmica)

Merror_OII = zeros((2,size(el_OII)))   # matriz que guarda los errores
Merror_OII[0,:]= el_OII                 # en la primera fila se pone el inferior
Merror_OII[1,:]= eh_OII                 # en la segunda fila se pone el superior

############# DATOS EXPERIMENTALES DE Hbeta ###############

ox_Hb, oy_Hb, el_Hb, eh_Hb = jc.read_jc_lf(obs_path,redshift,
                               infile='H1_4862-data-summary-Planck15.txt')
# 'ox': luminosidad y 'oy': funcion luminosidad (ambos en escala logaritmica)

Merror_Hb = zeros((2,size(el_Hb)))   # matriz que guarda los errores
Merror_Hb[0,:]= el_Hb                 # en la primera fila se pone el inferior
Merror_Hb[1,:]= eh_Hb                # en la segunda fila se pone el superior

# PLOT DE LOS TRES SAMs #

hu = 0.6777 # parámetro Hubble
lbox = 1000 # lado de la caja en Mpc/h
vol = lbox**3 # volumen en Mpc^3/h^3

# intervalo de luminosidad
lmin = 38.
lmax = 44.
dl = 0.1
lbins = arange(lmin,lmax,dl)
lhist = lbins
vbins = append(lbins,lmax) 

# array con los nombres de las columnas de datos
types = ['unatt', 'att']
SAM_names =['SAGE','SAG','Galacticus']
colors=['r','b','g']
lines = ['-','--']

# parámetros para usar las funciones de atenuación de Cardelli y Calzetti
Mcold_disc = 0.2        # (en masas solares)
rhalf_mass_disc = 100   # (en Mpc)
Z_disc = 0.9

#######################

Lclass = 'LOII'    # Tipo de luminosidad para plot: LOII, LOIII, LHa o LHb

if Lclass == 'LOII':
    waveA = 3728        # Longitud de onda correspondiente (en Angstroms)
elif Lclass == 'LOIII':
    waveA = 5007
elif Lclass == 'LHa':
    waveA = 6563
elif Lclass == 'LHb':
    waveA = 4861
            
######################

figure(1)   # plot sin atenuar y atenuado con Calzetti

for isam, SAM in enumerate(SAM_names):
    infile = SAM +'.fits'
    Data = get_pkg_data_filename(infile)
    datos = Table.read(Data, hdu=1)
  
    v_datos = zeros((size(types),size(datos)))
    
    #v_datos = zeros((size(types),size(datos)))    
    print(size(types),size(datos))
    
    # la fila 1 guarda los datos sin atenuar
    
    v_datos[0,:] = datos[Lclass]+2*log10(hu)
    
    # la fila 2 mete atenuación por Calzetti
    
    attenuation_disc = att_Cal_1(waveA, Mcold_disc, rhalf_mass_disc, Z_disc)
    
    for j in range(size(datos)):   # que tome como input vectores
    
        v_datos[1,j] = att_Cal_2(v_datos[0,j], attenuation_disc)
        
    for i in range(size(types)):
        H, bins_edges = histogram(v_datos[i],bins=vbins)
        y = log10(H/(vol*dl))
         
        plot(lhist,y,color=colors[isam],linestyle=lines[i],
             label=SAM+'('+types[i]+')')    
    
leg = plt.legend(loc=0)

# PLOT DATOS EXPERIMENTALES

title('z = %.2f' %redshift)   # el título es el redshift que has puesto

if Lclass == 'LOII':
    errorbar(ox_OII, oy_OII, Merror_OII, fmt='o', color='grey', capsize=2)
    xlabel("$log_{10}$(L[OII][erg $s^{-1}$ $h^{-2}$]])",size=15)
elif Lclass == 'LHb':
    errorbar(ox_Hb, oy_Hb, Merror_Hb, fmt='o', color='grey', capsize=2)
    xlabel("$log_{10}$(L[H $\u03B2$ ][erg $s^{-1}$ $h^{-2}$]])",size=15)

ylabel("$log_{10}$($\Phi$[$Mpc^{-3}$ $dex^{-1}$ $h^{-3}$])",size=15)