import read_jc_obs as jc
from matplotlib import pyplot as plt  # para usar legend
from pylab import figure, plot, xlabel, ylabel, errorbar, legend, xlim
from numpy import zeros, size, histogram, log10, append, arange

from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table

obs_path = '/Users/pablofernandez/Documents/5º curso/TFG2/Datos/Comparat_2016/'
redshift = 0.94

ox, oy, el, eh = jc.read_jc_lf(obs_path,redshift,
                               infile='O2_3728-data-summary-Planck15.txt')

# 'ox' es la luminosidad y 'oy' la funcion luminosidad (ambos en escala logaritmica)
# 'el' es el error inferior y 'eh' es el error superior de los datos de 'oy'

# Creo una matriz que guarde estos errores para poder usar la funcion errorbar

Merror = zeros((2,size(el)))
Merror[0,:]=el                 # en la primera fila se pone el inferior
Merror[1,:]=eh                 # en la segunda fila se pone el superior

# Dibujamos la funion luminosidad con las barras de error

figure(1)
errorbar(ox, oy, Merror, fmt='o', color='grey', capsize=2)
# capsize pone remates a las barras de error

xlabel("$log_{10}$(L[OII][erg $s^{-1}$ $h^{-2}$]])")
ylabel("$log_{10}$($\Phi$[$Mpc^{-3}$ $dex^{-1}$ $h^{-3}$])")

###############################################################################
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
types = ['LOII', 'LOII_att']
SAM_names =['SAGE','SAG','Galacticus']
colors=['r','b','g']
lines = ['-','--']

figure(2)

for isam, SAM in enumerate(SAM_names):
    infile = SAM +'.fits'
    Data = get_pkg_data_filename(infile)
    datos = Table.read(Data, hdu=1)
    v_datos = zeros((size(types),size(datos)))
    print(size(types),size(datos))
    
    i = -1
    for itype in types:
        i += 1
        Data = get_pkg_data_filename(infile)
        datos = Table.read(Data, hdu=1)
        v_datos[i,:] = datos[itype]+2*log10(hu)
        
    for i in range(size(types)):
        H, bins_edges = histogram(v_datos[i],bins=vbins)
        y = log10(H/(vol*dl))
         
        plot(lhist,y,color=colors[isam],linestyle=lines[i],
             label=SAM+'('+types[i]+')')    
    
leg = plt.legend(loc=0)

errorbar(ox, oy, Merror, fmt='o', color='grey', capsize=2)

xlabel("log$_{10}$(L[OII][erg s$^{-1}$ h$^{-2}$]])",size=15)
ylabel("log$_{10}$($\Phi$[Mpc$^{-3}$ dex$^{-1}$ h$^{-3}$])",size=15)
xlim(40.7,43.5)


