# To attenuate
#------------------------------------------------------------------------------------
#   Emission lines:
#------------------------------------------------------------------------------------

waveHa=6563. #A

#------------------------------------------------------------------------------------
#   Cosmology:
#------------------------------------------------------------------------------------
h0 = 0.6777
omega0=0.307
omegab = 0.0482
lambda0=0.693
cosmo = Cosmology(h=h0,omegaM=omega0, omegaL=lambda0)

#------------------------------------------------------------------------------------
#   Conversion factors and constants:
#------------------------------------------------------------------------------------
kg_to_Msun=1./(1.989e30)
Msun_to_kg=1.989e30
Mpc_to_cm=3.086e24


#------------------------------------------------------------------------------------
#   Inputs to be changed:
#------------------------------------------------------------------------------------

path = '/dust-master/' #working directory
redshift=0.09 #redshift of the SAM catalogue
#------------------------------------------------------------------------------------
#   Read input file and implement attenuation:
#   the input file must include: 
#   - Ha luminosity logscale:        log(L/ergs^-1) 
#   - Cold gas mass of the disc:     Mcold_disk [Msun] 
#   - Half-mass radius of the disc:  rhalf_mass_disc [Mpc], 
#   - Metallicity of the disc:       Z_disc
#
#   !!WARNING!! SAG provides 'rhalf_disc' which is the scale radius r0 with wrong name.
#   In order to get rhalf_mass_disc, apply the correction:
#   rhalf_mass_disc=1.68*r0,  see Gonzalez et al. 2009
#------------------------------------------------------------------------------------

hdulista = pyfits.open(path+'SAGinput_z=0.09.fits')
Mcold_disc=hdulista[1].data.field('Mcold_disc')  #Msun
rhalf_mass_disc=hdulista[1].data.field('rhalf_mass_disc')  #Mpc
logL=hdulista[1].data.field('logL') 
Z_disc=hdulista[1].data.field('Z_disc') 

logLatt, logFatt, logF, attenuation_disc = attenuatedQ(logL, Mcold_disc, rhalf_mass_disc, Z_disc)

#------------------------------------------------------------------------------------
#   Write output file with columns:  
#   - non-attenuated ELG luminosity: logL 
#   - attenuated ELG luminosity:     logLatt 
#   - non-attenuated ELG flux:       logF 
#   - attenuated ELG flux:           logFatt
#   - attenuation coefficient:       attenuation_disc
#------------------------------------------------------------------------------------
out=np.array([logL, logLatt, logF, logFatt, attenuation_disc])
out=np.transpose(out)      
np.savetxt(path+'SAGoutput_z=0.09_dust.txt', out, fmt='%12.9f', header='#log10(L/ergs^-1)    log10(Latt/ergs^-1)  log10(F/ergs^-1cm^-2)     log10(Fatt/ergs^-1cm^-2)    att_coeff') 

