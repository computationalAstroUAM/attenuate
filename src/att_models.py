#----------------------------------------------------------------------------
#   Calzetti's law
#   At the moment this function is for removing attenuation from obs
#----------------------------------------------------------------------------
def calzetti_mag(omag,esbv,efflAA):
    #Calzettis's law relates the intrinsic flux, Fi, 
    #with the attenuated one Fo through:
    # Fi(lam)=f0(lam)10**(0.4*es*k)
    # Observed magnitude - Intrinsic magnitude = -2.5*log(Fo/Fi)=
    #      =2.5*log(Fi/Fo)=esbv*k
    # esbv = Colour excess of the stellar continuum (expected input)
    # ebv = color excess derived from the nebular gas emission 
    # esbv = 0.44*ebv
    # effl expected in AA

    Rv = 4.05  # Eq. 5 from Calzetti+2000

    effl = efflAA*1e-4 #Convert to microns

    if (effl >= 0.12 and effl < 0.63):
        k = 2.659*(-2.156 + 1.509/effl - 0.198/(effl**2) +\
                        0.011/(effl**3)) + Rv
    elif (effl >= 0.63 and effl <= 2.2):
        k = 2.659*(-1.857 + 1.040/effl) + Rv
    else:
        print('STOP calzetti_mag: {} out of allowed wavelength range'.format(effl))
        print('STOP calzetti_mag: wavelength (AA)')
        sys.exit()

    # Intrinsic magnitude
    imag = omag - esbv*k
    return imag

#----------------------------------------------------------------------------
#   Cardelli et al. 1989 extinction laws in FIR and IR/OPT:
#   From Favole et al. 2020
#----------------------------------------------------------------------------
def cardelliIR(waveA):
    Rv=3.1 #Mathis et al. 1989  
    wave=waveA/10000.
    x=1./wave
    ax=0.574*(x**1.61) 
    bx=-0.527*(x**1.61)
    return ax+bx/Rv

def cardelliOPT(waveA):
    Rv=3.1 
    wave=waveA/10000.
    x=1./wave
    y=x-1.82
    ax=1.+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4+0.01979*y**5-0.77530*y**6+0.32999*y**7 
    bx=1.411338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5+5.30260*y**6-2.09002*y**7
    return ax+bx/Rv


#-------------------------------------------------------------------------------------
#  Attenuated Luminosity and Flux (log scale)
#-------------------------------------------------------------------------------------
def attenuatedQ(logL, Mcold_disc, rhalf_mass_disc, Z_disc):
    Zsun=0.0134 #Asplund et al. 2009
    s=1.6 #if lambda>2000A
    a_disc=1.68
    mp=1.67e-27 #kg
    Al_Av=cardelliOPT(waveHa)
    costheta=0.30 #scattering forward oriented
    sectheta=1./costheta
    albedo=0.56

    mean_col_dens_disc_log=np.log10(Mcold_disc)+np.log10(Msun_to_kg)-np.log10(1.4*mp*np.pi)-2.*np.log10(a_disc*rhalf_mass_disc*Mpc_to_cm)
    tau_disc=np.log10(Al_Av)+np.log10((Z_disc/Zsun)**s)+mean_col_dens_disc_log-np.log10(2.1e+21)
    tau_disc=10.**tau_disc
    al_disc=(np.sqrt(1.-albedo))*tau_disc
    attenuation_disc=-2.5*np.log10((1.-np.exp(-al_disc*sectheta))/(al_disc*sectheta))
    logL_attenuated=logL-0.4*attenuation_disc
    logflux_att=logL_attenuated-np.log10(4.*np.pi)-2.*np.log10(cosmo.Dl(redshift)*Mpc_to_cm) 
    logflux=logL-np.log10(4.*np.pi)-2.*np.log10(cosmo.Dl(redshift)*Mpc_to_cm) 
    return logL_attenuated, logflux_att, logflux, attenuation_disc       


