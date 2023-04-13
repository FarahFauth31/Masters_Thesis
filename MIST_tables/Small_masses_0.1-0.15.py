#%% First try

# import numpy as np
# import astropy.units as u
# import astropy.constants as const
# import sys
# import pickle

# sys.path.append( '/home/farah/Documents/Project/Data/MIST_tables/' )

# #Create mass array to find nearest neighbour
# MASSES = np.arange(0.1,2.05, 0.05) #all the MIST masses we have

# #Find nearest neighbour in array
# def find_nearest(array, value):
#     array = np.asarray(array)
#     idx = (np.abs(array - value)).argmin()
#     return array[idx]

# def load_mist_tables(Mstar=1., filepath='/home/farah/Documents/Project/Data/MIST_tables/'):
#         """
#         Load in the MIST tables.
#         Mstar: Stellar masses in units of solar masses
#         filepath: Path where the MIST tables are stored
#         """
#         import read_mist_models
#         import astropy.units as u

#         print(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)))

#         eep = read_mist_models.EEP(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)), verbose=False)
#         AGE_mist = eep.eeps['star_age']*u.yr # stellar age in years
#         TAU_mist = (eep.eeps['conv_env_turnover_time_g']*u.s).to(u.d) # convective turnover time in days
#         MOI_mist = eep.eeps['moment_of_inertia']*u.g*u.cm**2. # moment of inertia in cgs

#         return AGE_mist, TAU_mist, MOI_mist


# def define_cgs_constants():
#     import astropy.constants as const
    
#     M_sun = const.M_sun.cgs     # solar mass [gr]
#     R_sun = const.R_sun.cgs     # solar radius [cm]
#     G = const.G.cgs             # Newton's constant in cm^3 g^-1 s^-2

#     return M_sun, R_sun, G


# def spin_down_evol(Prot_init, Mstar_init, t0=1., tdisc=13., a=0.02, b=2., J0=1e41, n_min=1., complexity=True):

#     """
#     Evolves the stellar rotation period based on the spin-evolution model published in Garraffo et al. (2015, 2016, 2018)
#     Prot_init: Starting rotation periods
#     Mstar_init: Corresponding stellar masses
#     t0: Time in Myr at which the rotational evolution should be started
#     a, b, J0, n_min: Prefactors for the n-vs-Ro relation in Garraffo+2018 
#     """
#     import astropy.units as u
#     from scipy.interpolate import InterpolatedUnivariateSpline
#     from scipy.interpolate import UnivariateSpline
#     #from scipy.signal import savgol_filter
    
#     M_sun, R_Sun, G = define_cgs_constants()

#     # DEFINE PARAMETERS
#     Prot0 = Prot_init*u.d               # initial rotation periods
#     Mstar = Mstar_init*u.Msun           # initial stellar masses
#     N = len(Prot_init)                  # sample size (number of stars)

#     age_max = 1e10*u.yr                 # oldest cluster's age [million years]
#     tf = age_max.to(u.d)                # final time (cluster's age)  [d]
#     ti = (t0*u.Myr).to(u.d)             # initial age [d]  --> age of ONC

#     # define empty arrays to which the data will be saved
#     Prot_evol = []
#     Prot_interp = []

#     tdisc *= u.Myr
#     t0 *= u.Myr

#     ### EVOLVE ALGORITHM
#     for j in range(N): 
#         ### LOAD MIST table
#         AGE_mist, TAU_mist, MOI_mist = load_mist_tables(Mstar=Mstar[j].value) 
                
#         # time at which each star will start its rotational evolution [Myrs]
#         t_DL = (tdisc-t0).to(u.yr)

#         # find index in MIST table at which t = t_DL
#         ind0 = min(range(len(AGE_mist)), key= lambda x:abs((AGE_mist[x]-t_DL).value))

#         # Rotation Periods [d]
#         Prot = np.tile(Prot0[j], len(AGE_mist))

#         # Angular Momentum [g cm^2/d]
#         J = MOI_mist*2.*np.pi/Prot  

#         # Initial Angular Velocity [1/d]
#         Omega = 2.*np.pi/Prot[0]

#         ## modulate the angular momentum loss and evolve Prot
#         ## go through each time step in the MIST table starting from t_DL
                
#         if Prot_init[j]>=3:
            
#             for i in range(ind0+1, (AGE_mist[AGE_mist < age_max]).shape[0]-1): 
    
#                 # define timestep [d]
#                 dt = (AGE_mist[i+1]-AGE_mist[i]).to(u.d) 
    
#                 # angular momentum loss rate, Eq. 2 from Garraffo+2018 [g cm^2 / d^2] 
#                 Jdot_dip = J0*u.g*u.cm**2. * Omega**3. * TAU_mist[i] 
    
#                 # Rossby number
#                 Ro = Prot[i]/TAU_mist[i]
    
#                 # complexity function n(Ro), Eq. 5 in Garraffo et al 2018
#                 n = a/Ro + b*Ro + n_min 
#                 # n = 0.03/Ro + 0.5*Ro + 1. #Gossage+2018
                
#                 # In the following, consider only stars with magnetic complexity 1 < n < 8
#                 if (n < 1) or (complexity == False):
#                     # Turn this on to ignore complexity: Skumanich evolution
#                     n = 1
#                 elif n > 8:
#                     n = 8 
    
#                 # Magnetic complexity modulation of AML Q_J(n), Eq. 4 in Garraffo et al 2018
#                 B = 100. #Gauss
#                 QJ = 4.05*np.exp(-1.4*n)# + (n[0]-1.)/(60.*B*n[0]) #2nd term becomes only important for n > 7 # dimensionless
    
#                 # Modulate angular momentum loss
#                 Jdot =  QJ*Jdot_dip
#                 Omega = J[min(abs(i), abs(i-1))]/MOI_mist[i]
#                 Prot[min(i+1, len(AGE_mist)-1)] = 2.*np.pi/Omega
#                 Omega = Omega - dt*Jdot/MOI_mist[i]
#                 J[i] = MOI_mist[i]* Omega
                
#             # stellar mass [Msun], time [Myr], evolved Prot [d]
#             Prot[Prot < 0.] = np.nan
#             Prot_evol.append([Mstar[j], AGE_mist.to(u.Myr), Prot])

#             # interpolate all evolved Prots along the same age-array for the direct comparison
#             # (as different masses were evolved for different long times in the MIST tables)
#             if j==0:
#                 # define a common age-array for the interpolation
#                 age_zero = AGE_mist/1e6
#             spl_Prot = InterpolatedUnivariateSpline((AGE_mist.to(u.Myr)).value, Prot.value, ext=0)
#             Prot_interp.append(spl_Prot(age_zero))
            
#         if Prot_init[j]<3:
        
#             if Mstar[j].value!=MASSES[0] and Mstar[j].value!=MASSES[1]: # for all masses but 0.1 and 0.15 M_sol
            
#                 for i in range(ind0+1, (AGE_mist[AGE_mist < age_max]).shape[0]-1): 
        
#                     # define timestep [d]
#                     dt = (AGE_mist[i+1]-AGE_mist[i]).to(u.d) 
        
#                     # angular momentum loss rate, Eq. 2 from Garraffo+2018 [g cm^2 / d^2] 
#                     Jdot_dip = J0*u.g*u.cm**2. * Omega**3. * TAU_mist[i] 
        
#                     # Rossby number
#                     Ro = Prot[i]/TAU_mist[i]
        
#                     # complexity function n(Ro), Eq. 5 in Garraffo et al 2018
#                     n = a/Ro + b*Ro + n_min 
#                     # n = 0.03/Ro + 0.5*Ro + 1. #Gossage+2018
                    
#                     # In the following, consider only stars with magnetic complexity 1 < n < 8
#                     if (n < 1) or (complexity == False):
#                         # Turn this on to ignore complexity: Skumanich evolution
#                         n = 1
#                     elif n > 8:
#                         n = 8 
        
#                     # Magnetic complexity modulation of AML Q_J(n), Eq. 4 in Garraffo et al 2018
#                     B = 100. #Gauss
#                     QJ = 4.05*np.exp(-1.4*n)# + (n[0]-1.)/(60.*B*n[0]) #2nd term becomes only important for n > 7 # dimensionless
        
#                     # Modulate angular momentum loss
#                     Jdot =  QJ*Jdot_dip
#                     Omega = J[min(abs(i), abs(i-1))]/MOI_mist[i]
#                     Prot[min(i+1, len(AGE_mist)-1)] = 2.*np.pi/Omega
#                     Omega = Omega - dt*Jdot/MOI_mist[i]
#                     J[i] = MOI_mist[i]* Omega
                    
#                 # stellar mass [Msun], time [Myr], evolved Prot [d]
#                 Prot[Prot < 0.] = np.nan
#                 Prot_evol.append([Mstar[j], AGE_mist.to(u.Myr), Prot])

#                 # interpolate all evolved Prots along the same age-array for the direct comparison
#                 # (as different masses were evolved for different long times in the MIST tables)
#                 if j==0:
#                     # define a common age-array for the interpolation
#                     age_zero = AGE_mist/1e6
#                 spl_Prot = InterpolatedUnivariateSpline((AGE_mist.to(u.Myr)).value, Prot.value, ext=0)
#                 Prot_interp.append(spl_Prot(age_zero))

#             if Mstar[j].value==MASSES[1]: # for mass = 0.15 M_sol
#                 # age at which code should stop
#                 c1=AGE_mist.value
#                 crit_age=550*1e6
#                 d1=find_nearest(c1, crit_age)
#                 ind=int(np.where(c1==d1)[0])
                
#                 with open(f'/home/farah/Documents/Project/Data/Prot_i_3d_evolution_0.15M','rb') as f: file = pickle.load(f)
                
#                 Prot_i_3_015=file*u.d 
                
#                 for i in range(ind0+1, ind): 
        
#                     # define timestep [d]
#                     dt = (AGE_mist[i+1]-AGE_mist[i]).to(u.d) 
        
#                     # angular momentum loss rate, Eq. 2 from Garraffo+2018 [g cm^2 / d^2] 
#                     Jdot_dip = J0*u.g*u.cm**2. * Omega**3. * TAU_mist[i] 
        
#                     # Rossby number
#                     Ro = Prot[i]/TAU_mist[i]
        
#                     # complexity function n(Ro), Eq. 5 in Garraffo et al 2018
#                     n = a/Ro + b*Ro + n_min 
#                     # n = 0.03/Ro + 0.5*Ro + 1. #Gossage+2018
                    
#                     # In the following, consider only stars with magnetic complexity 1 < n < 8
#                     if (n < 1) or (complexity == False):
#                         # Turn this on to ignore complexity: Skumanich evolution
#                         n = 1
#                     elif n > 8:
#                         n = 8 
        
#                     # Magnetic complexity modulation of AML Q_J(n), Eq. 4 in Garraffo et al 2018
#                     B = 100. #Gauss
#                     QJ = 4.05*np.exp(-1.4*n)# + (n[0]-1.)/(60.*B*n[0]) #2nd term becomes only important for n > 7 # dimensionless
        
#                     # Modulate angular momentum loss
#                     Jdot =  QJ*Jdot_dip
#                     Omega = J[min(abs(i), abs(i-1))]/MOI_mist[i]
#                     Prot[min(i+1, len(AGE_mist)-1)] = 2.*np.pi/Omega
#                     Omega = Omega - dt*Jdot/MOI_mist[i]
#                     J[i] = MOI_mist[i]* Omega
                    
#                 # stellar mass [Msun], time [Myr], evolved Prot [d]
#                 Prot[Prot < 0.] = np.nan
#                 Prot[ind:248]=Prot_i_3_015[ind:248]
#                 Prot_evol.append([Mstar[j], AGE_mist.to(u.Myr), Prot])
    
#                 # interpolate all evolved Prots along the same age-array for the direct comparison
#                 # (as different masses were evolved for different long times in the MIST tables)
#                 if j==0:
#                     # define a common age-array for the interpolation
#                     age_zero = AGE_mist/1e6
#                 spl_Prot = InterpolatedUnivariateSpline((AGE_mist.to(u.Myr)).value, Prot.value, ext=0)
#                 Prot_interp.append(spl_Prot(age_zero))
                
#             if Mstar[j].value==MASSES[0]: # for mass = 0.1 M_sol
#                 # age at which code should stop
#                 c1=AGE_mist.value
#                 crit_age=390*1e6
#                 d1=find_nearest(c1, crit_age)
#                 ind=int(np.where(c1==d1)[0])
                
#                 with open(f'/home/farah/Documents/Project/Data/Prot_i_3d_evolution_0.1M','rb') as f: file = pickle.load(f)
                
#                 Prot_i_3_01=file*u.d 
                
#                 for i in range(ind0+1, ind): 
        
#                     # define timestep [d]
#                     dt = (AGE_mist[i+1]-AGE_mist[i]).to(u.d) 
        
#                     # angular momentum loss rate, Eq. 2 from Garraffo+2018 [g cm^2 / d^2] 
#                     Jdot_dip = J0*u.g*u.cm**2. * Omega**3. * TAU_mist[i] 
        
#                     # Rossby number
#                     Ro = Prot[i]/TAU_mist[i]
        
#                     # complexity function n(Ro), Eq. 5 in Garraffo et al 2018
#                     n = a/Ro + b*Ro + n_min 
#                     # n = 0.03/Ro + 0.5*Ro + 1. #Gossage+2018
                    
#                     # In the following, consider only stars with magnetic complexity 1 < n < 8
#                     if (n < 1) or (complexity == False):
#                         # Turn this on to ignore complexity: Skumanich evolution
#                         n = 1
#                     elif n > 8:
#                         n = 8 
        
#                     # Magnetic complexity modulation of AML Q_J(n), Eq. 4 in Garraffo et al 2018
#                     B = 100. #Gauss
#                     QJ = 4.05*np.exp(-1.4*n)# + (n[0]-1.)/(60.*B*n[0]) #2nd term becomes only important for n > 7 # dimensionless
        
#                     # Modulate angular momentum loss
#                     Jdot =  QJ*Jdot_dip
#                     Omega = J[min(abs(i), abs(i-1))]/MOI_mist[i]
#                     Prot[min(i+1, len(AGE_mist)-1)] = 2.*np.pi/Omega
#                     Omega = Omega - dt*Jdot/MOI_mist[i]
#                     J[i] = MOI_mist[i]* Omega
                    
#                 # stellar mass [Msun], time [Myr], evolved Prot [d]
#                 Prot[Prot < 0.] = np.nan
#                 Prot[ind:239]=Prot_i_3_01[ind:239]
#                 Prot_evol.append([Mstar[j], AGE_mist.to(u.Myr), Prot])
    
#                 # interpolate all evolved Prots along the same age-array for the direct comparison
#                 # (as different masses were evolved for different long times in the MIST tables)
#                 if j==0:
#                     # define a common age-array for the interpolation
#                     age_zero = AGE_mist/1e6
#                 spl_Prot = InterpolatedUnivariateSpline((AGE_mist.to(u.Myr)).value, Prot.value, ext=0)
#                 Prot_interp.append(spl_Prot(age_zero))

#     Prot_evol = np.array(Prot_evol)   
#     Prot_interp = np.array(Prot_interp)
    
#     return Prot_evol, age_zero, Prot_interp


#%% Smooth plots for stars with mass 0.1 and 0.15 M_sol

'Smooth plots for stars with mass 0.1 and 0.15 M_sol'


import numpy as np
import astropy.units as u
import astropy.constants as const
import sys

sys.path.append( '/home/farah/Documents/Project/Data/MIST_tables/' )

#Create mass array to find nearest neighbour
MASSES = np.arange(0.1,2.05, 0.05) #all the MIST masses we have

#Find nearest neighbour in array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def load_mist_tables(Mstar=1., filepath='/home/farah/Documents/Project/Data/MIST_tables/'):
        """
        Load in the MIST tables.
        Mstar: Stellar masses in units of solar masses
        filepath: Path where the MIST tables are stored
        """
        import read_mist_models
        import astropy.units as u

        print(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)))

        eep = read_mist_models.EEP(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)), verbose=False)
        AGE_mist = eep.eeps['star_age']*u.yr # stellar age in years
        TAU_mist = (eep.eeps['conv_env_turnover_time_g']*u.s).to(u.d) # convective turnover time in days
        MOI_mist = eep.eeps['moment_of_inertia']*u.g*u.cm**2. # moment of inertia in cgs

        return AGE_mist, TAU_mist, MOI_mist


def define_cgs_constants():
    import astropy.constants as const
    
    M_sun = const.M_sun.cgs     # solar mass [gr]
    R_sun = const.R_sun.cgs     # solar radius [cm]
    G = const.G.cgs             # Newton's constant in cm^3 g^-1 s^-2

    return M_sun, R_sun, G


def spin_down_evol(Prot_init, Mstar_init, t0=1., tdisc=13., a=0.02, b=2., J0=1e41, n_min=1., complexity=True):

    """
    Evolves the stellar rotation period based on the spin-evolution model published in Garraffo et al. (2015, 2016, 2018)
    Prot_init: Starting rotation periods
    Mstar_init: Corresponding stellar masses
    t0: Time in Myr at which the rotational evolution should be started
    a, b, J0, n_min: Prefactors for the n-vs-Ro relation in Garraffo+2018 
    """
    import astropy.units as u
    from scipy.interpolate import InterpolatedUnivariateSpline
    from scipy.interpolate import UnivariateSpline
    #from scipy.signal import savgol_filter
    
    M_sun, R_Sun, G = define_cgs_constants()

    # DEFINE PARAMETERS
    Prot0 = Prot_init*u.d               # initial rotation periods
    Mstar = Mstar_init*u.Msun           # initial stellar masses
    N = len(Prot_init)                  # sample size (number of stars)

    age_max = 2e10*u.yr                 # oldest cluster's age [million years]
    tf = age_max.to(u.d)                # final time (cluster's age)  [d]
    ti = (t0*u.Myr).to(u.d)             # initial age [d]  --> age of ONC

    # define empty arrays to which the data will be saved
    Prot_evol = []
    Prot_interp = []

    tdisc *= u.Myr
    t0 *= u.Myr

    ### EVOLVE ALGORITHM
    for j in range(N): 
        ### LOAD MIST table
        AGE_mist, TAU_mist, MOI_mist = load_mist_tables(Mstar=Mstar[j].value) 
                
        # time at which each star will start its rotational evolution [Myrs]
        t_DL = (tdisc-t0).to(u.yr)

        # find index in MIST table at which t = t_DL
        ind0 = min(range(len(AGE_mist)), key= lambda x:abs((AGE_mist[x]-t_DL).value))

        # Rotation Periods [d]
        Prot = np.tile(Prot0[j], len(AGE_mist))

        # Angular Momentum [g cm^2/d]
        J = MOI_mist*2.*np.pi/Prot  

        # Initial Angular Velocity [1/d]
        Omega = 2.*np.pi/Prot[0]

        ## modulate the angular momentum loss and evolve Prot
        ## go through each time step in the MIST table starting from t_DL
        for i in range(ind0+1, (AGE_mist[AGE_mist < age_max]).shape[0]-1): 
    
            # define timestep [d]
            dt = (AGE_mist[i+1]-AGE_mist[i]).to(u.d) 

            # angular momentum loss rate, Eq. 2 from Garraffo+2018 [g cm^2 / d^2] 
            Jdot_dip = J0*u.g*u.cm**2. * Omega**3. * TAU_mist[i] 

            # Rossby number
            Ro = Prot[i]/TAU_mist[i]

            # complexity function n(Ro), Eq. 5 in Garraffo et al 2018
            n = a/Ro + b*Ro + n_min 
            # n = 0.03/Ro + 0.5*Ro + 1. #Gossage+2018
            
            # In the following, consider only stars with magnetic complexity 1 < n < 8
            if (n < 1) or (complexity == False):
                # Turn this on to ignore complexity: Skumanich evolution
                n = 1
            elif n > 8:
                n = 8 

            # Magnetic complexity modulation of AML Q_J(n), Eq. 4 in Garraffo et al 2018
            B = 100. #Gauss
            QJ = 4.05*np.exp(-1.4*n)# + (n[0]-1.)/(60.*B*n[0]) #2nd term becomes only important for n > 7 # dimensionless

            # Modulate angular momentum loss
            Jdot =  QJ*Jdot_dip
            Omega = J[min(abs(i), abs(i-1))]/MOI_mist[i]
            Prot[min(i+1, len(AGE_mist)-1)] = 2.*np.pi/Omega
            Omega = Omega - dt*Jdot/MOI_mist[i]
            J[i] = MOI_mist[i]* Omega

        # stellar mass [Msun], time [Myr], evolved Prot [d]
        Prot[Prot < 0.] = np.nan
        Prot_evol.append([Mstar[j], AGE_mist.to(u.Myr), Prot])
        
        if Mstar[j].value==MASSES[0] and Prot_init[j]<3: # for mass = 0.1 M_sol and Prot < 3 days, we smooth the plot 
            # age at which code should stop
            c1=AGE_mist.value
            crit_age=390*1e6 #Age at which stars start to behave weirdly
            d1=find_nearest(c1, crit_age)
            ind=int(np.where(c1==d1)[0]) #Find index of value in AGE_mist array closest to the critical age
            
            #Open evolution of P_rot,i = 3 days
            with open(f'/home/farah/Documents/Project/Data/Prot_i_3d_evolution_0.1M','rb') as f: file = pickle.load(f)
            Prot_i_3_01=file*u.d 
            
            Prot[ind:]=Prot_i_3_01[ind:] #Smooth plot by using the evolution of P_rot,i = 3 days after the critical age
            
        if Mstar[j].value==MASSES[1] and Prot_init[j]<3: # for mass = 0.15 M_sol and Prot < 3 days, we smooth the plot
            # age at which code should stop
            c1=AGE_mist.value
            crit_age=550*1e6 #Age at which stars start to behave weirdly
            d1=find_nearest(c1, crit_age)
            ind=int(np.where(c1==d1)[0]) #Find index of value in AGE_mist array closest to the critical age
            
            #Open evolution of P_rot,i = 3 days
            with open(f'/home/farah/Documents/Project/Data/Prot_i_3d_evolution_0.15M','rb') as f: file = pickle.load(f)
            Prot_i_3_015=file*u.d 
            
            Prot[ind:]=Prot_i_3_015[ind:] #Smooth plot by using the evolution of P_rot,i = 3 days after the critical age


        # interpolate all evolved Prots along the same age-array for the direct comparison
        # (as different masses were evolved for different long times in the MIST tables)
        if j==0:
            # define a common age-array for the interpolation
            age_zero = AGE_mist/1e6
        spl_Prot = InterpolatedUnivariateSpline((AGE_mist.to(u.Myr)).value, Prot.value, ext=0)
        Prot_interp.append(spl_Prot(age_zero))
        
    Prot_evol = np.array(Prot_evol)   
    Prot_interp = np.array(Prot_interp)
    
    return Prot_evol, age_zero, Prot_interp


#%% Smooth plots for stars with mass 0.1 and 0.15 M_sol - Better

'Smooth plots for stars with mass 0.1 and 0.15 M_sol - Better'


import numpy as np
import astropy.units as u
import astropy.constants as const
import sys

sys.path.append( '/home/farah/Documents/Project/Data/MIST_tables/' )

#Create mass array to find nearest neighbour
MASSES = np.arange(0.1,2.05, 0.05) #all the MIST masses we have

#Find nearest neighbour in array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def load_mist_tables(Mstar=1., filepath='/home/farah/Documents/Project/Data/MIST_tables/'):
        """
        Load in the MIST tables.
        Mstar: Stellar masses in units of solar masses
        filepath: Path where the MIST tables are stored
        """
        import read_mist_models
        import astropy.units as u

        print(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)))

        eep = read_mist_models.EEP(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)), verbose=False)
        AGE_mist = eep.eeps['star_age']*u.yr # stellar age in years
        TAU_mist = (eep.eeps['conv_env_turnover_time_g']*u.s).to(u.d) # convective turnover time in days
        MOI_mist = eep.eeps['moment_of_inertia']*u.g*u.cm**2. # moment of inertia in cgs

        return AGE_mist, TAU_mist, MOI_mist


def define_cgs_constants():
    import astropy.constants as const
    
    M_sun = const.M_sun.cgs     # solar mass [gr]
    R_sun = const.R_sun.cgs     # solar radius [cm]
    G = const.G.cgs             # Newton's constant in cm^3 g^-1 s^-2

    return M_sun, R_sun, G


def spin_down_evol(Prot_init, Mstar_init, t0=1., tdisc=13., a=0.02, b=2., J0=1e41, n_min=1., complexity=True):

    """
    Evolves the stellar rotation period based on the spin-evolution model published in Garraffo et al. (2015, 2016, 2018)
    Prot_init: Starting rotation periods
    Mstar_init: Corresponding stellar masses
    t0: Time in Myr at which the rotational evolution should be started
    a, b, J0, n_min: Prefactors for the n-vs-Ro relation in Garraffo+2018 
    """
    import astropy.units as u
    from scipy.interpolate import InterpolatedUnivariateSpline
    from scipy.interpolate import UnivariateSpline
    #from scipy.signal import savgol_filter
    
    M_sun, R_Sun, G = define_cgs_constants()

    # DEFINE PARAMETERS
    Prot0 = Prot_init*u.d               # initial rotation periods
    Mstar = Mstar_init*u.Msun           # initial stellar masses
    N = len(Prot_init)                  # sample size (number of stars)

    age_max = 2e10*u.yr                 # oldest cluster's age [million years]
    tf = age_max.to(u.d)                # final time (cluster's age)  [d]
    ti = (t0*u.Myr).to(u.d)             # initial age [d]  --> age of ONC

    # define empty arrays to which the data will be saved
    Prot_evol = []
    Prot_interp = []

    tdisc *= u.Myr
    t0 *= u.Myr

    ### EVOLVE ALGORITHM
    for j in range(N): 
        ### LOAD MIST table
        AGE_mist, TAU_mist, MOI_mist = load_mist_tables(Mstar=Mstar[j].value) 
                
        # time at which each star will start its rotational evolution [Myrs]
        t_DL = (tdisc-t0).to(u.yr)

        # find index in MIST table at which t = t_DL
        ind0 = min(range(len(AGE_mist)), key= lambda x:abs((AGE_mist[x]-t_DL).value))
        
        # find index at which we want to stop the rotational evolution
        prot_ind = (AGE_mist[AGE_mist < age_max]).shape[0]

        # Rotation Periods [d]
        Prot = np.tile(Prot0[j], prot_ind)

        # Angular Momentum [g cm^2/d]
        J = MOI_mist[:prot_ind]*2.*np.pi/Prot  

        # Initial Angular Velocity [1/d]
        Omega = 2.*np.pi/Prot[0]

        ## modulate the angular momentum loss and evolve Prot
        ## go through each time step in the MIST table starting from t_DL
        for i in range(ind0+1, prot_ind-1): 
    
            # define timestep [d]
            dt = (AGE_mist[i+1]-AGE_mist[i]).to(u.d) 

            # angular momentum loss rate, Eq. 2 from Garraffo+2018 [g cm^2 / d^2] 
            Jdot_dip = J0*u.g*u.cm**2. * Omega**3. * TAU_mist[i] 

            # Rossby number
            Ro = Prot[i]/TAU_mist[i]

            # complexity function n(Ro), Eq. 5 in Garraffo et al 2018
            n = a/Ro + b*Ro + n_min 
            # n = 0.03/Ro + 0.5*Ro + 1. #Gossage+2018
            
            # In the following, consider only stars with magnetic complexity 1 < n < 8
            if (n < 1) or (complexity == False):
                # Turn this on to ignore complexity: Skumanich evolution
                n = 1
            elif n > 8:
                n = 8 

            # Magnetic complexity modulation of AML Q_J(n), Eq. 4 in Garraffo et al 2018
            B = 100. #Gauss
            QJ = 4.05*np.exp(-1.4*n)# + (n[0]-1.)/(60.*B*n[0]) #2nd term becomes only important for n > 7 # dimensionless

            # Modulate angular momentum loss
            Jdot =  QJ*Jdot_dip
            Omega = J[min(abs(i), abs(i-1))]/MOI_mist[i]
            Prot[min(i+1, len(AGE_mist)-1)] = 2.*np.pi/Omega
            Omega = Omega - dt*Jdot/MOI_mist[i]
            J[i] = MOI_mist[i]* Omega

        # stellar mass [Msun], time [Myr], evolved Prot [d]
        Prot[Prot < 0.] = np.nan
        Prot_evol.append([Mstar[j], AGE_mist[:prot_ind].to(u.Myr), Prot])
        
        if Mstar[j].value==MASSES[0] and Prot_init[j]<3: # for mass = 0.1 M_sol and Prot < 3 days, we smooth the plot 
            # age at which code should stop
            c1=AGE_mist.value
            crit_age=390*1e6 #Age at which stars start to behave weirdly
            d1=find_nearest(c1, crit_age)
            ind=int(np.where(c1==d1)[0]) #Find index of value in AGE_mist array closest to the critical age
            
            #Open evolution of P_rot,i = 3 days
            with open(f'/home/farah/Documents/Project/Data/Prot_i_3d_evolution_0.1M','rb') as f: file = pickle.load(f)
            Prot_i_3_01=file*u.d 
            
            Prot[ind:]=Prot_i_3_01[ind:] #Smooth plot by using the evolution of P_rot,i = 3 days after the critical age
            
        if Mstar[j].value==MASSES[1] and Prot_init[j]<3: # for mass = 0.15 M_sol and Prot < 3 days, we smooth the plot
            # age at which code should stop
            c1=AGE_mist.value
            crit_age=550*1e6 #Age at which stars start to behave weirdly
            d1=find_nearest(c1, crit_age)
            ind=int(np.where(c1==d1)[0]) #Find index of value in AGE_mist array closest to the critical age
            
            #Open evolution of P_rot,i = 3 days
            with open(f'/home/farah/Documents/Project/Data/Prot_i_3d_evolution_0.15M','rb') as f: file = pickle.load(f)
            Prot_i_3_015=file*u.d 
            
            Prot[ind:]=Prot_i_3_015[ind:] #Smooth plot by using the evolution of P_rot,i = 3 days after the critical age


        # interpolate all evolved Prots along the same age-array for the direct comparison
        # (as different masses were evolved for different long times in the MIST tables)
        if j==0:
            # define a common age-array for the interpolation
            age_zero = AGE_mist/1e6
        spl_Prot = InterpolatedUnivariateSpline((AGE_mist[:prot_ind].to(u.Myr)).value, Prot.value, ext=0)
        Prot_interp.append(spl_Prot(age_zero))
        
    Prot_evol = np.array(Prot_evol)   
    Prot_interp = np.array(Prot_interp)
    
    return Prot_evol, age_zero, Prot_interp




#%% 'Unsmooth plots for stars with mass 0.1 and 0.15 M_sol'

'Unsmooth plots for stars with mass 0.1 and 0.15 M_sol'


import numpy as np
import astropy.units as u
import astropy.constants as const
import sys

sys.path.append( '/home/farah/Documents/Project/Data/MIST_tables/' )

#Create mass array to find nearest neighbour
MASSES = np.arange(0.1,2.05, 0.05) #all the MIST masses we have

#Find nearest neighbour in array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def load_mist_tables(Mstar=1., filepath='/home/farah/Documents/Project/Data/MIST_tables/'):
        """
        Load in the MIST tables.
        Mstar: Stellar masses in units of solar masses
        filepath: Path where the MIST tables are stored
        """
        import read_mist_models
        import astropy.units as u

        print(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)))

        eep = read_mist_models.EEP(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)), verbose=False)
        AGE_mist = eep.eeps['star_age']*u.yr # stellar age in years
        TAU_mist = (eep.eeps['conv_env_turnover_time_g']*u.s).to(u.d) # convective turnover time in days
        MOI_mist = eep.eeps['moment_of_inertia']*u.g*u.cm**2. # moment of inertia in cgs

        return AGE_mist, TAU_mist, MOI_mist


def define_cgs_constants():
    import astropy.constants as const
    
    M_sun = const.M_sun.cgs     # solar mass [gr]
    R_sun = const.R_sun.cgs     # solar radius [cm]
    G = const.G.cgs             # Newton's constant in cm^3 g^-1 s^-2

    return M_sun, R_sun, G


def spin_down_evol(Prot_init, Mstar_init, t0=1., tdisc=13., a=0.02, b=2., J0=1e41, n_min=1., complexity=True):

    """
    Evolves the stellar rotation period based on the spin-evolution model published in Garraffo et al. (2015, 2016, 2018)
    Prot_init: Starting rotation periods
    Mstar_init: Corresponding stellar masses
    t0: Time in Myr at which the rotational evolution should be started
    a, b, J0, n_min: Prefactors for the n-vs-Ro relation in Garraffo+2018 
    """
    import astropy.units as u
    from scipy.interpolate import InterpolatedUnivariateSpline
    from scipy.interpolate import UnivariateSpline
    #from scipy.signal import savgol_filter
    
    M_sun, R_Sun, G = define_cgs_constants()

    # DEFINE PARAMETERS
    Prot0 = Prot_init*u.d               # initial rotation periods
    Mstar = Mstar_init*u.Msun           # initial stellar masses
    N = len(Prot_init)                  # sample size (number of stars)

    age_max = 2e10*u.yr                 # oldest cluster's age [million years]
    tf = age_max.to(u.d)                # final time (cluster's age)  [d]
    ti = (t0*u.Myr).to(u.d)             # initial age [d]  --> age of ONC

    # define empty arrays to which the data will be saved
    Prot_evol = []
    Prot_interp = []

    tdisc *= u.Myr
    t0 *= u.Myr

    ### EVOLVE ALGORITHM
    for j in range(N): 
        ### LOAD MIST table
        AGE_mist, TAU_mist, MOI_mist = load_mist_tables(Mstar=Mstar[j].value) 
                
        # time at which each star will start its rotational evolution [Myrs]
        t_DL = (tdisc-t0).to(u.yr)

        # find index in MIST table at which t = t_DL
        ind0 = min(range(len(AGE_mist)), key= lambda x:abs((AGE_mist[x]-t_DL).value))

        # Rotation Periods [d]
        Prot = np.tile(Prot0[j], len(AGE_mist))

        # Angular Momentum [g cm^2/d]
        J = MOI_mist*2.*np.pi/Prot  

        # Initial Angular Velocity [1/d]
        Omega = 2.*np.pi/Prot[0]

        ## modulate the angular momentum loss and evolve Prot
        ## go through each time step in the MIST table starting from t_DL
        for i in range(ind0+1, (AGE_mist[AGE_mist < age_max]).shape[0]-1): 
    
            # define timestep [d]
            dt = (AGE_mist[i+1]-AGE_mist[i]).to(u.d) 

            # angular momentum loss rate, Eq. 2 from Garraffo+2018 [g cm^2 / d^2] 
            Jdot_dip = J0*u.g*u.cm**2. * Omega**3. * TAU_mist[i] 

            # Rossby number
            Ro = Prot[i]/TAU_mist[i]

            # complexity function n(Ro), Eq. 5 in Garraffo et al 2018
            n = a/Ro + b*Ro + n_min 
            # n = 0.03/Ro + 0.5*Ro + 1. #Gossage+2018
            
            # In the following, consider only stars with magnetic complexity 1 < n < 8
            if (n < 1) or (complexity == False):
                # Turn this on to ignore complexity: Skumanich evolution
                n = 1
            elif n > 8:
                n = 8 

            # Magnetic complexity modulation of AML Q_J(n), Eq. 4 in Garraffo et al 2018
            B = 100. #Gauss
            QJ = 4.05*np.exp(-1.4*n)# + (n[0]-1.)/(60.*B*n[0]) #2nd term becomes only important for n > 7 # dimensionless

            # Modulate angular momentum loss
            Jdot =  QJ*Jdot_dip
            Omega = J[min(abs(i), abs(i-1))]/MOI_mist[i]
            Prot[min(i+1, len(AGE_mist)-1)] = 2.*np.pi/Omega
            Omega = Omega - dt*Jdot/MOI_mist[i]
            J[i] = MOI_mist[i]* Omega

        # stellar mass [Msun], time [Myr], evolved Prot [d]
        Prot[Prot < 0.] = np.nan
        Prot_evol.append([Mstar[j], AGE_mist.to(u.Myr), Prot])

        # interpolate all evolved Prots along the same age-array for the direct comparison
        # (as different masses were evolved for different long times in the MIST tables)
        if j==0:
            # define a common age-array for the interpolation
            age_zero = AGE_mist/1e6
        spl_Prot = InterpolatedUnivariateSpline((AGE_mist.to(u.Myr)).value, Prot.value, ext=0)
        Prot_interp.append(spl_Prot(age_zero))
        
    Prot_evol = np.array(Prot_evol)   
    Prot_interp = np.array(Prot_interp)
    
    return Prot_evol, age_zero, Prot_interp



#%% For a mass of 0.1 and 0.15 solar mass, smooth rotation period evolution plots

'For a mass of 0.1 and 0.15 solar mass, smooth rotation period evolution plots'

#Create initial rotation period data with hPer data

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import pickle

# define initial rotation periods
P0 = np.array([3])

mass = [0.15]

for i in mass:
    # define initial stellar masses
    M0 = np.ones((P0.shape))*find_nearest(MASSES, mass)
    
    # t0: time [Myr] at which the model should start 
    # tdisc: disk-locking time [Myr]
    Prot_evol, age_zero, Prot_interp = spin_down_evol(Prot_init=P0, 
                                                        Mstar_init=M0, 
                                                        t0=1., tdisc=13.)
    
    fig = plt.figure(figsize=(8, 5))
    fig.add_axes([0.1, 0.1, 0.8, 0.8]) #real plot
        
    plt.xlim([10,2.5*10**4])
    plt.ylim([ 0.01, 10000])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('P [days]', fontsize=15)
    plt.title(f'Stellar spindown plot for {i}$M_\odot$')
    for el in range(len(Prot_evol)):
        plt.plot( Prot_evol[el,1], Prot_evol[el,2])
       
    Prot_evolution=Prot_evol[el,2].value
    print(Prot_evolution)
    with open(f'/home/farah/Documents/Project/Data/Prot_i_3d_evolution_0.15M','wb') as f: pickle.dump(Prot_evolution, f) #Save data for 0.15 mass and 3 days Prot_i as a pickle file
    with open(f'/home/farah/Documents/Project/Data/Prot_i_3d_evolution_0.15M','rb') as f: arrayname1 = pickle.load(f)
    
    b=np.array_equal(Prot_evolution,arrayname1) #Sanity check
    print(b)

    
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)
    #plt.text(20, 80, r'Initial Mass={}'.format(1.0)+'$M_{\odot}$',  fontsize=20)
    
    plt.savefig('/home/farah/Documents/Project/Data/015M_Smooth_plot.png')

    
    
    
    
#%%

a=[2,3,4,5,6,78,100]
b=[2,4,55,66,77,88,99]

b[2:]=a[2:]

print(b)

