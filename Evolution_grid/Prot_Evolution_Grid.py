#%%
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

    age_max = 1e10*u.yr                 # oldest cluster's age [million years]
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
    
    return Prot_evol, age_zero, Prot_interp, spl_Prot

#%% Better code

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

    age_max = 1e10*u.yr                 # oldest cluster's age [million years]
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
        
        #Age at which the MIST table should stop depending on mass of star
        index=int(np.where(MASSES==Mstar_init[j])[0])
        stop_MIST=[20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,17469,15012,13013,11368,10000,7000,5500,4500,4000,3500]
        max_age=stop_MIST[index]*1e6*u.yr 
                
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
        for i in range(ind0+1, (AGE_mist[AGE_mist < max_age]).shape[0]-1): 
                
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
    
    return Prot_evol, age_zero, Prot_interp, spl_Prot

#%% Even better

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
        
        #Age at which the MIST table should stop depending on mass of star
        index=int(np.where(MASSES==Mstar_init[j])[0])
        stop_MIST=[20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,17469,15012,13013,11368,10000,7000,5500,4500,4000,3500]
        max_age=stop_MIST[index]*1e6*u.yr 
                
        # time at which each star will start its rotational evolution [Myrs]
        t_DL = (tdisc-t0).to(u.yr)

        # find index in MIST table at which t = t_DL
        ind0 = min(range(len(AGE_mist)), key= lambda x:abs((AGE_mist[x]-t_DL).value))
        
        # find index at which we want to stop the rotational evolution
        prot_ind = (AGE_mist[AGE_mist < max_age]).shape[0]

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
    
    return Prot_evol, age_zero, Prot_interp, spl_Prot

#%% Rotation period and age evolution of stars - Grid attempt

'Rotation period and age evolution of stars - Grid attempt'

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import pickle

mass = np.arange(0.1,1.3,0.05) #0.05 steps

name_mass = [0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1, 1.05, 1.1, 1.15, 1.2, 1.25] #exact masses

RotationP=np.arange(0.1,12.1,0.1) #0.1 steps
              
m_sample=len(mass) #How many masses we want to look at
p_sample=len(RotationP) #How many rotation periods we want to look at

#Empty data one unit arrays so they can be used by the model
Prot_new=np.zeros(1) 
M_new=np.zeros(1)

#Creating a loop to create the evolutionary data we need by taking each initial period and running the model for each mass
for c in range(m_sample): #For each mass
    M_new[0]=find_nearest(MASSES, mass[c]) #Put mass number in the single unit mass array
    a=load_mist_tables(Mstar=find_nearest(MASSES, mass[c])) #Look at how many time steps the MIST table has
    #Age at which the MIST table should stop depending on mass of star
    index=int(np.where(MASSES==M_new[0])[0])
    stop_MIST=[20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,17469,15012,13013,11368,10000,7000,5500,4500,4000,3500]
    max_age=stop_MIST[index]*1e6*u.yr  
    # find index at which we want to stop the rotational evolution
    prot_ind = (a[0][a[0] < max_age]).shape[0]
    Final_set=np.zeros((p_sample,2,len(a[0][:prot_ind]))) #Empty array with final data with rotation period evolution over time

    for i in range(p_sample): #For each initial rotation period
        Prot_new[0]=RotationP[i] #Put rotation period number in the single unit Prot array
        # t0: time [Myr] at which the model should start 
        # tdisc: disk-locking time [Myr]
        Prot_evol, age_zero, Prot_interp, spl_Prot = spin_down_evol(Prot_init=Prot_new, 
                                                                Mstar_init=M_new, 
                                                                t0=1., tdisc=13.)

        Final_set[i][0]=Prot_evol[0,1] #Evolution of time
        Final_set[i][1]=Prot_evol[0,2] #Evolution of rotation period
        
    'Scatter plot to study time steps'

    #Plot Age vs Rotation period
    fig = plt.figure(figsize=(15, 10))
    fig.add_axes([0.1, 0.1, 0.8, 0.8]) #real plot
            
    plt.xlim([10,2.5*10**4])
    plt.ylim([ 0.01, 10000])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('P [days]', fontsize=15)
    plt.title(f'Evolution of the rotation period of stars - {name_mass[c]}$M_\odot$', fontsize= 25)

    for el in range(p_sample):
        plt.scatter(Final_set[el][0], Final_set[el][1])
            
            
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)
       #plt.text(20, 80, r'Initial Mass={}'.format(1.0)+'$M_{\odot}$',  fontsize=20)
    
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/{name_mass[c]}_evolution','wb') as f: pickle.dump(Final_set, f) #Save data for one mass as a pickle file
    
    'Open pickle file'

    # with open(f'/home/farah/Documents/Project/Data/Evolution_grid/{name_mass[c]}_evolution','rb') as f: yes = pickle.load(f)
        
    # b=np.array_equal(Final_set,yes) #Sanity check
    # print(b)
    
    plt.savefig(f'/home/farah/Documents/Project/Data/Evolution_grid/Scatter_Time_evolution_study_{name_mass[c]}.png')


#print(Final_set) #Print final array with all data for one mass

#%% Open pickle file

'Open pickle file'

with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.1_evolution','rb') as f: arrayname1 = pickle.load(f)
    
b=np.array_equal(Final_set,arrayname1) #Sanity check
print(b)


#%% Scatter plot to study time steps

'Scatter plot to study time steps'

#Plot Age vs Rotation period
fig = plt.figure(figsize=(15, 10))
fig.add_axes([0.1, 0.1, 0.8, 0.8]) #real plot
        
plt.xlim([10,6*10**4])
plt.ylim([ 0.01, 300])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Age [Myrs]', fontsize=15)
plt.ylabel('P [days]', fontsize=15)
plt.title('Evolution of the rotation period of stars', fontsize= 25)

for el in range(p_sample):
    plt.scatter(Final_set[el][0], Final_set[el][1])
        
        
plt.tick_params(labelsize=14)
plt.tick_params(labelsize=14)
plt.tick_params(labelsize=14)
plt.tick_params(axis='y', which='both')
plt.tick_params(bottom=True, left=True, right=True)
plt.tick_params(which='both',labelright=True)
   #plt.text(20, 80, r'Initial Mass={}'.format(1.0)+'$M_{\odot}$',  fontsize=20)
   
#%%Stop of MIST table for each mass (0.1 - 1.25 m)

'Stop of MIST table for each mass (0.1 - 1.25 m)'

#stop_MIST=[9975.20925027,9900.60243034,9614.39463046,9986.37580404,9815.1120863,9957.58642756,9996.14347392,9856.07927672,9650.83978166,9733.36059057,9847.03027128,9927.86228216,9703.72491413,9987.46324024,9987.05593671,9969.62478491,9924.96314283,9936.39692114,9991.34041282,7000,5500,4500,4000,3500]
stop_MIST=[20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,17469,15012,13013,11368,10000,7000,5500,4500,4000,3500]
n_steps=[239,248,254,262,270,299,306,310,314,320,326,333,339,347,356,380,401,426,465]






