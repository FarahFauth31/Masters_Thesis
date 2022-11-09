import numpy as np
import astropy.units as u
import astropy.constants as const
import sys
import pickle

sys.path.append( '/File/path/where/you/have/the/MIST_tables/' ) #You can download MIST tables from https://github.com/cgarraffo/Spin-down-model/tree/master/Mist_data

#Create mass array with all the MIST masses we have
MASSES = np.arange(0.1,2.05, 0.05)

#Function to find nearest neighbour in an array
def find_nearest(array, value):
    """
    Find nearest neighbour in an array given a value
    array: Contains all values
    value: Variable from which we want to find the nearest neighbour in the array
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

#Function that loads all the MIST tables
def load_mist_tables(Mstar=1., filepath='/File/path/where/you/have/the/MIST_tables/'):
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

#Function that defines cgs constants
def define_cgs_constants():
    import astropy.constants as const
    
    M_sun = const.M_sun.cgs     # solar mass [gr]
    R_sun = const.R_sun.cgs     # solar radius [cm]
    G = const.G.cgs             # Newton's constant in cm^3 g^-1 s^-2

    return M_sun, R_sun, G

#Function that evolves the stellar rotation period based on the spin-evolution model published in Garraffo et al. (2015, 2016, 2018)
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

    # age_max = 1e10*u.yr                 # oldest cluster's age [million years]
    # tf = age_max.to(u.d)                # final time (cluster's age)  [d]
    # ti = (t0*u.Myr).to(u.d)             # initial age [d]  --> age of ONC
    
    stop_MIST=[20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,20000,17469,15012,13013,11368,10000,7000,5500,4500,4000,3500] # age at which we want to stop the rotation period evolution for each individual mass [million years]

    # define empty arrays to which the data will be saved
    Prot_evol = []
    Prot_interp = []

    tdisc *= u.Myr
    t0 *= u.Myr

    ### EVOLVE ALGORITHM
    for j in range(N): 
        ### LOAD MIST table
        AGE_mist, TAU_mist, MOI_mist = load_mist_tables(Mstar=Mstar[j].value) 
        
        #Age from MIST table at which the code should stop depending on mass of star
        index=int(np.where(MASSES==Mstar_init[j])[0]) #Find index of mass in the MASSES array
        max_age=stop_MIST[index]*1e6*u.yr #Use that index to find the age at which the code should stop depending on mass of star
                
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
        
        ### For masses m=0.1 and m=0.15 [Solar masses], the rotation period evolution behaves weirdly with initial rotation periods P_rot,i < 3 days. That's why we smooth the plots
        if Mstar[j].value==MASSES[0] and Prot_init[j]<3: # for mass = 0.1 M_sol and Prot < 3 days, we smooth the plot 
            # age at which code should stop
            c1=AGE_mist.value
            crit_age=390*1e6 #Age at which stars start to behave weirdly
            d1=find_nearest(c1, crit_age)
            ind=int(np.where(c1==d1)[0]) #Find index of value in AGE_mist array closest to the critical age
            
            #Open evolution of P_rot,i = 3 days (Create that file running the code for the Prot_init=np.array([3]) and Mstar_init=np.array([0.1]) and saving the ouput as a pickloe file)
            with open(f'/home/farah/Documents/Project/Data/Prot_i_3d_evolution_0.1M','rb') as f: file = pickle.load(f)
            Prot_i_3_01=file*u.d 
            
            Prot[ind:]=Prot_i_3_01[ind:] #Smooth plot by using the evolution of P_rot,i = 3 days after the critical age
           
        if Mstar[j].value==MASSES[1] and Prot_init[j]<3: # for mass = 0.15 M_sol and Prot < 3 days, we smooth the plot
            # age at which code should stop
            c1=AGE_mist.value
            crit_age=550*1e6 #Age at which stars start to behave weirdly
            d1=find_nearest(c1, crit_age)
            ind=int(np.where(c1==d1)[0]) #Find index of value in AGE_mist array closest to the critical age
            
            #Open evolution of P_rot,i = 3 days (Create that file running the code for the Prot_init=np.array([3]) and Mstar_init=np.array([0.15]) and saving the ouput as a pickloe file)
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
