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


#%%

'Rotation evolution of a single star'

#Create initial rotation period data with hPer data

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

data = pd.read_csv("/home/farah/Documents/Project/Data/hPer_Data.csv") #Read hPer csv file
star_data = pd.read_csv("/media/farah/T7 Shield/GaiaUniverseModel_0000.csv") #Read GUMS csv file

#Create an initial rotation period distribution from hPer data
Full_m=data.Mass #Select the data set you want, in this case the mass (y data)
Full_Prot=data.Per #Select the data set you want, in this case the rotation periods

#Create new lists. The mass data is limited between 0.3 < m < 1.2 solar masses and the rotation period data is limited between 0.01 < Prot < 12 days.
m=[]
Prot=[]

for i in range(len(Full_m)):
    if 0.3 <= Full_m[i] <= 1.2 and 0.01 <= Full_Prot[i] <= 12:
        m.append(Full_m[i])
        Prot.append(Full_Prot[i])
        
#Resample data limiting mass and rotation period again        
n_sample=1 #Sample number
step=1 #One data point at a time
c=0 #Count

#Empty resampled data arrays with sample number size
Prot_resampled=np.zeros(n_sample) 
M_resampled=np.zeros(n_sample)

#Choosing the resampled data points that fall inside the mass and rotation period boundaries we set
while c < n_sample:
    values = np.vstack([Prot, m]) # create 2D array that contains the properties you want to resample
    kde_ProtMass = stats.gaussian_kde(values) # calculate 2D KDE of Prot-vs-Mstar
    
    #Create a resampled data point
    re_Prot = kde_ProtMass.resample(step)[0,:]
    re_M = kde_ProtMass.resample(step)[1,:]
    
    if 0.3 <= re_M <= 1.2 and 0.01 <= re_Prot <= 12:
        M_resampled[c]=re_M
        Prot_resampled[c]=re_Prot
        c+=1

#Mass of the object we are studying using star data
n_star=1
solar_mass=star_data.mass[n_star-1]
used_mass=find_nearest(MASSES, solar_mass)
M0 = np.ones((Prot_resampled.shape))*used_mass


# t0: time [Myr] at which the model should start 
# tdisc: disk-locking time [Myr]
Prot_evol, age_zero, Prot_interp, spl_Prot = spin_down_evol(Prot_init=Prot_resampled, 
                                                        Mstar_init=M0, 
                                                        t0=1., tdisc=13.)
#Select the rotation period for the age of the star
import numpy as np
#Select the age step closest to star's age
def closest(lst, K):
      
     idx = (np.abs(lst - K)).argmin()
     return lst[idx],idx #Time step closest to star's age, index of the age in the array

#K=1e8*u.yr #Age of star
K=star_data.age[n_star-1]*(10**9)*u.yr #Age of star (it is given in Gyr)
Calculated_Age,idx=closest(Prot_evol[0,1], K)
Calculated_Prot=Prot_evol[0,2][idx]

print(Calculated_Prot)

#Plot Age vs Rotation period
fig = plt.figure(figsize=(15, 10))
fig.add_axes([0.1, 0.1, 0.8, 0.8]) #real plot
    
plt.xlim([10,1*10**4])
plt.ylim([ 0.01, 100])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Age [Myrs]', fontsize=15)
plt.ylabel('P [days]', fontsize=15)
plt.scatter(Calculated_Age,Calculated_Prot, c='black', marker='$\odot$', s=500)
plt.title(f'Evolution of the rotation period of {solar_mass}$M_\odot$ star', fontsize= 20)
for el in range(len(Prot_evol)):
    plt.plot( Prot_evol[el,1], Prot_evol[el,2], c='orange')
    
    
plt.tick_params(labelsize=14)
plt.tick_params(labelsize=14)
plt.tick_params(labelsize=14)
plt.tick_params(axis='y', which='both')
plt.tick_params(bottom=True, left=True, right=True)
plt.tick_params(which='both',labelright=True)
#plt.text(20, 80, r'Initial Mass={}'.format(1.0)+'$M_{\odot}$',  fontsize=20)

plt.savefig('/home/farah/Documents/Project/Data/Star_Evolution_Prot_Approx.png')

#%%

import numpy as np
def closest(lst, K):
      
     idx = (np.abs(lst - K)).argmin()
     return lst[idx],idx

K=1e8*u.yr
lst,idx=closest(Prot_evol[0,1], K)
Calculated_Prot=Prot_evol[0,2][idx]
print(lst,idx)
print(Prot_evol[0,2][idx])
print(Prot_evol[0,1])
print(Prot_evol[0,2])

#%%

'Rotation evolution of a single star'

#Create initial rotation period data with hPer data

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

data = pd.read_csv("/home/farah/Documents/Project/Data/hPer_Data.csv") #Read hPer csv file6.24087732e+04
star_data = pd.read_csv("/media/farah/T7 Shield/GaiaUniverseModel_0000.csv") #Read GUMS csv file

#Create an initial rotation period distribution from hPer data
Full_m=data.Mass #Select the data set you want, in this case the mass (y data)
Full_Prot=data.Per #Select the data set you want, in this case the rotation periods

#Create new lists. The mass data is limited between 0.3 < m < 1.2 solar masses and the rotation period data is limited between 0.01 < Prot < 12 days.
m=[]
Prot=[]

for i in range(len(Full_m)):
    if 0.3 <= Full_m[i] <= 1.2 and 0.01 <= Full_Prot[i] <= 12:
        m.append(Full_m[i])
        Prot.append(Full_Prot[i])
        
#Resample data limiting mass and rotation period again        
n_sample=1 #Sample number
step=1 #One data point at a time
c=0 #Count

#Empty resampled data arrays with sample number size
Prot_resampled=np.zeros(n_sample) 
M_resampled=np.zeros(n_sample)

#Choosing the resampled data points that fall inside the mass and rotation period boundaries we set
while c < n_sample:
    values = np.vstack([Prot, m]) # create 2D array that contains the properties you want to resample
    kde_ProtMass = stats.gaussian_kde(values) # calculate 2D KDE of Prot-vs-Mstar
    
    #Create a resampled data point
    re_Prot = kde_ProtMass.resample(step)[0,:]
    re_M = kde_ProtMass.resample(step)[1,:]
    
    if 0.3 <= re_M <= 1.2 and 0.01 <= re_Prot <= 12:
        M_resampled[c]=re_M
        Prot_resampled[c]=re_Prot
        c+=1

#Mass of the object we are studying using star data
n_star=1
solar_mass=star_data.mass[n_star-1]
used_mass=find_nearest(MASSES, solar_mass)
M0 = np.ones((Prot_resampled.shape))*used_mass


# t0: time [Myr] at which the model should start 
# tdisc: disk-locking time [Myr]
Prot_evol, age_zero, Prot_interp, spl_Prot = spin_down_evol(Prot_init=Prot_resampled, 
                                                        Mstar_init=M0, 
                                                        t0=1., tdisc=13.)
#Select the rotation period for the age of the star
age=star_data.age[n_star-1]*(10**3) #Age of star (it is given in Gyr and we pass it to Myr)
Calculated_Prot=spl_Prot(age)

#Plot Age vs Rotation period
fig = plt.figure(figsize=(15, 10))
fig.add_axes([0.1, 0.1, 0.8, 0.8]) #real plot
    
plt.xlim([10,1*10**4])
plt.ylim([ 0.01, 100])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Age [Myrs]', fontsize=15)
plt.ylabel('P [days]', fontsize=15)
plt.scatter(age,Calculated_Prot, c='black', marker='$\odot$', s=500)
plt.title(f'Evolution of the rotation period of {solar_mass}$M_\odot$ star', fontsize= 20)
for el in range(len(Prot_evol)):
    plt.plot( Prot_evol[el,1], Prot_evol[el,2], c='orange')
    
    
plt.tick_params(labelsize=14)
plt.tick_params(labelsize=14)
plt.tick_params(labelsize=14)
plt.tick_params(axis='y', which='both')
plt.tick_params(bottom=True, left=True, right=True)
plt.tick_params(which='both',labelright=True)
#plt.text(20, 80, r'Initial Mass={}'.format(1.0)+'$M_{\odot}$',  fontsize=20)

plt.savefig('/home/farah/Documents/Project/Data/Star_Evolution_Prot.png')

print('This star with a mass of ',solar_mass,' solar masses has an age of ',age,' Myrs. The calculated rotation period for that age is ',Calculated_Prot,' days.')


#%%

'Individual star rotation evolution plots and the calculated Prot at individual ages'

#Create initial rotation period data with hPer data

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

data = pd.read_csv("/home/farah/Documents/Project/Data/hPer_Data.csv") #Read hPer csv file6.24087732e+04
star_data = pd.read_csv("/media/farah/T7 Shield/GaiaUniverseModel_0000.csv") #Read GUMS csv file

#Create an initial rotation period distribution from hPer data
Full_m=data.Mass #Select the data set you want, in this case the mass (y data)
Full_Prot=data.Per #Select the data set you want, in this case the rotation periods

#Create new lists. The mass data is limited between 0.3 < m < 1.2 solar masses and the rotation period data is limited between 0.01 < Prot < 12 days.
m=[]
Prot=[]

for i in range(len(Full_m)):
    if 0.3 <= Full_m[i] <= 1.2 and 0.01 <= Full_Prot[i] <= 12:
        m.append(Full_m[i])
        Prot.append(Full_Prot[i])
        
#Resample data limiting mass and rotation period again        
n_sample=1 #Sample number
step=1 #One data point at a time
c=0 #Count

#Empty resampled data arrays with sample number size
Prot_resampled=np.zeros(n_sample) 
M_resampled=np.zeros(n_sample)

#Choosing the resampled data points that fall inside the mass and rotation period boundaries we set
while c < n_sample:
    values = np.vstack([Prot, m]) # create 2D array that contains the properties you want to resample
    kde_ProtMass = stats.gaussian_kde(values) # calculate 2D KDE of Prot-vs-Mstar
    
    #Create a resampled data point
    re_Prot = kde_ProtMass.resample(step)[0,:]
    re_M = kde_ProtMass.resample(step)[1,:]
    
    if 0.3 <= re_M <= 1.2 and 0.01 <= re_Prot <= 12:
        M_resampled[c]=re_M
        Prot_resampled[c]=re_Prot
        c+=1

stars=5
Prots=[]

for n_star in range(stars):
    #Mass of the object we are studying using star data
    solar_mass=star_data.mass[n_star]
    used_mass=find_nearest(MASSES, solar_mass)
    M0 = np.ones((Prot_resampled.shape))*used_mass
    
    
    # t0: time [Myr] at which the model should start 
    # tdisc: disk-locking time [Myr]
    Prot_evol, age_zero, Prot_interp, spl_Prot = spin_down_evol(Prot_init=Prot_resampled, 
                                                            Mstar_init=M0, 
                                                            t0=1., tdisc=13.)
    #Select the rotation period for the age of the star
    age=star_data.age[n_star]*(10**3) #Age of star (it is given in Gyr and we pass it to Myr)
    print(age)
    Calculated_Prot=spl_Prot(age)
    Prots.append(Calculated_Prot)
    
    #Plot Age vs Rotation period
    fig = plt.figure(figsize=(15, 10))
    fig.add_axes([0.1, 0.1, 0.8, 0.8]) #real plot
        
    plt.xlim([10,6*10**4])
    plt.ylim([ 0.01, Calculated_Prot+50])
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('P [days]', fontsize=15)
    plt.scatter(age,Calculated_Prot, c='black', marker='$\odot$', s=500)
    plt.title(f'Evolution of the rotation period of {solar_mass}$M_\odot$ star', fontsize= 25)
    for el in range(len(Prot_evol)):
        plt.plot( Prot_evol[el,1], Prot_evol[el,2], c='orange')
        
        
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)
    #plt.text(20, 80, r'Initial Mass={}'.format(1.0)+'$M_{\odot}$',  fontsize=20)
    
    plt.savefig('/home/farah/Documents/Project/Data/Star_Evolution_Prot.png')

print(Prots)
