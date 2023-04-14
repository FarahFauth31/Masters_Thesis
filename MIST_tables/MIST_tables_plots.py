import numpy as np
import astropy.units as u
import astropy.constants as const
import sys
import matplotlib.pyplot as plt

'MIST table plots'

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
        AGE_mist = (eep.eeps['star_age']*u.yr).to(u.Myr) # stellar age in Myears
        TAU_mist = (eep.eeps['conv_env_turnover_time_g']*u.s).to(u.d) # convective turnover time in days
        MOI_mist = eep.eeps['moment_of_inertia']*u.g*u.cm**2. # moment of inertia in cgs
        MASS_mist = eep.eeps['star_mass']
        log_RADIUS_mist = eep.eeps['log_R']
        RADIUS_mist = 10**log_RADIUS_mist

        return AGE_mist, TAU_mist, MOI_mist, MASS_mist, RADIUS_mist


mass = np.arange(0.1,1.3,0.05) #0.05 steps
name_mass = [0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15, 1.2, 1.25] #exact masses

n_steps=len(mass)

for i in range(n_steps):
    AGE_mist, TAU_mist, MOI_mist, MASS_mist, RADIUS_mist=load_mist_tables(Mstar=find_nearest(MASSES, mass[i])) #Look at how many time steps the MIST table has
    
    #Tau plot
    fig = plt.figure(figsize=(15, 10))
    plt.xlim([10,3*10**4])
    # plt.ylim([ 0, 500])
    plt.yscale('log')
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('$\u03C4$ [days]', fontsize=15)
    plt.xscale('log')
    plt.title(f'Evolution of $\u03C4$ of {name_mass[i]}$M_\odot$ star', fontsize= 20)
    plt.scatter(AGE_mist, TAU_mist, c='blue')
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)
    
    #MoI plot
    fig = plt.figure(figsize=(15, 10))
    plt.xlim([10,3*10**4])
    plt.ylim([ 0, 2])
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('Moment of Inertia [$x 10^{54} g*s^2$]', fontsize=15)
    plt.xscale('log')
    plt.title(f'Evolution of MoI of {name_mass[i]}$M_\odot$ star', fontsize= 20)
    plt.scatter(AGE_mist, MOI_mist/(10**(54)), c='red')
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)
    
    #Mass plot
    fig = plt.figure(figsize=(15, 10))
    plt.xlim([10,3*10**4])
    #plt.ylim([ 0, 10000000000])
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('Mass [$M_\odot$]', fontsize=15)
    plt.xscale('log')
    plt.title(f'Evolution of mass of {name_mass[i]}$M_\odot$ star', fontsize= 20)
    plt.scatter(AGE_mist, MASS_mist, c='green')
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)
    
    #Radius plot
    fig = plt.figure(figsize=(15, 10))
    plt.xlim([10,3*10**4])
    plt.ylim([ 0, 2])
    plt.xlabel('Age [Myrs]', fontsize=15)
    plt.ylabel('Radius [$R_\odot$]', fontsize=15)
    plt.xscale('log')
    plt.title(f'Evolution of radius of {name_mass[i]}$M_\odot$ star', fontsize= 20)
    plt.scatter(AGE_mist, RADIUS_mist, c='yellow')
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(labelsize=14)
    plt.tick_params(axis='y', which='both')
    plt.tick_params(bottom=True, left=True, right=True)
    plt.tick_params(which='both',labelright=True)