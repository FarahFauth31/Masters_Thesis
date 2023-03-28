# Calculate Xray emission for all documents

'Calculate Xray emission for all documents'

import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
import sys
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy import stats
import time
import os

#/File/path/where/you/have/the/MIST_tables/
sys.path.append( '/home/farah/Documents/Project/Data/MIST_tables/' ) #You can download MIST tables from https://github.com/cgarraffo/Spin-down-model/tree/master/Mist_data


def constants():
    
    # 'Find constant C'
    
    # Wright et al. 2018 equation: (Lx/Lbol)sat = C*(R_{o,sat}^beta) --- C = (Lx/Lbol)sat / (R_{o,sat}^beta)
    
    Lx_Lbol_sat = 10**(-3.05)
    Ro_sat = 0.14
    beta = -2.3
    
    C = Lx_Lbol_sat/(Ro_sat**beta)
    
    return Lx_Lbol_sat, Ro_sat, beta, C


def empirical_tau(mass):
    t = 10**((-1.101593022796932*mass) + 2.240877938410134)
    return t


# 'Calculate Ro'

# Ro = Prot/tau

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

# Function that finds a file and returns its path
def find(name, path):
    for root, dirs, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)


#Function that loads all the MIST tables
def load_mist_tables(Mstar=1., filepath='/home/farah/Documents/Project/Data/MIST_tables/'):
        """
        Load in the MIST tables.
        Mstar: Stellar masses in units of solar masses
        filepath: Path where the MIST tables are stored
        """
        import read_mist_models
        import astropy.units as u

        # print(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)))

        eep = read_mist_models.EEP(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)), verbose=False)
        AGE_mist = (eep.eeps['star_age']*u.yr).to(u.Myr) # stellar age in Myears
        TAU_mist = (eep.eeps['conv_env_turnover_time_g']*u.s).to(u.d) # convective turnover time in days
        # LBOL_mist = 10**(eep.eeps['log_L']) # Bolometric luminosity in units of solar luminosity
        # log_RADIUS_mist = eep.eeps['log_R']
        # RADIUS_mist = 10**log_RADIUS_mist

        # return AGE_mist, TAU_mist, LBOL_mist, RADIUS_mist
        return AGE_mist, TAU_mist
    
def open_all_mist_files():
    #Open grid of the masses
    main_array_ages=[]
    main_array_tau=[]
    
    arrayname1a, arrayname1t = load_mist_tables(Mstar=MASSES[0])
    arrayname2a, arrayname2t = load_mist_tables(Mstar=MASSES[1])
    arrayname3a, arrayname3t = load_mist_tables(Mstar=MASSES[2])
    arrayname4a, arrayname4t = load_mist_tables(Mstar=MASSES[3])
    arrayname5a, arrayname5t = load_mist_tables(Mstar=MASSES[4])
    arrayname6a, arrayname6t = load_mist_tables(Mstar=MASSES[5])
    arrayname7a, arrayname7t = load_mist_tables(Mstar=MASSES[6])
    arrayname8a, arrayname8t = load_mist_tables(Mstar=MASSES[7])
    arrayname9a, arrayname9t = load_mist_tables(Mstar=MASSES[8])
    arrayname10a, arrayname10t = load_mist_tables(Mstar=MASSES[9])
    arrayname11a, arrayname11t = load_mist_tables(Mstar=MASSES[10])
    arrayname12a, arrayname12t = load_mist_tables(Mstar=MASSES[11])
    arrayname13a, arrayname13t = load_mist_tables(Mstar=MASSES[12])
    arrayname14a, arrayname14t = load_mist_tables(Mstar=MASSES[13])
    arrayname15a, arrayname15t = load_mist_tables(Mstar=MASSES[14])
    arrayname16a, arrayname16t = load_mist_tables(Mstar=MASSES[15])
    arrayname17a, arrayname17t = load_mist_tables(Mstar=MASSES[16])
    arrayname18a, arrayname18t = load_mist_tables(Mstar=MASSES[17])
    arrayname19a, arrayname19t = load_mist_tables(Mstar=MASSES[18])
    arrayname20a, arrayname20t = load_mist_tables(Mstar=MASSES[19])
    arrayname21a, arrayname21t = load_mist_tables(Mstar=MASSES[20])
    arrayname22a, arrayname22t = load_mist_tables(Mstar=MASSES[21])
    arrayname23a, arrayname23t = load_mist_tables(Mstar=MASSES[22])
    arrayname24a, arrayname24t = load_mist_tables(Mstar=MASSES[23])

    
    main_array_ages.append(arrayname1a)
    main_array_ages.append(arrayname2a)
    main_array_ages.append(arrayname3a)
    main_array_ages.append(arrayname4a)
    main_array_ages.append(arrayname5a)
    main_array_ages.append(arrayname6a)
    main_array_ages.append(arrayname7a)
    main_array_ages.append(arrayname8a)
    main_array_ages.append(arrayname9a)
    main_array_ages.append(arrayname10a)
    main_array_ages.append(arrayname11a)
    main_array_ages.append(arrayname12a)
    main_array_ages.append(arrayname13a)
    main_array_ages.append(arrayname14a)
    main_array_ages.append(arrayname15a)
    main_array_ages.append(arrayname16a)
    main_array_ages.append(arrayname17a)
    main_array_ages.append(arrayname18a)
    main_array_ages.append(arrayname19a)
    main_array_ages.append(arrayname20a)
    main_array_ages.append(arrayname21a)
    main_array_ages.append(arrayname22a)
    main_array_ages.append(arrayname23a)
    main_array_ages.append(arrayname24a)
    
    main_array_tau.append(arrayname1t)
    main_array_tau.append(arrayname2t)
    main_array_tau.append(arrayname3t)
    main_array_tau.append(arrayname4t)
    main_array_tau.append(arrayname5t)
    main_array_tau.append(arrayname6t)
    main_array_tau.append(arrayname7t)
    main_array_tau.append(arrayname8t)
    main_array_tau.append(arrayname9t)
    main_array_tau.append(arrayname10t)
    main_array_tau.append(arrayname11t)
    main_array_tau.append(arrayname12t)
    main_array_tau.append(arrayname13t)
    main_array_tau.append(arrayname14t)
    main_array_tau.append(arrayname15t)
    main_array_tau.append(arrayname16t)
    main_array_tau.append(arrayname17t)
    main_array_tau.append(arrayname18t)
    main_array_tau.append(arrayname19t)
    main_array_tau.append(arrayname20t)
    main_array_tau.append(arrayname21t)
    main_array_tau.append(arrayname22t)
    main_array_tau.append(arrayname23t)
    main_array_tau.append(arrayname24t)

    
    return main_array_ages, main_array_tau


def calculate_xray(RA_steps, DEC_steps, GUMS_file_path):

    main_array_ages, main_array_tau = open_all_mist_files()
    
    Lx_Lbol_sat, Ro_sat, beta, C = constants()
    
    LX=[]
    
    # Loop that looks into each file
    for k in RA_steps:
                            
        for b in DEC_steps:
    
            a=1
            
            while a<500:
                            
                found = find(f'RA_{k}_{k+4}_DEC_{b}_{b+4}_target.csv_Part{a}',GUMS_file_path)
                
                if found == None:
                    
                    a=500
                
                else:
                    
                    print(f'RA_{k}_{k+4}_DEC_{b}_{b+4}_target.csv_Part{a}')
                    
                    #Open csv file containing star data
                    prot_data = pd.read_csv(found) #Read csv file
    
                    'BODY OF THE CODE'
    
                    #Create empty lists that we will use to create a new file with the calculated data
                    # M=[]
                    # AGE=[]
                    PROT=[]
    
                    #Do this loop for each star in file
                    n_star=len(prot_data['ra']) #Number of stars in file we want to evaluate
                    
                    for i in range(n_star):
                        
                        #Mass of star
                        mass = prot_data.mass[i]
                        
                        if mass >= 0.4:
                            
                            nearest_mass = find_nearest(MASSES, mass)
                            index_nearest_mass = int(np.where(MASSES==nearest_mass)[0])
                        
                            #Load convection turnover time data for that mass
                            AGE_mist = main_array_ages[index_nearest_mass]
                            TAU_mist = main_array_tau[index_nearest_mass]
                                
                            #Age of star
                            age = prot_data.age[i] #Unit [Myears]
                            nearest_age = find_nearest(AGE_mist, age)
                            index = int(np.where(AGE_mist.value == nearest_age)[0])
                            
                            #Tau at specific age
                            tau = TAU_mist[index].value #Units [days]
                            
                            
                        else:
                            
                            tau = empirical_tau(mass)
    
    
                            
                        #Rotation period of star
                        Prot = prot_data.Prot[i] #Units [days]
                        
                        Ro = Prot/tau #Unitless
                        
                        Mbol = prot_data.mbol #Absolute bolometric magnitude data
                        
                        Lbol = 10**(0.4*(4.85-Mbol[i])) #Units of solar luminosity
                                            
                        if Ro < Ro_sat:
                            
                            Lx_Lbol = Lx_Lbol_sat
                            
                            Lx = Lbol*Lx_Lbol # units of solar luminosity [L_\odot]
                            
                        else:
    
                            Lx_Lbol = C*(Ro**beta)
                            
                            Lx = Lbol*Lx_Lbol # units of solar luminosity [L_\odot]
                            
                        LX.append(Lx)
                        
                        if Lx_Lbol > Lx_Lbol_sat:
                            print('yes')
                        
    
                    dictionary = {'Lx': LX}  
                    dataframe = pd.DataFrame(dictionary) 
                    prot_data['Lx'] = dataframe
                    prot_data.to_csv(found, index=False)
                                    
                    a+=1