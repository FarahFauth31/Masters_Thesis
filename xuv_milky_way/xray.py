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
from xuv_milky_way import common

#/File/path/where/you/have/the/MIST_tables/
sys.path.append( '/home/farah/Documents/Project/Data/MIST_tables/' ) #You can download MIST tables from https://github.com/cgarraffo/Spin-down-model/tree/master/Mist_data


def constants():
    """
    This function calculates the constants needed for Xray calculations using the Wright 2011 coronal-Rossby relation.

    return: saturated Lx_Lbol value, saturated Rossby value, beta (slope) and C (constant).
    """
    
    # Wright et al. 2018 equation: (Lx/Lbol)sat = C*(R_{o,sat}^beta) --- C = (Lx/Lbol)sat / (R_{o,sat}^beta)
    
    Lx_Lbol_sat = 10**(-3.05)
    Ro_sat = 0.14
    beta = -2.3
    
    C = Lx_Lbol_sat/(Ro_sat**beta)
    
    return Lx_Lbol_sat, Ro_sat, beta, C


def empirical_tau(mass):
    """
    This function calculates the convective turnover time from an empirical prescription in Drake et al. (in prep).

    mass: stellar mass.
    return: calculated convective turnover time.
    """
    
    t = 10**((-1.101593022796932*mass) + 2.240877938410134)
    return t


#Create mass array with all the MIST masses we have
MASSES = np.arange(0.1,2.05, 0.05)


#Function that loads all the MIST tables
def load_mist_tables(Mstar=1., filepath='/home/farah/Documents/Project/Data/MIST_tables/'):
        """
        This function loads in the MIST tables.

        Mstar: Stellar masses in units of solar masses.
        filepath: Path where the MIST tables are stored (string).
        return: arrays of the MIST tables that contain AGE and TAU (convective turnover time).
        """
        from MIST_tables import read_mist_models
        import astropy.units as u

        # print(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)))

        eep = read_mist_models.EEP(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)), verbose=False)
        AGE_mist = (eep.eeps['star_age']*u.yr).to(u.Gyr) # stellar age in Myears
        TAU_mist = (eep.eeps['conv_env_turnover_time_g']*u.s).to(u.d) # convective turnover time in days
        LBOL_mist = 10**(eep.eeps['log_L']) # Bolometric luminosity in units of solar luminosity
        # log_RADIUS_mist = eep.eeps['log_R']
        # RADIUS_mist = 10**log_RADIUS_mist

        # return AGE_mist, TAU_mist, LBOL_mist, RADIUS_mist
        return AGE_mist, TAU_mist, LBOL_mist

def age_MS(mass):
    """
    This function assigns the main squence age to each star depending on its mass.

    mass: stellar mass in units of solar mass.
    return: main sequence age value.
    """
    t_sol=10 #Main sequence lifetime of the Sun (in Gyrs)
    m_sol=1 #Mass of Sun in solar masses
    
    if mass<=0.95:
        
        age_main_sequence=t_sol*(m_sol/mass**(2.5))
        
    if mass>0.95 and mass<=1:

        age_main_sequence=9.5 # in Gyrs
    
    if mass>1 and mass<=1.05:

        age_main_sequence=7 # in Gyrs
        
    if mass>1.05 and mass<=1.1:

        age_main_sequence=5.5 # in Gyrs
        
    if mass>1.1 and mass<=1.15:

        age_main_sequence=4.5 # in Gyrs
        
    if mass>1.15 and mass<=1.2:

        age_main_sequence=4 # in Gyrs
        
    if mass>1.2 and mass<=1.25:

        age_main_sequence=3.5 # in Gyrs    
        
    return age_main_sequence

    
def open_all_mist_files():
    """
    This function opens all MIST tables and puts them in a main array for easy access and time efficiency.

    return: main array for all MIST ages, main array for all MIST convective turnover times, main array for all MIST bolometric luminosities.
    """

    #Open grid of the masses
    main_array_ages=[]
    main_array_tau=[]
    main_array_lbol=[]
    
    arrayname1a, arrayname1t, arrayname1l = load_mist_tables(Mstar=MASSES[0])
    arrayname2a, arrayname2t, arrayname2l = load_mist_tables(Mstar=MASSES[1])
    arrayname3a, arrayname3t, arrayname3l = load_mist_tables(Mstar=MASSES[2])
    arrayname4a, arrayname4t, arrayname4l = load_mist_tables(Mstar=MASSES[3])
    arrayname5a, arrayname5t, arrayname5l = load_mist_tables(Mstar=MASSES[4])
    arrayname6a, arrayname6t, arrayname6l = load_mist_tables(Mstar=MASSES[5])
    arrayname7a, arrayname7t, arrayname7l = load_mist_tables(Mstar=MASSES[6])
    arrayname8a, arrayname8t, arrayname8l = load_mist_tables(Mstar=MASSES[7])
    arrayname9a, arrayname9t, arrayname9l = load_mist_tables(Mstar=MASSES[8])
    arrayname10a, arrayname10t, arrayname10l = load_mist_tables(Mstar=MASSES[9])
    arrayname11a, arrayname11t, arrayname11l = load_mist_tables(Mstar=MASSES[10])
    arrayname12a, arrayname12t, arrayname12l = load_mist_tables(Mstar=MASSES[11])
    arrayname13a, arrayname13t, arrayname13l = load_mist_tables(Mstar=MASSES[12])
    arrayname14a, arrayname14t, arrayname14l = load_mist_tables(Mstar=MASSES[13])
    arrayname15a, arrayname15t, arrayname15l = load_mist_tables(Mstar=MASSES[14])
    arrayname16a, arrayname16t, arrayname16l = load_mist_tables(Mstar=MASSES[15])
    arrayname17a, arrayname17t, arrayname17l = load_mist_tables(Mstar=MASSES[16])
    arrayname18a, arrayname18t, arrayname18l = load_mist_tables(Mstar=MASSES[17])
    arrayname19a, arrayname19t, arrayname19l = load_mist_tables(Mstar=MASSES[18])
    arrayname20a, arrayname20t, arrayname20l = load_mist_tables(Mstar=MASSES[19])
    arrayname21a, arrayname21t, arrayname21l = load_mist_tables(Mstar=MASSES[20])
    arrayname22a, arrayname22t, arrayname22l = load_mist_tables(Mstar=MASSES[21])
    arrayname23a, arrayname23t, arrayname23l = load_mist_tables(Mstar=MASSES[22])
    arrayname24a, arrayname24t, arrayname24l = load_mist_tables(Mstar=MASSES[23])

    
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
    
    main_array_lbol.append(arrayname1l)
    main_array_lbol.append(arrayname2l)
    main_array_lbol.append(arrayname3l)
    main_array_lbol.append(arrayname4l)
    main_array_lbol.append(arrayname5l)
    main_array_lbol.append(arrayname6l)
    main_array_lbol.append(arrayname7l)
    main_array_lbol.append(arrayname8l)
    main_array_lbol.append(arrayname9l)
    main_array_lbol.append(arrayname10l)
    main_array_lbol.append(arrayname11l)
    main_array_lbol.append(arrayname12l)
    main_array_lbol.append(arrayname13l)
    main_array_lbol.append(arrayname14l)
    main_array_lbol.append(arrayname15l)
    main_array_lbol.append(arrayname16l)
    main_array_lbol.append(arrayname17l)
    main_array_lbol.append(arrayname18l)
    main_array_lbol.append(arrayname19l)
    main_array_lbol.append(arrayname20l)
    main_array_lbol.append(arrayname21l)
    main_array_lbol.append(arrayname22l)
    main_array_lbol.append(arrayname23l)
    main_array_lbol.append(arrayname24l)

    
    return main_array_ages, main_array_tau, main_array_lbol


def calculate_xray(RA_steps, DEC_steps, GUMS_file_path):
    """
    This function calculates the Xray emission of stars in document from basic stellar parameters using an empirical description of the convective turnover time and the bolometric luminosity from MIST tables.

    RA_steps: Steps of RA for file names.
    DEC_steps: Steps of DEC for file names.
    GUMS_file_directory: Path of file with GUMS data or stellar data.
    return: saves original csv file with an extra column for calculated Xray emissions.
    """

    main_array_ages, main_array_tau, main_array_lbol = open_all_mist_files()
    
    Lx_Lbol_sat, Ro_sat, beta, C = constants()
    
    LX=[]
    
    # Loop that looks into each file
    for k in RA_steps:
                            
        for b in DEC_steps:

            a=1
            
            # Loop that looks into the parts of files in each section of sky
            while a<500:
                
                # Trying to find the file
                found = common.find(f'RA_{k}_{k+4}_DEC_{b}_{b+4}_target.csv_Part{a}',GUMS_file_path)

                # If the file is not found the loop looks into next section of sky
                if found == None:
                    
                    a=500
                
                # If file is found we calculate the Xray emission of stars for each document
                else:
                                                    
                    print(f'RA_{k}_{k+4}_DEC_{b}_{b+4}_target.csv_Part{a}')
                    
                    #Open csv file containing star data
                    prot_data = pd.read_csv(found) #Read csv file
                    
                    'BODY OF THE CODE'

                    #Create empty lists that we will use to create a new column
                    # PROT=[]
                    LX=[]

                    n_star=len(prot_data['ra']) #Number of stars in file we want to evaluate
    
                    #Do this loop for each star in file
                    for i in range(n_star):
                        
                        #Mass of star
                        mass = prot_data.mass[i]
                        
                        #Age of star
                        age = prot_data.uniform_ages[i] #Unit [Gyears]
                        
                        #Rotation period of star
                        Prot = prot_data.Prot[i] #Units [days]

                        # Calculate age MS of star                    
                        age_main_sequence=age_MS(mass)
                                            
                        #If stars age is bigger than Main Sequence turn off then Lx is nan
                        if age >= age_main_sequence:
                                
                            Lx = np.nan
                            Prot = np.nan
                        
                        #Else calculate Lx
                        else: 
                            
                            #Find indexes for the mass documents with values nearest to the star's mass
                            nearest_mass, nearest_mass2 = common.find_2_nearest(MASSES, mass)
                            index_nearest_mass = int(np.where(MASSES==nearest_mass)[0])
                            index_nearest_mass2 = int(np.where(MASSES==nearest_mass2)[0])
                            
                            #Load age data for that mass
                            AGE_mist = main_array_ages[index_nearest_mass]
                                                
                            #Find index for the age of the star in MIST table
                            nearest_age = common.find_nearest(AGE_mist, age)
                            index = int(np.where(AGE_mist.value == nearest_age)[0])
                            
                            #Calculate convective turnover time
                            tau = empirical_tau(mass)
                            
                            #Rossby number
                            Ro = Prot/tau #Unitless
                            
                            #Load Lbol data from the two mass documents found before
                            LBOL_mist = main_array_lbol[index_nearest_mass]
                            LBOL_mist2 = main_array_lbol[index_nearest_mass2]
                            
                            #Get Lbol of star from the two mass documents at the age required
                            Lbol1 = LBOL_mist[index] #Units of solar luminosity
                            Lbol2 = LBOL_mist2[index] #Units of solar luminosity
                            
                            #Create a list with the two nearest masses and two calculated Lbols for that star
                            two_masses=[nearest_mass,nearest_mass2]
                            two_lbol=[Lbol1,Lbol2]
                            
                            #Calculate Lbol by interpolating between mass documents
                            Lbol = common.interpolation(two_masses,two_lbol,mass)
        
                            #If Rossby number is smaller than the saturated Rossby number then we use the saturated Lx/Lbol value to calculate Lx              
                            if Ro < Ro_sat:
                                
                                Lx_Lbol = Lx_Lbol_sat
                                
                                Lx = Lbol*Lx_Lbol # units of solar luminosity [L_\odot]
                            
                            #If Ro is bigger than Ro_sat then Lx/Lbol follows the relationship from Wright et al. 2011
                            else:
        
                                Lx_Lbol = C*(Ro**beta)
                                
                                Lx = Lbol*Lx_Lbol # units of solar luminosity [L_\odot]
                                                                                    
                            
                        LX.append(Lx)
                        # PROT.append(Prot)
                        # LXLBOL.append(Lx_Lbol)
                        # RO.append(Ro)
                        # TAU.append(tau)
                        
                        if Lx_Lbol > Lx_Lbol_sat:
                            print('yes')
                        

                    #dictionary = {'Lx': LX} 
                    #dataframe = pd.DataFrame(dictionary)
                    #prot_data['Lx'] = dataframe
                    #prot_data.to_csv(found, index=False)
                    print(LX)
                                    
                    a+=1