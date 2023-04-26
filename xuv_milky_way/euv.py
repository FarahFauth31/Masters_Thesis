# Calculate surface X-ray flux and corresponding EUV spectrum

'Calculate surface X-ray flux and corresponding EUV spectrum'

import numpy as np
import pandas as pd
import pickle
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
from scipy import stats
import time
import os
from scipy.io.idl import readsav
import sys
from scipy.integrate import simpson
from numpy import trapz
import astropy.units as u
from xuv_milky_way import common


#/File/path/where/you/have/the/MIST_tables/
sys.path.append( '/home/farah/Documents/Project/Data/MIST_tables/' ) #You can download MIST tables from https://github.com/cgarraffo/Spin-down-model/tree/master/Mist_data


#Create mass array with all the MIST masses we have
MASSES = np.arange(0.1,2.05, 0.05)


#Function that loads all the MIST tables
def load_mist_tables(Mstar=1., filepath='/home/farah/Documents/Project/Data/MIST_tables/'):
        """
        This function loads in the MIST tables.

        Mstar: Stellar masses in units of solar masses.
        filepath: Path where the MIST tables are stored (string).
        return: arrays of the MIST tables that contain AGE, TAU (convective turnover time), LBOL, RADIUS (in that order).
        """
        
        from MIST_tables import read_mist_models
        import astropy.units as u

        print(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)))

        eep = read_mist_models.EEP(filepath+'00{:03d}M.track.eep'.format(int(Mstar*100)), verbose=False)
        AGE_mist = (eep.eeps['star_age']*u.yr).to(u.Myr) # stellar age in Myears
        TAU_mist = (eep.eeps['conv_env_turnover_time_g']*u.s).to(u.d) # convective turnover time in days
        LBOL_mist = 10**(eep.eeps['log_L']) # Bolometric luminosity in units of solar luminosity
        log_RADIUS_mist = eep.eeps['log_R']
        RADIUS_mist = 10**log_RADIUS_mist

        return AGE_mist, TAU_mist, LBOL_mist, RADIUS_mist


def wood_spectra(IDL_file_path):
    """
    This function uploads all the information contained in the IDL file.

    IDL_file_path: file path of IDL file.
    return: returns the arrays containing data for wavelength (array[1100]), wavelength vs log surface flux matrix (array[1100, 36]) and vector of logarithmic surface fluxes (array[36]).
    """ 
    
    # Equation used F = L/(4*pi*R^2) --- Fx/Lx = Fbol/Lbol --- Fx = Lx * (Fbol/Lbol)
    
    #Open data files (IDL and .dat)
    s = readsav(IDL_file_path)
    data = np.loadtxt('/home/farah/Documents/Project/Data/EUV/log_FX.dat')
    #log_T = [4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5]
    
    #Data from IDL file
    wavelength = s['sp_wvl'] #Wavelength, Array[1100]
    matrix = s['woodspec'] #[wavelength,log surface flux] matrix, Array[1100, 36]
    log_F = s['wfxi'] #Vector of logarithmic surface fluxes, Array[36]
    
    return wavelength, matrix, log_F


def xray_wav_interval(wavelength):
    """
    This function gives you the array indices for the wavelength interval of the Xray regime.

    wavelength: array with wavelength values.
    return: returns the index of the wavelength value of the start and end of the Xray regime in the given wavelength array.
    """    

    # The y values
    # Fnd index of the wavelength array where value equals the start and end of the xray regime
    index1=int(np.where(wavelength==1)[0])
    index2=int(np.where(wavelength==30)[0])
    
    return index1, index2


def which_spectra(GUMS_file_path, IDL_file_path, n_star):
    """
    This function finds the index of the spectrum that corresponds to a star and the Lx of the star.

    GUMS_file_path: file path of the star data (string).
    IDL_file_path: file path of the IDL file containing the spectra data (string).
    n_star: number of stars in the document from which you want to calculate the Lx and spectrum.
    return: returns the index of the best fitting spectrum and the Lx value
    """

    #Open csv file containing star data
    star_data = pd.read_csv(GUMS_file_path) #Read csv file
    
    R_star = star_data.radius #Absolute bolometric magnitude data
    
    wavelength, matrix, log_F = wood_spectra(IDL_file_path)
    
    index1, index2 = xray_wav_interval(wavelength)

    IDX=[]
    LX=[]
    
    # print(index1, index2, wavelength[index1], wavelength[index2], wavelength[index1:index2+1])
    
    for i in range(n_star):
        
        if np.isnan(R_star[i]):
            
            #Mass of star
            mass = star_data.mass[i]
            
            nearest_mass = common.find_nearest(MASSES, mass)
            
            #Load convection turnover time data for that mass
            AGE_mist, TAU_mist, LBOL_mist, RADIUS_mist = load_mist_tables(Mstar=nearest_mass)
            
            #Age of star
            age = star_data.age[i] #Unit [Myears]
            nearest_age = common.find_nearest(AGE_mist, age)
            index = int(np.where(AGE_mist.value == nearest_age)[0])
            
            #Radius of star
            R_solar = RADIUS_mist[index] #Units [R_\odot]
                    
        else:
            
            R_solar = R_star[i] # [Solar radius, R_\odot]
            
        R = R_solar*(6.957*(10**8)) # [m]
        
        Lx_solar = star_data.Lx[i] #Solar luminosity units [L_\odot]
    
        Lx = Lx_solar*(3.826*(10**33)) # [erg s^(-1)]
            
        Fx_solar = Lx/(4*np.pi*(R**2)) #Units of solar luminosity divided by meters squared [erg s^(-1) m^(-2)]
    
        Fx = Fx_solar/10000 # [erg s^(-1) cm^(-2)]
    
        log_Fx = np.log10(Fx)
        
        nearest_log_Fx = common.find_nearest(log_F,log_Fx)
        idx = int(np.where(log_F==nearest_log_Fx)[0])

        IDX.append(idx)
        LX.append(Lx)    
    
    return matrix, wavelength, IDX, index1, index2, LX
        

def norm_spectra(matrix, wavelength, IDX, index1, index2, LX, n_star):
    """
    This function normalizes the spectrum of the star.

    matrix: matrix from IDL file that contains data of all spectra depending on Fx.
    wavelength: array with all wavelength values.
    IDX: array with index of best fitting spectrum.
    index1: index for the start of Xray regime.
    index2: index for end of Xray regime.
    LX: array with Xray luminosity of stars.
    n_star: number of stars in the document from which you want to calculate the normalized spectrum.
    return: returns Xray wavelength interval, the new points of the spectrum after normalization and the old points of the spectrum before normalization.
    """

    Y=[]
    NEW_Y=[]

    for i in range(n_star):

        # Compute the area using the composite Simpson's rule.
        y=matrix[IDX[i]][index1:index2+1]
        xray_w=wavelength[index1:index2+1]
        area = simpson(y, dx=0.05)

        d=LX[i]/area

        new_y=y*d

        Y.append(y)
        NEW_Y.append(new_y)
    
    return xray_w, NEW_Y, Y


def plot_spectra(xray_w, NEW_Y, Y, n_star):
    """
    This function plots the original and normalized spectra.

    xray_w: Xray wavelength interval.
    NEW_Y: array with the new points of the spectrum after normalization of each star.
    Y: array with the old points of the spectrum before normalization of each star.
    n_star: number of stars in the document from which you want to plot the original and normalized spectra.
    return: plots.
    """

    for i in range(n_star):
        #Original EUV spectrum
        fig = plt.figure(figsize=(10, 7))
        plt.xlim([0.3,30])
        plt.xlabel('Wavelength [$\AA$]', fontsize=15)
        plt.ylabel('$L_{x,\u03BB}$ [erg $s^{-1}$]', fontsize=15)
        plt.title(f'Original EUV Spectra - Star {i}', fontsize= 20)
        plt.plot(xray_w,Y[i])
        ax = plt.gca()
        plt.tick_params(labelsize=14)
        plt.tick_params(labelsize=14)
        plt.tick_params(labelsize=14)
        plt.tick_params(axis='y', which='both')

        #Normalized EUV spectrum
        fig = plt.figure(figsize=(10, 7))
        plt.xlim([0.3,30])
        plt.xlabel('Wavelength [$\AA$]', fontsize=15)
        plt.ylabel('$L_{x,\u03BB}$ [erg $s^{-1}$]', fontsize=15)
        plt.title(f'Normalized EUV Spectra - Star {i}', fontsize= 20)
        plt.plot(xray_w,NEW_Y[i])
        ax = plt.gca()
        plt.tick_params(labelsize=14)
        plt.tick_params(labelsize=14)
        plt.tick_params(labelsize=14)
        plt.tick_params(axis='y', which='both')