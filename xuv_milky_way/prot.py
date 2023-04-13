#Calculate Prot of GUMS stars for all documents - More efficient

import numpy as np
import pandas as pd
import pickle
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
from scipy import stats
import time
import os
from xuv_milky_way import common

'Calculate Prot of GUMS stars for all documents - More efficient'
        
def kde(file_directory):
    """
    This function creates a 2D Gaussian kde distribution from cluster data.

    file_directory: Path of the cluster data document.
    return: Prot-Mass 2D Gaussian distribution.
    """
    
    #Open csv file containing initial rotation period distribution data
    data = pd.read_csv(file_directory) #Read hPer csv file
            
    df = pd.DataFrame(data)
            
    #Create new lists. The mass data is limited between 0.1 < m < 1.25 solar masses and the rotation period data is limited between 0.01 < Prot < 12 days.
    red_data = df[(df['Mass']>=0.1) & (df['Mass']<=1.25) & (df['Per']>=0.1) & (df['Per']<=12)]
            
    #Create an initial rotation period distribution from hPer data
    Full_m=red_data.Mass #Select the data set you want, in this case the mass (y data)
    Full_Prot=red_data.Per #Select the data set you want, in this case the rotation periods
    values = np.vstack([Full_Prot, Full_m]) # create 2D array that contains the properties you want to resample
    kde_ProtMass = stats.gaussian_kde(values) # calculate 2D KDE of Prot-vs-Mstar
    
    return kde_ProtMass
    
        
def initial_Prot(kde_ProtMass):
    """
    This function resamples a initial rotation period value from a 2D Gaussian kde distribution.

    kde_ProtMass: Prot-Mass 2D Gaussian distribution.
    return: resampled rotation period value.
    """ 

    #Resample data limiting mass and rotation period again        
    step=1 #One data point at a time
    c=0 #Count

    #Choosing the resampled data points that fall inside the mass and rotation period boundaries we set
    while c < 1:
        #Create a resampled data point
        re_Prot = kde_ProtMass.resample(step)[0,:]
        re_M = kde_ProtMass.resample(step)[1,:]
        
        if 0.1 <= re_M <= 1.25 and 0.1 <= re_Prot <= 12:
            c+=1
    
    return re_Prot

def open_mass_files():
    """
    This function opens all mass files from the Rotation Period Evolution Grid.

    return: main array with all loaded mass documents.
    """
        
    #Open grid of the masses
    main_array=[]
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.1_evolution','rb') as af: arrayname1 = pickle.load(af)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.15_evolution','rb') as ag: arrayname2 = pickle.load(ag)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.2_evolution','rb') as ah: arrayname3 = pickle.load(ah)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.25_evolution','rb') as ai: arrayname4 = pickle.load(ai)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.3_evolution','rb') as aj: arrayname5 = pickle.load(aj)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.35_evolution','rb') as ak: arrayname6 = pickle.load(ak)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.4_evolution','rb') as al: arrayname7 = pickle.load(al)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.45_evolution','rb') as am: arrayname8 = pickle.load(am)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.5_evolution','rb') as an: arrayname9 = pickle.load(an)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.55_evolution','rb') as ao: arrayname10 = pickle.load(ao)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.6_evolution','rb') as ap: arrayname11 = pickle.load(ap)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.65_evolution','rb') as aq: arrayname12 = pickle.load(aq)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.7_evolution','rb') as ar: arrayname13 = pickle.load(ar)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.75_evolution','rb') as ass: arrayname14 = pickle.load(ass)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.8_evolution','rb') as at: arrayname15 = pickle.load(at)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.85_evolution','rb') as au: arrayname16 = pickle.load(au)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.9_evolution','rb') as av: arrayname17 = pickle.load(av)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/0.95_evolution','rb') as aw: arrayname18 = pickle.load(aw)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/1.0_evolution','rb') as ax: arrayname19 = pickle.load(ax)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/1.05_evolution','rb') as ay: arrayname20 = pickle.load(ay)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/1.1_evolution','rb') as az: arrayname21 = pickle.load(az)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/1.15_evolution','rb') as aa: arrayname22 = pickle.load(aa)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/1.2_evolution','rb') as ab: arrayname23 = pickle.load(ab)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/1.25_evolution','rb') as ac: arrayname24 = pickle.load(ac)
    
    main_array.append(arrayname1)
    main_array.append(arrayname2)
    main_array.append(arrayname3)
    main_array.append(arrayname4)
    main_array.append(arrayname5)
    main_array.append(arrayname6)
    main_array.append(arrayname7)
    main_array.append(arrayname8)
    main_array.append(arrayname9)
    main_array.append(arrayname10)
    main_array.append(arrayname11)
    main_array.append(arrayname12)
    main_array.append(arrayname13)
    main_array.append(arrayname14)
    main_array.append(arrayname15)
    main_array.append(arrayname16)
    main_array.append(arrayname17)
    main_array.append(arrayname18)
    main_array.append(arrayname19)
    main_array.append(arrayname20)
    main_array.append(arrayname21)
    main_array.append(arrayname22)
    main_array.append(arrayname23)
    main_array.append(arrayname24)
    
    return main_array
    
    
    
def calculate_prot(Proti_file_directory, RA_steps, DEC_steps, GUMS_file_directory):
    """
    This function calculates the rotation period of stars in document from basic stellar parameters.

    Proti_file_directory: File directory of document containing young cluster data aimed for the resampling of initial rotation periods.
    RA_steps: Steps of RA for file names.
    DEC_steps: Steps of DEC for file names.
    GUMS_file_directory: Path of file with GUMS data or stellar data.
    return: saves original csv file with an extra column for calculated Prots.
    """

    'INPUT'
    
    main_array=open_mass_files()
    
    kde_ProtMass = kde(Proti_file_directory)
    
    #Create array with initial rotation periods present in the grid
    Initial_Prots=np.arange(0.1,12.1,0.1) #0.1 steps
    #Create mass array with all the MIST masses we have
    MASSES = [0.1,0.15, 0.2, 0.25, 0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1,1.05,1.1,1.15, 1.2, 1.25] #all the MIST masses we have

    
    # Loop that looks into each file
    for k in RA_steps:
                            
        for b in DEC_steps:
    
            a=1
            
            while a<500:
                            
                found = common.find(f'RA_{k}_{k+4}_DEC_{b}_{b+4}_target.csv_Part{a}',GUMS_file_directory)
                
                if found == None:
                    
                    a=500
                
                else:
                    
                    print(f'RA_{k}_{k+4}_DEC_{b}_{b+4}_target.csv_Part{a}')
                    
                    #Open csv file containing star data
                    mrdata = pd.read_csv(found) #Read csv file
    
                    'BODY OF THE CODE'
    
                    #Create empty lists that we will use to create a new file with the calculated data
                    # M=[]
                    # AGE=[]
                    PROT=[]
    
                    #Do this loop for each star in file
                    n_star=len(mrdata['ra']) #Number of stars in file we want to evaluate
                    
                    for i in range(n_star):
                        
                        s=0
                        
                        while s<1:
                        
                            Prot_i=initial_Prot(kde_ProtMass)[0] #Initial rotation period
                            Grid_Prot_i1,Grid_Prot_i2=common.find_2_nearest(Initial_Prots, Prot_i) #Nearest initial rotation periods in grid
                            index1=int(np.where(Initial_Prots==Grid_Prot_i1)[0]) #Index of 1st nearest initial rotation period
                            index2=int(np.where(Initial_Prots==Grid_Prot_i2)[0]) #Index of 2nd nearest initial rotation period
                            mass=mrdata.mass[i] #Mass in solar masses
                            age=mrdata.uniform_ages[i]*1000 #Age in Myrs
                                    
                            Mass_1,Mass_2=common.find_2_nearest(MASSES, mass) #Find 2 nearest masses from MIST tables
                            
                            main_index1=int(np.where(MASSES==Mass_1)[0]) #Index of 1st nearest initial rotation period
                            main_index2=int(np.where(MASSES==Mass_2)[0]) #Index of 2nd nearest initial rotation period
                            
                            arrayname1=main_array[main_index1]       
                            arrayname2=main_array[main_index2]                        
                            
                            #Interpolate rotation period throughout whole life of star for the two MIST tables
                            spl_Prot_1 = InterpolatedUnivariateSpline(arrayname1[index1][0], arrayname1[index1][1], ext=0) #Interpolated line rotation period of Mass_1 and P_rot,i_1
                            spl_Prot_2 = InterpolatedUnivariateSpline(arrayname2[index1][0], arrayname2[index1][1], ext=0) #Interpolated line rotation period of Mass_2 and P_rot,i_1
                            spl_Prot_3 = InterpolatedUnivariateSpline(arrayname1[index2][0], arrayname1[index2][1], ext=0) #Interpolated line rotation period of Mass_1 and P_rot,i_2
                            spl_Prot_4 = InterpolatedUnivariateSpline(arrayname2[index2][0], arrayname2[index2][1], ext=0) #Interpolated line rotation period of Mass_2 and P_rot,i_2
                        
                                
                            #Calculate rotation period at specific age of star for the two MIST tables for different mass and different initial rotation period
                            interp_Prot_1=float(spl_Prot_1(age)) #Interpolated rotation period of Mass_1 and P_rot,i_1
                            interp_Prot_2=float(spl_Prot_2(age)) #Interpolated rotation period of Mass_2 and P_rot,i_1
                            interp_Prot_3=float(spl_Prot_3(age)) #Interpolated rotation period of Mass_1 and P_rot,i_2
                            interp_Prot_4=float(spl_Prot_4(age)) #Interpolated rotation period of Mass_2 and P_rot,i_2
                            
                            if np.isnan(interp_Prot_1) or np.isnan(interp_Prot_2) or np.isnan(interp_Prot_3) or np.isnan(interp_Prot_4):
                                s=0
                            else:
                                s=2
    
                        #Create arrays with the two initial rotation periods we are looking at and the rotation periods calculated for each mass
                        two_Prot_i=[Grid_Prot_i1,Grid_Prot_i2] #Initial rotation periods
                        diff_Prot1=[interp_Prot_1,interp_Prot_3] #For Mass_1
                        diff_Prot2=[interp_Prot_2,interp_Prot_4] #For Mass_2
                        
                        #Calculate rotation period for each mass for Prot_i by interpolating between the two calculated rotation periods for the two nearest initial rotation periods in grid
                        Med_Prot1 = common.interpolation(two_Prot_i, diff_Prot1, Prot_i)
                        Med_Prot2 = common.interpolation(two_Prot_i, diff_Prot2, Prot_i)
                        
                        #Create lists with the two nearest masses and the calculated rotation periods
                        two_masses=[Mass_1,Mass_2]
                        two_periods=[Med_Prot1,Med_Prot2]
                        
                        #Calculate the final rotation period by interpolating between the two previous results
                        Final_Prot = common.interpolation(two_masses, two_periods, mass)
                                
                        # AGE.append(age)
                        # M.append(mass)
                        PROT.append(Final_Prot)
                        
                        
                        if np.isnan(Final_Prot):
                            print('yes')
                            break
    
                    dictionary = {'Prot': PROT}  
                    dataframe = pd.DataFrame(dictionary) 
                    mrdata['Prot'] = dataframe
                    mrdata.to_csv(found, index=False)
                                    
                    a+=1