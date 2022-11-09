import numpy as np
import pandas as pd
import pickle
from scipy.interpolate import InterpolatedUnivariateSpline
import matplotlib.pyplot as plt

'Calculate Prot of GUMS stars'

"""
Code that uses the evolution grid to calculate the rotation period of stars with an initial rotation period in a specific age.

mass: Array with all masses of MIST tables you want to look at
name_mass: List of all exact values of MIST tables to use for the file names of the grid
RotationP: Array with all initial rotation periods you want to look at
"""

#Find nearest neighbour in array
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

#Find nearest 2 neighbours in array
def find_2_nearest(array, value):
    array = np.asarray(array)
    dif = np.abs(array - value)
    sort = np.argsort(dif)
    idx = sort[0] #Nearest neighbour
    idx2 = sort[1] #2nd nearest neighbour
    return array[idx], array[idx2]

#Interpolation between two data points
def interpolation(a, b, x):
    output = b[0] + (x - a[0]) * ((b[1] - b[0])/(a[1] - a[0]))
    return output


'INPUT'

#Open csv file containing star data
mrdata = pd.read_csv("/File/path/where/you/have/your/data/filename.csv") #Read csv file 


'BODY OF THE CODE'

#Create empty lists that we will use to create a new file with the calculated data
M=[]
AGE=[]
PROT=[]

#Create array with initial rotation periods present in the grid
Initial_Prots=np.arange(0.1,12.1,0.1) #0.1 steps
#Create mass array with all the MIST masses we have
MASSES = np.arange(0.1,2.05, 0.05)

#Do this loop for each star in file
n_star=1 #Number of stars in file we want to evaluate
for i in range(n_star):
    Prot_i=0.4 #Initial rotation period
    Grid_Prot_i1,Grid_Prot_i2=find_2_nearest(Initial_Prots, Prot_i) #Nearest initial rotation periods in grid
    index1=int(np.where(Initial_Prots==Grid_Prot_i1)[0]) #Index of 1st nearest initial rotation period
    index2=int(np.where(Initial_Prots==Grid_Prot_i2)[0]) #Index of 2nd nearest initial rotation period
    number_star=i
    mass=0.1#mrdata.mass[number_star] #Mass in solar masses
    age=mrdata.age[number_star]*1000 #Age in Myrs
        
    Mass_1,Mass_2=find_2_nearest(MASSES, mass) #Find 2 nearest masses from MIST tables
    
    #Open grid of the masses
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/{Mass_1}_evolution','rb') as f: arrayname1 = pickle.load(f)
    with open(f'/home/farah/Documents/Project/Data/Evolution_grid/{Mass_2}_evolution','rb') as g: arrayname2 = pickle.load(g)
    
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

    #Create arrays with the two initial rotation periods we are looking at and the rotation periods calculated for each mass
    two_Prot_i=[Grid_Prot_i1,Grid_Prot_i2] #Initial rotation periods
    diff_Prot1=[interp_Prot_1,interp_Prot_3] #For Mass_1
    diff_Prot2=[interp_Prot_2,interp_Prot_4] #For Mass_2
    
    #Calculate rotation period for each mass for Prot_i by interpolating between the two calculated rotation periods for the two nearest initial rotation periods in grid
    Med_Prot1 = interpolation(two_Prot_i, diff_Prot1, Prot_i)
    Med_Prot2 = interpolation(two_Prot_i, diff_Prot2, Prot_i)
    
    #Create lists with the two nearest masses and the calculated rotation periods
    two_masses=[Mass_1,Mass_2]
    two_periods=[Med_Prot1,Med_Prot2]
    
    #Calculate the final rotation period by interpolating between the two previous results
    Final_Prot = interpolation(two_masses, two_periods, mass)
            
    AGE.append(age)
    M.append(mass)
    PROT.append(Final_Prot)


'OUTPUT'
#Create a new csv file with the properties calculated
dictionary = {'age': AGE,'mass': M, 'Prot': PROT}  
dataframe = pd.DataFrame(dictionary) 
dataframe.to_csv('/File/path/where/you/want/to/save/the/new/file/filename.csv')
