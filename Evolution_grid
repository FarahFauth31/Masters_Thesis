"""
Code that creates an evolution grid. It saves files for each mass that contain the rotation period evolution for all initial rotation periods given

mass: Array with all masses of MIST tables you want to look at
name_mass: List of all exact values of MIST tables to use for the file names of the grid
RotationP: Array with all initial rotation periods you want to look at
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import pickle

'INPUT'

#Array with all masses of MIST tables you want to look at
mass = np.arange(0.1,1.3,0.05) #0.05 steps
name_mass = [0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1, 1.05, 1.1, 1.15, 1.2, 1.25] #exact masses

#Array with all initial rotation periods you want to look at
RotationP=np.arange(0.1,12.1,0.1) #0.1 steps


'BODY OF THE CODE'      
     
m_sample=len(mass) #Number of masses we want to look at from the mass array
p_sample=len(RotationP) #Number of rotation periods we want to look at from the RotationP array

#Empty data one unit arrays so they can be used by the model
Prot_new=np.zeros(1) 
M_new=np.zeros(1)

#Creating a loop to create the evolutionary data we need by taking each initial period and running the model for each mass
for c in range(m_sample): #For each mass
    M_new[0]=find_nearest(MASSES, mass[c]) #Put mass number in the single unit mass array
    a=load_mist_tables(Mstar=find_nearest(MASSES, mass[c])) #Look at how many time steps the MIST table has
    Final_set=np.zeros((p_sample,2,len(a[0]))) #Empty array with final data with rotation period evolution over time

    for i in range(p_sample): #For each initial rotation period
        Prot_new[0]=RotationP[i] #Put rotation period number in the single unit Prot array
        # t0: time [Myr] at which the model should start 
        # tdisc: disk-locking time [Myr]
        Prot_evol, age_zero, Prot_interp, spl_Prot = spin_down_evol(Prot_init=Prot_new, 
                                                                Mstar_init=M_new, 
                                                                t0=1., tdisc=13.)
        Final_set[i][0]=Prot_evol[0,1] #Evolution of time
        Final_set[i][1]=Prot_evol[0,2] #Evolution of rotation period
    
    'OUTPUT'
    #Save data for each mass
    with open(f'/File/path/where/you/want/to/save/the/grid/files/{name_mass[c]}_evolution','wb') as f: pickle.dump(Final_set, f) #Save data for one mass as a pickle file
