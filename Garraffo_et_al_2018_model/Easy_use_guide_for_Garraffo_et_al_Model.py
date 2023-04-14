'INPUT'

#Your data
Prot_init=np.array([3]) #Create an array of all initial rotation periods from which you want to start your rotation period evolution for a targeted mass
Mstar_init=[0.1] #Just one mass at a time


'BODY OF THE CODE'

#Run the evolution code
"""
Run the evolution code
Prot_evol: Array containing 1)Value of mass used, 2) Time steps of evolution, 3) Calculated rotation period at each time step
age_zero: Common age-array for the interpolation
Prot_interp: All rotation periods calculated by interpolation using age_zero time steps
spl_Prot: 1-D interpolating spline for the given Prot set of data points
"""
Prot_evol, age_zero, Prot_interp, spl_Prot = spin_down_evol(Prot_init, Mstar_init)


'OUTPUT'

#Evolution data
Masses_used=Prot_evol[0,0] #Mass used
Time_steps=Prot_evol[0,1] #All time steps of the evolution
Prot_evolution=Prot_evol[0,2] #The value of the rotation period for each time step

#Interpolated data
Interpolated_time_steps=age_zero.value #Common age-array for the interpolation
Interpolated_Prot_evolution=Prot_interp #ll rotation periods calculated by interpolation using age_zero time steps
Interpolating_function=spl_Prot #1-D interpolating spline for the given Prot set of data points
age=2500 #Age in Myrs
Prot_at_any_age=spl_Prot(age) #Calculated rotation period at any specific age using the interpolation function
