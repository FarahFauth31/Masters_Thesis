## Evolution grid directory

This directory contains the files of the Evolution Grid created with the Garraffo et al. 2018 model.
Each document corresponds to a mass from the MIST tables.
Each document contains the evolution of the rotation period of a star of a specific mass value for different initial rotation periods.
The masses range from 0.1 to 1.25 M_sol in steps of 0.05 M_sol.
The initial rotation periods range from 0.1 to 12 days in steps of 0.1 days.
The directory also contains graphs depicting the content of each file.

There are also two python files. Evolution_grip.py creates an evolution grid. It saves files for each mass that contain the rotation period evolution for all initial rotation periods given. Using_the_evolution_grid.py uses the evolution grid to calculate the rotation period of stars with an initial rotation period in a specific age. It can be ignored since there is a better example in /xuv_milky_way/Examples.ipynb .
