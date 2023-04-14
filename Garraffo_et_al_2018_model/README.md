## Scripts

Easy_use_guide_for_Garraffo_et_al_Model

The code of the Garraffo et al. 2018 model can be very tricky to understand. That's why I made an easy-to-use guide. This script can be used for creating the rotation period evolution of different initial rotation periods of a star of a particular mass. The only thing you have to do is input your data for the mass of the star and an array with different initial rotation periods. To use this script, the Garraffo et al. 2018 model code and MIST data have to be downloaded first. The original version of the code can be found here. The MIST data can be found here. An updated version of the code can be found here.

Evolution_grid

This code is used to create an evolution grid. It saves files for each possible stellar mass that contain the rotation period evolution for all initial rotation periods given. To use this code you need the Garraffo et al. 2018 model code and MIST data (found here and here).

Garraffo_et_al_2018_Model

    Original version of the Garraffo et al. 2018 model. To use it you need MIST data found here.

    Prot_i_3d_evolution_0.1M and Prot_i_3d_evolution_0.15M

    Files used by Variation_of_Garraffo_et_al_Model. The code uses them to smooth some tracks for very low mass stars (0.1-0.15M).

    Using_the_evolution_grid

    Code that uses the evolution grid to calculate the rotation period of stars with an initial rotation period in a specific age. Your input has to include the mass, initial rotation period and age value of each mass. You can do that by creating a csv file with all the information for each star and run the code.

    Variation_of_Garraffo_et_al_Model

    Updated version of the Garraffo et al. 2018 model that includes some changes to make the output less confusing and more reliable (e.g. it outputs the rotation period evolution of the star until the specified age, smoothing of evolution tracks, etc.). To run the code you need to download the MIST data and the Prot_i_3d_evolution_0.1M and Prot_i_3d_evolution_0.15M files.

