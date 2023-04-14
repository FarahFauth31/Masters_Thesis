## Scripts

1. Easy_use_guide_for_Garraffo_et_al_Model
    
   The code of the Garraffo et al. 2018 model can be very tricky to understand. That's why I made an easy-to-use guide. This script can be used for creating the rotation period evolution of different initial rotation periods of a star of a particular mass. The only thing you have to do is input your data for the mass of the star and an array with different initial rotation periods. To use this script, the Garraffo et al. 2018 model code and MIST data have to be downloaded first. The original version of the code can be found [here](https://github.com/FarahFauth31/Year4Project/blob/main/Garraffo_et_al_2018_Model). The MIST data can be found [here](https://github.com/cgarraffo/Spin-down-model/tree/master/Mist_data). An updated version of the code can be found [here](https://github.com/FarahFauth31/Year4Project/blob/main/Variation_of_Garraffo_et_al_Model).
   
2. Garraffo_et_al_2018_Model

   Original version of the Garraffo et al. 2018 model. To use it you need MIST data found [here](https://github.com/cgarraffo/Spin-down-model/tree/master/Mist_data).
   
3. Prot_i_3d_evolution_0.1M and Prot_i_3d_evolution_0.15M

   Files used by Variation_of_Garraffo_et_al_Model. The code uses them to smooth some tracks for very low mass stars (0.1-0.15M).
   
4. Variation_of_Garraffo_et_al_Model

   Updated version of the Garraffo et al. 2018 model that includes some changes to make the output less confusing and more reliable (e.g. it outputs the rotation period evolution of the star until the specified age, smoothing of evolution tracks, etc.). To run the code you need to download the MIST data and the Prot_i_3d_evolution_0.1M and Prot_i_3d_evolution_0.15M files.

