Year4Project

# Simulating the Stellar Galactic X-ray Emission of Cool Stars for Exoplanet Habitability

*Note*: This is my Master's project. I have been working on it from the beginning of August. It's still an ongoing process. Here I present some of the scripts I have been creating over the past months.

## Abstract

Stellar coronal X-ray and EUV emissions heat and evaporate planetary atmospheres. Understanding their evolution through
time can give us new insight into exoplanet habitability. In this work, we make a first attempt at simulating the X-ray emission of
a realistic synthetic stellar galactic population by employing the Gaia Universe Model Snapshot (GUMS) model, together with
the Garraffo et al. (2018) spin-down model and the Wright et al. (2011) model of stellar X-ray emission (𝐿 𝑥 ). Using observed
𝐿 𝑥 data, we will test the stellar rotation and X-ray emission models, both statistically (comparing GUMS to X-ray survey data)
and through matching Gaia data of stars with determined masses and ages with Chandra observations (using machine learning
techniques). This takes us a step forward in building a trustworthy model that describes these stars’ rotation rate and the resulting
X-ray emission.

## Scripts

1. Easy_use_guide_for_Garraffo_et_al_Model
    
   The code of the Garraffo et al. 2018 model can be very tricky to understand. That's why I made an easy-to-use guide. This script can be used for creating the rotation period evolution of different initial rotation periods of a star of a particular mass. The only thing you have to do is input your data for the mass of the star and an array with different initial rotation periods. To use this script, the Garraffo et al. 2018 model code and MIST data have to be downloaded first. The original version of the code can be found [here](https://github.com/FarahFauth31/Year4Project/blob/main/Garraffo_et_al_2018_Model). The MIST data can be found [here](https://github.com/cgarraffo/Spin-down-model/tree/master/Mist_data). An updated version of the code can be found [here](https://github.com/FarahFauth31/Year4Project/blob/main/Variation_of_Garraffo_et_al_Model).
    
2. Evolution_grid

   This code is used to create an evolution grid. It saves files for each possible stellar mass that contain the rotation period evolution for all initial rotation periods given. To use this code you need the Garraffo et al. 2018 model code and MIST data (found [here](https://github.com/FarahFauth31/Year4Project/blob/main/Garraffo_et_al_2018_Model) and [here](https://github.com/cgarraffo/Spin-down-model/tree/master/Mist_data)).
   
3. Garraffo_et_al_2018_Model

   Original version of the Garraffo et al. 2018 model. To use it you need MIST data found [here](https://github.com/cgarraffo/Spin-down-model/tree/master/Mist_data).
   
4. Prot_i_3d_evolution_0.1M and Prot_i_3d_evolution_0.15M

   Files used by Variation_of_Garraffo_et_al_Model. The code uses them to smooth some tracks for very low mass stars (0.1-0.15M).
   
5. Using_the_evolution_grid

   Code that uses the evolution grid to calculate the rotation period of stars with an initial rotation period in a specific age. Your input has to include the mass, initial rotation period and age value of each mass. You can do that by creating a csv file with all the information for each star and run the code.
   
6. Variation_of_Garraffo_et_al_Model

   Updated version of the Garraffo et al. 2018 model that includes some changes to make the output less confusing and more reliable (e.g. it outputs the rotation period evolution of the star until the specified age, smoothing of evolution tracks, etc.). To run the code you need to download the MIST data and the Prot_i_3d_evolution_0.1M and Prot_i_3d_evolution_0.15M files.
   
## Roadmap

I am aiming to publish some new scripts to be able to calculate the X-ray emission of stars using the rotation period found with the Garraffo et al. 2018 model. The X-ray emission wil be modelled using the Wright et al. 2011 model. Another extension to this project will be to use the results obtained to create the full Xray-EUV spectrum of the studied stars.

## Authors and acknowledgment

Farah Fauth Puigdomenech
Cecilia Garraffo
Jeremy Drake
Kristina Monsch
Vinay Kashyap
Rafael Martinez Galarza


