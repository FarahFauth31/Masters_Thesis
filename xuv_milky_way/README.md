## XUV Milky Way python package

This directory contains the files with all the functions necessary to calculate Prot, Xray and EUV emissions of stars from basic stellar parameters.

common.py contains the functions used by all the other functions in the oher python files.

prot.py contains functions to treat stellar data and calculate the rotation period of star from mass, age and initial rotation period.

xray.py contains functions to calculate the Xray emission of stars from rotation period, convective turnover time and bolometric luminosity.

euv.py contains functions to calculate the EUV emissions of stars and assign to each star a spectra depending on its Xray surface flux.

Examples.ipynb gives an example of how to use all the python functions of the python scripts described above to calculate Prot, Xray and EUV from basic stellar parameters found in a .csv file.
            
