## paper_Winds_AtmPatterns_Seasonal_Interannual_TTS_MultipleInlet
This repository contains scripts and data for reproducing the figures from the manuscript:

> Fajardo-Urbina, J. M., Arts, G., Grawe, U., Herman, J. H., Gerkema, T., and Duran-Matute, M. (2022): Atmospherically Driven Seasonal and Interannual Variability in the Lagrangian Transport Time Scales of a Multiple-inlet Coastal System. 

### Data
There are in total 8 NetCDF files which can be used to reproduce the 3 figures from the main manuscript and the 5 figures from the supporting information. They are provided inside the folder [data](https://github.com/JeancarloFU/paper_Winds_AtmPatterns_Seasonal_Interannual_TTS_MultipleInlet/blob/main/data).

### Sofware
The environment employed for the analysis is based on Python v3.8 and can be found in the file [environment.yml](https://github.com/JeancarloFU/paper_Winds_AtmPatterns_Seasonal_Interannual_TTS_MultipleInlet/blob/main/environment.yml).
The wavelet analysis is based on the Python package [Pycwt v0.3.0a22](https://anaconda.org/conda-forge/pycwt), but we added a script to perform the bias correction ([Liu et al., 2007](https://journals.ametsoc.org/view/journals/atot/24/12/2007jtecho511_1.xml)) and a wavelet filter ([Torrence & Compo, 1998](https://journals.ametsoc.org/view/journals/bams/79/1/1520-0477_1998_079_0061_apgtwa_2_0_co_2.xml)). This script (called [get_results.py](https://github.com/JeancarloFU/paper_Winds_AtmPatterns_Seasonal_Interannual_TTS_MultipleInlet/tree/main/functions/pycwt_wavelet_package/get_results.py)) and the core functions from Pycwt are in [/functions/pycwt_wavelet_package](https://github.com/JeancarloFU/paper_Winds_AtmPatterns_Seasonal_Interannual_TTS_MultipleInlet/tree/main/functions/pycwt_wavelet_package).

### Running the Notebooks
All the scripts (notebooks) necessary for reproducing the figures are in the folder [notebooks](https://github.com/JeancarloFU/paper_Winds_AtmPatterns_Seasonal_Interannual_TTS_MultipleInlet/blob/main/notebooks). They can be run using any of the following instructions:
- Mybinder.org: click on the binder icon to open a jupyter-lab session. 
- Google Colab: follow the instructions from the notebook [clone_repo_using_google_colab.ipynb](https://github.com/JeancarloFU/paper_Winds_AtmPatterns_Seasonal_Interannual_TTS_MultipleInlet/blob/main/clone_repo_using_google_colab.ipynb)
- Clone or download the repository in your PC: install the packages of the file [environment.yml](https://github.com/JeancarloFU/paper_Winds_AtmPatterns_Seasonal_Interannual_TTS_MultipleInlet/blob/main/environment.yml).

**Jup-lab:** [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/JeancarloFU/paper_Winds_AtmPatterns_Seasonal_Interannual_TTS_MultipleInlet/main?urlpath=lab)

### Information about raw numerical data 
The 8 netCDF files used to generate the figures, were obtained from the following raw data:

- Eulerian data from the GETM/GOTM model, and its set-up is described in:
    * Duran-Matute et al. (2014): https://doi.org/10.5194/os-10-611-2014
    * Grawe et al. (2016): https://doi.org/10.1002/2016JC011655
- Lagrangian data from the model Parcels v2.1.1, which can be installed from: 
    * https://anaconda.org/conda-forge/parcels
    * https://oceanparcels.org
- Monthly mean sea level pressure (MSLP) from NCEP-NCAR Reanalysis 1 (used for the large-scale analytical model):
    * https://psl.noaa.gov/data/gridded/data.ncep.reanalysis.html