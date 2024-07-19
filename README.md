This repository contains functions to allow for the calculation of precipitation type from a combination of a pysteps nowcast and data from an independent weather model, such as INCA or COSMO.
These functions can be found in the prtype python script.
The data required is a 2D array for each of the snow level, temperature, and ground temperature, as well as the metadata associated with the weather model.
In general, individual importers for each weather model will be required to extract this data from the model files.
An example of one such importer can be seen in the IncaGribImport script which can be used to extract data from .grb files produced by the INCA weather model.
An example of how to use the precipitation type calculator, as well as how to create plots and a gif from the output can be seen in the testing notebook.

There is a question mark over where to include this functionality in pysteps. It could be directly added to the postprocessing folder (alongside the ensemblestats and probmatching scripts), 
or it could be added as a plugin. However, there is currently no structure in place to create plugins in the postprocessing folder so this would have to be developed in order to add the precip-type script as a plugin.
