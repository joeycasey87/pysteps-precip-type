
import os
import datetime
import numpy as np
import matplotlib.pyplot as plt

from pysteps.utils.reprojection import reproject_grids
#from pysteps.io import import_netcdf_pysteps
from pysteps.visualization import plot_precip_field, quiver

from ptype_functions import *
from IncaGribImport import IncaGribImporter
import imageio


# -----------------------------------------------
# plot function

def worker(gridPlot, metadata, i, start_date, dirGif):
    title = 'Date: %s' % (start_date.strftime("%Y-%m-%d %H:%M") + ' - ' + str(i))
    fig = plt.figure(figsize=(15, 15))
    plot_precip_field(gridPlot, geodata=metadata)
    plt.suptitle(title)
    plt.tight_layout()
    filename = f'{i}.png'
    plt.savefig(os.path.join(dirGif, filename), dpi=72)
    plt.close()
    return filename


# PLOT Precipitation Type INCA ---------------------------------------------------------------------------

startdate = datetime.datetime.strptime('202201090210', "%Y%m%d%H%M")
filename = r'C:\Users\talen\Desktop\Student JOB\Data\INCA-BE example files\09\precip\PT_FC_INCA.grb'

dir_gif = r'C:\Users\talen\Desktop\Student JOB\Data\gifs\\membersTEST2'

projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'

# Import GRIB
importer = IncaGribImporter()
incaDictionary_PT = importer.retrieve_grib_data(filename=filename, metadata_keys=keys)

# Transform to a 3D matrix
R_inca_PT = inca_dictionary_to_3Dmatrix(incaDictionary_PT)

# Fill GRIB metadata
metadata_inca = {}
metadata_inca['projection'] = projectionString
metadata_inca['xpixelsize'] = incaDictionary_PT['Metadata']['DxInMetres']
metadata_inca['ypixelsize'] = incaDictionary_PT['Metadata']['DyInMetres']
metadata_inca['x1'] = 360000.0
metadata_inca['y1'] = 350000.0
metadata_inca['x2'] = 960000.0
metadata_inca['y2'] = 940000.0
metadata_inca['yorigin'] = 'upper'


# PLOT over time and create GIF
filenames = []
titleDate = startdate
# ts = 0

# Plot members
for ts in range(R_inca_PT.shape[0]):
    # Plot
    filenames.append(plot_ptype(R_inca_PT[ts, :, :], metadata_inca, ts, titleDate, dir_gif))
    titleDate = titleDate + datetime.timedelta(minutes=10)

# Build gif
kargs = {'duration': 0.4}
with imageio.get_writer(
        os.path.join(dir_gif, (
            r'INCA_PT' + '_' + startdate.strftime('%Y%m%d%H%M') + '.gif')), mode='I',
        **kargs) as writer:
    for filename in filenames:
        image = imageio.imread_v2(os.path.join(dir_gif, filename))
        writer.append_data(image)

# Close gif writer
writer.close()

# Remove temporary files
for filename in set(filenames):
    os.remove(os.path.join(dir_gif, filename))



 



# PLOT TOPO # ---------------------------------------------------------------------------

startdate = datetime.datetime.strptime('202201090210', "%Y%m%d%H%M")
topoFilename = r'C:\Users\talen\Desktop\Student JOB\Data\INCA-TOPO\inca_topo.asc'
dir_gif = r'C:\Users\talen\Desktop\Student JOB\Data\gifs\\membersTEST2'

inca_topo_grid = np.loadtxt(topoFilename)
inca_topo_grid = inca_topo_grid[::-1, :]  # Reorientation
print('Topography load done')

metadata_inca = {}
metadata_inca['projection'] = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
metadata_inca['xpixelsize'] = 1000
metadata_inca['ypixelsize'] = 1000
metadata_inca['x1'] = 360000.0
metadata_inca['y1'] = 350000.0
metadata_inca['x2'] = 960000.0
metadata_inca['y2'] = 940000.0
metadata_inca['yorigin'] = 'upper'

# plot (precipitation plot is been used here se a different plot function)
worker(inca_topo_grid, metadata_inca, 0, startdate, dir_gif)


