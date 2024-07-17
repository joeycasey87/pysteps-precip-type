# Import the precipitation type importer
from importer_rmi import *

# Specify the input parameters

# NetCDF4 file
filename = "C:\\Users\\josep\\OneDrive\\Documents\\Internship_RMI\\Hackathon_test_data\\blended_nowcast_date_shifted.nc"
# Start date and time of the data
startdate = datetime.datetime.strptime('202201080000', "%Y%m%d%H%M")
# Path to the directory where the INCA data files are stored
gribPath = "C:\\Users\\josep\\OneDrive\\Documents\\Internship_RMI\\Hackathon_test_data"
# INCA topography file path
topo = "C:\\Users\\josep\\OneDrive\\Documents\\Internship_RMI\\Hackathon_test_data\\inca_dem.asc"
# The directory in which the GIF shall be saved
dir_gif = "C:\\Users\\josep\\OneDrive\\Documents\\Internship_RMI\\Hackathon_test_data"
# A list of the names of the GRIB import selection keys
keys = ['Nx', 'Ny', 'latitudeOfFirstGridPointInDegrees', 'longitudeOfFirstGridPointInDegrees', 'earthIsOblate',
        'LoVInDegrees', 'DxInMetres', 'DyInMetres', 'iScansNegatively', 'jScansPositively', 'Latin1InDegrees',
        'Latin2InDegrees']
# The projection string for the INCA data
inca_projection_string = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
# The projection string for the pySTEPS data
nwc_projection_string =  '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
# The number of members to plot over time
members = 2
# The base length of time over which to plot
timeBase = 60
# The length of the time steps
timeStep = 5

array_list = precipitation_type_calculator(filename, startdate, gribPath, topo, keys, inca_projection_string, nwc_projection_string, members, timeBase, timeStep)

print(array_list)
# Build gif
# kargs = {'duration': 0.4}
# with imageio.get_writer(
#     os.path.join(dir_gif, (
#             r'INCA_mem_' + ('mean_' if member == mean_idx else '') + str(member) + '_' + startdate.strftime(
#         '%Y%m%d%H%M') + '.gif')), mode='I',
#     **kargs) as writer:
#         for filename in filenames:
#                 image = imageio.imread_v2(os.path.join(dir_gif, filename))
#                 writer.append_data(image)
#
# # Close gif writer
# writer.close()