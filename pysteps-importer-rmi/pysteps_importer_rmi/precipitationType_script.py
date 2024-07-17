# ---------------------------------------------------------------------------
# READ inca and pysteps files.
# Format, Reproject, interpolate and plot
# ---------------------------------------------------------------------------

import os
import datetime
import numpy as np
import imageio
import time

from ptype_functions import *
from IncaGribImport import IncaGribImporter

from pysteps.utils.reprojection import reproject_grids
from pysteps.io import import_netcdf_pysteps

# --------------------------------------------------------------------------
# Parameters

# Start date and time
startdate = datetime.datetime.strptime('202201080000', "%Y%m%d%H%M")

# File paths and names
gribPath = r"C:\\Users\\josep\\OneDrive\\Documents\\Internship_RMI\\Hackathon_test_data"
fn_zs = '_ZS_FC_INCA.grb'
fn_tt = '_TT_FC_INCA.grb'
fn_tg = '_TG_FC_INCA.grb'
#filename = NetCDFFilename_1 + NetCDFFilename_2
netCDF4Filename = r"C:\\Users\\josep\\OneDrive\\Documents\\Internship_RMI\\Hackathon_test_data\\blended_nowcast.nc"
topoFilename = r"C:\\Users\\josep\\OneDrive\\Documents\\Internship_RMI\\Hackathon_test_data\\inca_dem.asc"

# GIF output folder
dir_gif = r"C:\\Users\\josep\\OneDrive\\Documents\\Internship_RMI\\Hackathon_test_data"

# GRIB import selection keys
keys = ['Nx', 'Ny', 'latitudeOfFirstGridPointInDegrees', 'longitudeOfFirstGridPointInDegrees', 'earthIsOblate',
        'LoVInDegrees', 'DxInMetres', 'DyInMetres', 'iScansNegatively', 'jScansPositively', 'Latin1InDegrees',
        'Latin2InDegrees']

# Projection strings
inca_projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
nwc_projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '

# Nr of members to plot (over time)
members = 1

# Time period settings (these are default values in all function parameters)
timeBase = 60
timeStep = 5


# ---------------------------------------------------------------------------
# Import ZS Snow level

filename_ZS = os.path.join(gribPath, startdate.strftime("%Y%m%d%H%M") + fn_zs)
importer = IncaGribImporter()
incaDictionary_ZS = importer.retrieve_grib_data(filename=filename_ZS, metadata_keys=keys)

# Transform to a 3D matrix
R_inca_ZS = inca_dictionary_to_3Dmatrix(incaDictionary_ZS)
print('INCA Snow level load done')

del importer

# ---------------------------------------------------------------------------
# Import TT Temperature

filename_TT = os.path.join(gribPath, startdate.strftime("%Y%m%d%H%M") + fn_tt)
importer = IncaGribImporter()
incaDictionary_TT = importer.retrieve_grib_data(filename=filename_TT, metadata_keys=keys)

# Transform to a 3D matrix
R_inca_TT = inca_dictionary_to_3Dmatrix(incaDictionary_TT)

# Transform Kelvin to Celsius
R_inca_TT[:, :, :] = R_inca_TT[:, :, :] - 273.15
print('INCA temperature load done')

del importer

# ---------------------------------------------------------------------------
# Import TG Ground temperature

filename_TG = os.path.join(gribPath, startdate.strftime("%Y%m%d%H%M") + fn_tg)
importer = IncaGribImporter()
incaDictionary_TG = importer.retrieve_grib_data(filename=filename_TG, metadata_keys=keys)

# Transform to a 3D matrix
R_inca_TG = inca_dictionary_to_3Dmatrix(incaDictionary_TG)

# Transform Kelvin to Celsius
R_inca_TG[:, :, :] = R_inca_TG[:, :, :] - 273.15
print('INCA Ground temperature load done')

del importer

# ---------------------------------------------------------------------------
# Build inca metadata (this should be an output of the INCA importer)

metadata_inca = {}
metadata_inca['projection'] = inca_projectionString
metadata_inca['xpixelsize'] = incaDictionary_ZS['Metadata']['DxInMetres']
metadata_inca['ypixelsize'] = incaDictionary_ZS['Metadata']['DyInMetres']
metadata_inca['x1'] = 360000.0
metadata_inca['y1'] = 350000.0
metadata_inca['x2'] = 960000.0
metadata_inca['y2'] = 940000.0
metadata_inca['yorigin'] = 'upper'

# --------------------------------------------------------------------------
# Load INCA Topography (this might be a different file format in the future)

topo_grid = np.loadtxt(topoFilename)
topo_grid = topo_grid[::-1, :]  # Reorientation
print('Topography load done')

# Clean
del incaDictionary_ZS, incaDictionary_TT, incaDictionary_TG

# ---------------------------------------------------------------------------
# Load PYSTEPS data

# import netCDF file
r_nwc, metadata_nwc = import_netcdf_pysteps(netCDF4Filename)

# Set Metadata info
metadata_nwc['projection'] = nwc_projectionString
metadata_nwc['cartesian_unit'] = 'm'
print('netCDF4 load done')

# ---------------------------------------------------------------------------
# measure time
start_time = time.time()

# --------------------------------------------------------------------------
# Reproject

# INCA files reprojection over pySTEPS grid
R_inca_ZS, _ = reproject_grids(R_inca_ZS, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
R_inca_TT, _ = reproject_grids(R_inca_TT, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
R_inca_TG, _ = reproject_grids(R_inca_TG, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)

# No need to reproject topography file (590x600)
# topo_grid, _ = reproject_grids(np.array([topo_grid]), r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
print('Reprojection done')

# --------------------------------------------------------------------------
# Calculate interpolation matrices

# Calculate interpolations values for matching timestamps between INCA and pySTEPS
inca_interpolations_ZS, timestamps_idxs = generate_inca_interpolations(R_inca_ZS, metadata_nwc['timestamps'], startdate, timeStep, timeBase)
inca_interpolations_TT, _ = generate_inca_interpolations(R_inca_TT, metadata_nwc['timestamps'], startdate, timeStep, timeBase)
inca_interpolations_TG, _ = generate_inca_interpolations(R_inca_TG, metadata_nwc['timestamps'], startdate, timeStep, timeBase)
print("Interpolation done!")

# Clean (After interpolation, we don't need the reprojected data anymore)
del R_inca_ZS, R_inca_TT, R_inca_TG


# --------------------------------------------------------------------------
# Diagnose precipitation type per member over time, using mean mask

# WARNING (1): The grids have been sub-scripted to a smaller size (INCA 590x600). This requires the inca_metadata to be
# used for plotting. If the original PYSTEPS grid size is used (700x700) for plotting, the pysteps metadata_nwc should
# be used instead.
#
# WARNING (2): Topography original size is 590x600. This grid was not reprojected.

print("Calculate precipitation type per member over time...")

# Find subscript indexes for INCA grid (590x600)
x1, x2, y1, y2 = get_reprojected_indexes(inca_interpolations_ZS[0])

# Result list
ptype_list = np.zeros((r_nwc.shape[0] + 1, r_nwc.shape[1], x2 - x1, y2 - y1))

# loop over timestamps
for ts in range(len(timestamps_idxs)):
    print("Calculating precipitation types at: ", str(timestamps_idxs[ts]))

    # Members Mean matrix
    r_nwc_mean = calculate_members_mean(r_nwc[:, ts, x1:x2, y1:y2])

    # calculate precipitation type result with members mean
    ptype_mean = calculate_precip_type(incaZnow=inca_interpolations_ZS[ts, x1:x2, y1:y2],
                                       incaTemp=inca_interpolations_TT[ts, x1:x2, y1:y2],
                                       incaGroundTemp=inca_interpolations_TG[ts, x1:x2, y1:y2],
                                       precipGrid=r_nwc_mean,
                                       topographyGrid=topo_grid)

    # Intersect precipitation type by member using ptype_mean
    for member in range(r_nwc.shape[0]):
        res = np.copy(ptype_mean)
        res[r_nwc[member, ts, x1:x2, y1:y2] == 0] = 0
        ptype_list[member, ts, :, :] = res

    # Add mean result at the end
    ptype_list[-1, ts, :, :] = ptype_mean

print("--Script finished--")
print("--- %s seconds ---" % (time.time() - start_time))


# --------------------------------------------------------------------------
# PLOT (single member over time)

# measure time
start_time = time.time()

# Choose 1 member to plot
member = 0
# Members mean is always stored at the last index (used for the file name only)
mean_idx = ptype_list.shape[0] - 1

# Plot members
filenames = []
for ts in range(len(timestamps_idxs)):
    # Plot
    filenames.append(plot_ptype(ptype_list[member, ts, :, :], metadata_inca, ts, timestamps_idxs[ts], dir_gif))

# Build gif
kargs = {'duration': 0.4}
with imageio.get_writer(
        os.path.join(dir_gif, (
            r'INCA_mem_' + ('mean_' if member == mean_idx else '') + str(member) + '_' + startdate.strftime('%Y%m%d%H%M') + '.gif')), mode='I',
        **kargs) as writer:
    for filename in filenames:
        image = imageio.imread_v2(os.path.join(dir_gif, filename))
        writer.append_data(image)

# Close gif writer
writer.close()

# Remove temporary files
for filename in set(filenames):
    os.remove(os.path.join(dir_gif, filename))

print("--- finished plotting ---")
print("--- %s seconds ---" % (time.time() - start_time))


# test plots
# plot_ptype(np.array(ptype_mean), metadata_inca, 0, timestamps_idxs[0], dir_gif)
# plot_ptype(np.array(R_inca_ZS[0]), metadata_inca, 0, timestamps_idxs[0], dir_gif)
