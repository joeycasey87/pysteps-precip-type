# -*- coding: utf-8 -*-
"""
This is an importer for the precipitation type module.

This plugin is a conversion of the precipitation type calculator script adapted for pySTEPS from the INCA precipitation
type calculator.
"""

# Import the needed libraries
import numpy as np
from pysteps.decorators import postprocess_import
import os
import datetime
import imageio
import time

from ptype_functions import *
from IncaGribImport import IncaGribImporter

from pysteps.utils.reprojection import reproject_grids
from pysteps.io import import_netcdf_pysteps

@postprocess_import()
def importer_rmi_xxx(filename,
                     startdate,
                     gribPath,
                     topoFilename,
                     dir_gif,
                     keys=None,
                     inca_projectionString=None,
                     nwc_projectionString=None,
                     members=None,
                     timeBase=None,
                     timeStep=None,
                     ** kwargs):
    """
    Import a script calculating the precipitation types of hydrometeors from a combination of pySTEPS data and INCA
    data.

    Parameters
    ----------
    filename : str
        Name of the NetCDF file to import. The NetCDF files and the Grib files must be for the same date and time.

    startdate : datetime
        The time and date of the inca grib files in the format "%Y%m%d%H%M" e.g 202305010000 for midnight on the 1st of
        May 2023.

    gribPath : str
        The file path to the directory which stores the inca grib files. Grib files should be of the .grb format.
        Snow level (ZS), temperature (TT), and ground temperature (TG) files are required.

    topoFilename : str
        The path to the INCA topography file. The topography file is required to be in a text readable format such as
        .asc

    dir_gif : str
        The directory in which to store the output GIF.

    keys : list
        A list of the names of the GRIB import selection keys.

    inca_projectionString : str
        The INCA projection string.

    nwc_projectionString : str
        The NWC projection string.

    members : int
        The number of members to plot over time.

    timeBase : int
        The base for the time period. (min)

    timeStep : int
        The step for the time period. (min)

    {extra_kwargs_doc}

    Returns
    -------
    precipitation : 2D array, float32
        Precipitation field in mm/h. The dimensions are [latitude, longitude].
    quality : 2D array or None
        If no quality information is available, set to None.
    metadata : dict
        Associated metadata (pixel sizes, map projections, etc.).
    """

    ####################################################################################
    # Add the code to read the precipitation data here. Note that only cartesian grid
    # are supported by pysteps!
    # ---------------------------------------------------------------------------
    # READ inca and pysteps files.
    # Format, Reproject, interpolate and plot
    # ---------------------------------------------------------------------------
    # --------------------------------------------------------------------------
    # ---------------------------------------------------------------------------

    # Define default parameter values
    if keys is None:
        keys = ['Nx', 'Ny', 'latitudeOfFirstGridPointInDegrees', 'longitudeOfFirstGridPointInDegrees', 'earthIsOblate',
                'LoVInDegrees', 'DxInMetres', 'DyInMetres', 'iScansNegatively', 'jScansPositively', 'Latin1InDegrees',
                'Latin2InDegrees']
    if inca_projectionString is None:
        inca_projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'
    if nwc_projectionString is None:
        nwc_projectionString = '+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=50.797815 +lon_0=4.359215833333333 +x_0=649328 +y_0=665262 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs '
    if members is None:
        members = 1
    if timeBase is None:
        timeBase = 60
    if timeStep is None:
        timeStep = 5

    fn_zs = '_ZS_FC_INCA.grb'
    fn_tt = '_TT_FC_INCA.grb'
    fn_tg = '_TG_FC_INCA.grb'

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
    r_nwc, metadata_nwc = import_netcdf_pysteps(filename)

    # Set Metadata info
    metadata_nwc['projection'] = nwc_projectionString
    metadata_nwc['cartesian_unit'] = metadata_nwc['projection'][nwc_projectionString.find('units=') + 6:nwc_projectionString.find(' +no_defs', nwc_projectionString.find('units=') + 6)]
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
    # The topography file is not required to be interpolated if it is already of the shape matching the grib files,
    # i.e (591, 601)
    #topo_grid, _ = reproject_grids(np.array([topo_grid]), r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
    print('Reprojection done')

    # --------------------------------------------------------------------------
    # Calculate interpolation matrices

    # Calculate interpolations values for matching timestamps between INCA and pySTEPS
    inca_interpolations_ZS, timestamps_idxs = generate_inca_interpolations(R_inca_ZS, metadata_nwc['timestamps'],
                                                                           startdate, timeStep, timeBase)
    inca_interpolations_TT, _ = generate_inca_interpolations(R_inca_TT, metadata_nwc['timestamps'], startdate, timeStep,
                                                             timeBase)
    inca_interpolations_TG, _ = generate_inca_interpolations(R_inca_TG, metadata_nwc['timestamps'], startdate, timeStep,
                                                             timeBase)
    print("Interpolation done!")

    # Clean (After interpolation, we don't need the reprojected data anymore)
    del R_inca_ZS, R_inca_TT, R_inca_TG

    # --------------------------------------------------------------------------
    # Diagnose precipitation type per member over time, using mean mask

    # WARNING (1): The grids have been sub-scripted to a smaller size (INCA 590x600). This requires the inca_metadata to
    # be used for plotting. If the original PYSTEPS grid size is used (700x700) for plotting, the pysteps metadata_nwc
    # should be used instead.
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
                    r'INCA_mem_' + ('mean_' if member == mean_idx else '') + str(member) + '_' + startdate.strftime(
                '%Y%m%d%H%M') + '.gif')), mode='I',
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

    plot_ptype(np.array(ptype_mean), metadata_inca, 0, timestamps_idxs[0], dir_gif)

    precip = np.array(ptype_mean)

    quality = None

    projection_definition = inca_projectionString

    metadata = dict(
        xpixelsize=1,
        ypixelsize=1,
        cartesian_unit=metadata_nwc['cartesian_unit'],
        unit="mm/h",
        transform=None,
        zerovalue=0,
        institution="rmi",
        projection=projection_definition,
        yorigin="upper",
        threshold=0.03,
        x1=0,
        x2=100,
        y1=0,
        y2=100,
    )

    # IMPORTANT! The importers should always return the following fields:
    return precip, quality, metadata
