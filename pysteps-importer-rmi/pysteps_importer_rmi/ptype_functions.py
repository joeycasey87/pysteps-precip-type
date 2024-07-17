# Import
import os
import datetime
import time
import numpy as np
import matplotlib.pyplot as plt
from precipitationTypeFields import plot_precipType_field

from ptype_functions import *
from IncaGribImport import IncaGribImporter

from pysteps.utils.reprojection import reproject_grids
from pysteps.io import import_netcdf_pysteps


# Dictionary to matrix function (this can be added to the importer class)
def inca_dictionary_to_3Dmatrix(incaDict):
    resultMatrix = np.empty(shape=(
        len(incaDict['Messages'].keys()), incaDict['Messages'][1]['Grid'].shape[0],
        incaDict['Messages'][1]['Grid'].shape[1]))
    for i in range(len(resultMatrix)):
        resultMatrix[i, :, :] = incaDict['Messages'][i + 1]['Grid']
    return resultMatrix


def plot_ptype(ptype_grid, metadata, i, date_time, dir_gif, categoryNr=4):
    title = 'Precipitation type ' + date_time.strftime("%Y-%m-%d %H:%M")
    fig = plt.figure(figsize=(15, 15))
    #fig.add_subplot(1, 1, 1)
    plot_precipType_field(ptype_grid, geodata=metadata, title=title, colorscale="pysteps", categoryNr=categoryNr)
    #plt.suptitle('Precipitation Type', fontsize=30)
    plt.tight_layout()
    filename = f'{i}.png'
    #  filenames.append(filename)
    plt.savefig(os.path.join(dir_gif, filename), dpi=72)
    plt.close()
    return filename


def calculate_members_mean(membersData):
    """Function to calculate the members average over time

    membersData:
        3D matrix composed by [members, grid dimension 1, grid dimension 2]
    """

    if len(membersData.shape) != 3:
        raise ValueError("Invalid members data shape (expected [:,:,:]) " + str(membersData.shape))

    meanMatrix = np.zeros((membersData.shape[1], membersData.shape[2]))
    for member_idx in range(membersData.shape[0]):
        meanMatrix = meanMatrix + membersData[member_idx, :, :]
    meanMatrix = meanMatrix / membersData.shape[0]
    #print('Mean member matrix done!')

    return meanMatrix


def create_timestamp_indexing(nrOfIncaMessages, startDateTime, timeStep=5, timeBase=60):
    """create a timestamp array for INCA indexing

    nrOfIncaMessages:
        Number of INCA available messages

    startDateTime:
        Start date and time

    timeStep:
        Defines the size of the time step for interpolation

    timeBase:
        Time between messages in minutes (INCA have a message every hour: 60)

    ___
    Return:
          Array of timestamps similar to pysteps timestamps
    """

    if nrOfIncaMessages < 2:
        raise ValueError("Not enough interpolation messages, should be at least 2")

    result = []
    timestamp = startDateTime
    interPoints = np.arange(0, (timeBase + timeStep), timeStep)

    for i in range(nrOfIncaMessages-1):
        for j in interPoints[:-1]:
            result.append(timestamp)
            timestamp = timestamp + datetime.timedelta(minutes=timeStep)

    result.append(timestamp)
    return np.array(result)



def grid_interpolation(numpyGridStart, numpyGridEnd, timeStep=5, timeBase=60):
    """ Time interpolation between 2 2D grids

    numpyGridStart:
        Numpy 2-D grid of start values (INCA)
    numpyGridEnd:
        Numpy 2-D grid of end values (INCA)
    timeStep:
        Size of the time step for interpolation (every 5, 10, 15.. min)
    timeBase:
        Time period considered in minutes (e.g. over one hour = 60, 2 hours = 120)
    applyOver:
        Array with sub-indexes to calculate interpolation (inner grid)
    ----

    Return:
        Returns a list of 3D numpy interpolation matrix
    """
    if numpyGridStart.shape != numpyGridEnd.shape:
        raise ValueError("ERROR: Grids have different dimensions")

    interPoints = np.arange(0, (timeBase + timeStep), timeStep)
    interpolationGrid = np.zeros((len(interPoints), numpyGridStart.shape[0], numpyGridStart.shape[1]))
    interpolationGrid[:, :, :] = np.nan

    #print('Calculating linear interpolation..', end=' ')
    for i in range(len(interPoints)):
        interpolationGrid[i, :, :] = numpyGridStart + ((numpyGridEnd - numpyGridStart) / interPoints[-1]) * interPoints[
            i]
    #print('Done')

    return interpolationGrid



def generate_inca_interpolations(inca_reprojected_data, nwc_timestamps, startdate, timeStep=5, timeBase=60, dateFormat='%Y%m%d%H%M'):
    """Generate a sub-selection of INCA interpolation matrix for all messages available from INCA grib file

    inca_reprojected_data:
        INCA reprojected data.

    inca_timestamps:
        Array of timestamps every timeSteps period, for all grib messages available.

    nwc_timestamps:
        Array of timestamps available from PYSTEPS metadata ['timestamps']

    ----
    Return:
        3D matrix with depth equal to the common matching timestamps between INCA and PYSTEPS.

    """

    # Create a timestamp index array for INCA interpolation matrix
    inca_timestamps = create_timestamp_indexing(inca_reprojected_data.shape[0], startdate, timeStep=timeStep,
                                                timeBase=timeBase)
    # Convert metadata_nwc['timestamps'] to datetime
    nwc_ts = [datetime.datetime.strptime(ts.strftime(dateFormat), dateFormat) for ts in nwc_timestamps]

    inca_start = np.where(inca_timestamps == nwc_ts[0])[0][0]
    inca_end = np.where(inca_timestamps == nwc_ts[-1])[0][0] + 1
    timestamp_selection = inca_timestamps[inca_start:inca_end]  # to be returned

    # interpolation indexes
    resultMatrix = np.zeros((inca_start + len(timestamp_selection), inca_reprojected_data.shape[1], inca_reprojected_data.shape[2]))
    result_idx = 0

    # loop over the messages
    for m in range(1, inca_reprojected_data.shape[0]):
        if result_idx < resultMatrix.shape[0]:
            # calculate interpolations
            interpolationMatrix = grid_interpolation(inca_reprojected_data[m-1], inca_reprojected_data[m], timeStep=timeStep, timeBase=timeBase)
            interp_idx = 0
            # Add the interpolation values to the result matrix (this assignment can be done without looping...)
            while interp_idx < interpolationMatrix.shape[0] and (result_idx < resultMatrix.shape[0]):
                resultMatrix[result_idx, :, :] = interpolationMatrix[interp_idx, :, :]
                result_idx = result_idx + 1
                interp_idx = interp_idx + 1
            result_idx = result_idx - 1  # overwrite the last value

    return resultMatrix[inca_start:], timestamp_selection

# test
# inca_reprojected_data[3,65,65] == resultMatrix[36,65,65]

def calculate_precip_type(incaZnow, incaTemp, incaGroundTemp, precipGrid, topographyGrid, DZML=100., TT0=2., TG0=0.,
                          RRMIN=0):
    """Precipitation type algorithm, returns a 2D matrix with categorical values:
    # PT=0  no precip
    # PT=1  rain
    # PT=2  rain/snow mix
    # PT=3  snow
    # PT=4  freezing rain

    incaZnow:
        INCA snow level 2D grid
    incaTemp:
        INCA temperature 2D grid
    incaGroundTemp:
        INCA ground temperature 2D grid
    precipGrid:
        Precipitation (netCDF PYSTEPS) 2D grid
    topographyGrid:
        Topography grid 2D

    returns:
        2D matrix with categorical data for each type
    """

    # Result grid
    result = np.zeros((precipGrid.shape[0], precipGrid.shape[1]))
    topoZSDiffGrid = (incaZnow - topographyGrid)  # dzs
    precipMask = (precipGrid > RRMIN)

    # SNOW ((dzs<-1.5*DZML) || ( (ZH[i][j] <= 1.5*DZML) && (dzs<=0)))
    snowMask = (topoZSDiffGrid < (-1.5 * DZML)) | ((topographyGrid <= (1.5 * DZML)) & (topoZSDiffGrid <= 0))
    result[snowMask & precipMask] = 3

    # RAIN+SNOW DIAGNOSIS (dzs < 0.5 * DZML) = 2
    rainSnowMask = ~snowMask & (topoZSDiffGrid < (0.5 * DZML))
    result[rainSnowMask & precipMask] = 2

    # RAIN
    rainMask = ~snowMask & ~rainSnowMask
    result[rainMask & precipMask] = 1

    # FREEZING RAIN DIAGNOSIS 4
    # if ((PT[i][j]==1) && ( (tg_<TG0 && TT[i][j]<TT0) || TT[i][j]<TG0))
    freezingMask = (result == 1) & (((incaGroundTemp < TG0) & (incaTemp < TT0)) | (incaTemp < TG0))
    result[freezingMask] = 4

    return result


def get_reprojected_indexes(reprojectedGrid):
    """Reprojected INCA grids contains a frame of NAN values, this function returns the start and end indexes
    of the inca grid over the reprojected grid

    reprojectedGrid:
        INCA reprojected Grid

    ---
    Returns:
        x y indexes of inca reprojected grid over pysteps dimensions
    """

    x_start = np.where(~np.isnan(reprojectedGrid))[0][0]
    x_end = np.where(~np.isnan(reprojectedGrid))[0][-1] + 1
    y_start = np.where(~np.isnan(reprojectedGrid))[-1][0]
    y_end = np.where(~np.isnan(reprojectedGrid))[-1][-1] + 1

    return x_start, x_end, y_start, y_end


def precipitation_type_calculator(filename,
                                  startdate,
                                  gribPath,
                                  topoFilename,
                                  keys=None,
                                  inca_projectionString=None,
                                  nwc_projectionString=None,
                                  members=None,
                                  timeBase=None,
                                  timeStep=None,
                                  projection_bounds=None,
                                  topo_interpolation=None,
                                  desired_output=None):

    # Run checks to ensure correct input parameters
    if not isinstance(filename, str):
        raise TypeError("Filename must be a string containing the path to the netCDF file")
    if filename.rsplit('.', 1)[1] != 'nc':
        raise ValueError("Filename must be a netCDF file ending in '.nc'")

    if not isinstance(gribPath, str):
        raise TypeError("gribPath must be a string containing the path to the directory which contains the grib files")

    if not isinstance(topoFilename, str):
        raise TypeError("topoFilename must be a string containing the path to the topography file")

    if keys is not None and not isinstance(keys, (list, tuple)):
        raise TypeError("keys must be a list or tuple")

    if inca_projectionString is not None and not isinstance(inca_projectionString, str):
        raise TypeError("inca_projectionString must be a string containing the INCA projection string")

    if nwc_projectionString is not None and not isinstance(nwc_projectionString, str):
        raise TypeError("nwc_projectionString must be a string containing the NWC projection string")

    if members is not None and not isinstance(members, int):
        raise TypeError("members must be an integer")
    if not members > 0 or not members < 27:
        raise ValueError("members must be a positive integer between 1 and 26 inclusive")

    if timeBase is not None and not isinstance(timeBase, int):
        raise TypeError("timeBase must be a positive integer")
    if not timeBase > 0:
        raise ValueError("timeBase must be a positive integer")

    if timeStep is not None and not isinstance(timeStep, int):
        raise TypeError("timeStep must be a positive integer")
    if not timeStep > 0:
        raise ValueError("timeStep must be a positive integer")

    if projection_bounds is not None:
        if not isinstance(projection_bounds, (list, tuple)):
            raise TypeError("projection_bounds must be a list or tuple")
        if not len(projection_bounds) == 4:
            raise ValueError("projection_bounds must be a list or tuple of length 4")
        if not all(isinstance(x, (int, float)) for x in projection_bounds):
            raise ValueError("projection_bounds must be a list or tuple containing 4 ints or floats")

    if topo_interpolation is not None and not isinstance(topo_interpolation, bool):
        raise TypeError("topo_interpolation must be a boolean")

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

    keys : list
        A list of the names of the GRIB import selection keys.

    inca_projectionString : str
        The INCA projection string.

    nwc_projectionString : str
        The nowcast projection string.

    members : int
        The number of ensemble members which you would like to receive the results of. This number can be between 1 and 26.

    timeBase : int
        The base for the time period. (min)

    timeStep : int
        The step for the time period. (min)

    projection_bounds = list
        A list containing the four bounds of the INCA projection. For use in the metadata.

    topo_interpolation = boolean
        Whether the topography requires interpolation. Default no, but can be adjusted to yes.
        
    desired_output = str
        The desired output is an indicator used by the user to indicate their desired output parameter. 
        Currently the function only features the optional outputs of full arrays of the number of members. 
        This option is chosen by the string 'members'.
        Or the mean precipitation type across all members.
        This option is chosen by the string 'mean',
        
    Returns
    -------
    One or more 4D arrays containing the desired output of the user. 
    Output arrays take the form [member, timeStamp, X-coord, Y-coord].
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
    if projection_bounds is None:
        projection_bounds = [360000, 350000, 960000, 940000]
    if topo_interpolation is None:
        topo_interpolation = False
    if desired_output is None:
        desired_output = 'members'

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
    metadata_inca['x1'] = projection_bounds[0]
    metadata_inca['y1'] = projection_bounds[1]
    metadata_inca['x2'] = projection_bounds[2]
    metadata_inca['y2'] = projection_bounds[3]
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
    metadata_nwc['cartesian_unit'] = metadata_nwc['projection'][
                                     nwc_projectionString.find('units=') + 6:
                                     nwc_projectionString.find(
                                         ' +no_defs', nwc_projectionString.find(
                                             'units=') + 6)]
    print('netCDF4 load done')

    # ---------------------------------------------------------------------------
    # measure time
    start_time = time.time()

    # --------------------------------------------------------------------------
    # Reproject
    #     projection over pySTEPS grid
    R_inca_ZS, _ = reproject_grids(R_inca_ZS, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
    R_inca_TT, _ = reproject_grids(R_inca_TT, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
    R_inca_TG, _ = reproject_grids(R_inca_TG, r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
    # The topography file is not required to be interpolated if it is already of the shape matching the grib files,
    # i.e (591, 601)
    if topo_interpolation is True:
        topo_grid, _ = reproject_grids(np.array([topo_grid]), r_nwc[0, 0, :, :], metadata_inca, metadata_nwc)
    print('Re-projection done')

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
        for member in range(0, members):
            res = np.copy(ptype_mean)
            res[r_nwc[member, ts, x1:x2, y1:y2] == 0] = 0
            ptype_list[member, ts, :, :] = res

        # Add mean result at the end
        ptype_list[-1, ts, :, :] = ptype_mean

    if desired_output == 'members':
        output = ptype_list[:members]
    else:
        output = ptype_list[-1]

    print("--Script finished--")
    print("--- %s seconds ---" % (time.time() - start_time))

    return output
