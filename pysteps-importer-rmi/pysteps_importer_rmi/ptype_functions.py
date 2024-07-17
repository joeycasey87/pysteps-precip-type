# Import
import os
import datetime
import numpy as np
import matplotlib.pyplot as plt
from precipitationTypeFields import plot_precipType_field


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

