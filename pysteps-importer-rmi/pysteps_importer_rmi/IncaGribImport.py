import eccodes as ec
import numpy as np


class IncaGribReader(object):
    """Base class handling the eccodes module and storing the message index.
    Because eccodes is a stateful module (meaning that states and data are stored inside as class
    variables), the present class is therefore also a stateful class.

    Warning
    -------
    !!! Beware !!!
    Since eccodes is a stateful module (sadly), only one instance of MessagesCollection can
    exist at the same time in the same process.
    Two concurrent instances of this class will result in collisions and ultimately lead to inconsistencies.

    Attributes
    ----------
    filename: str
        Path to the GRIB file to consider.
    file: File Object
        Pointer to the file object to access the GRIB file.
    read_mode: str
        Read mode of the GRIB file. Default is 'rb'.
    gid: list(int)
        List of the messages index.
    currentMsgId: int
        Id of the current GRIB message
    """

    def __init__(self):

        self.filename = ''
        self.file = None
        self.read_mode = ''
        self.gid = list()
        self.currentMsgId = None

    def __len__(self):
        """Number of messages loaded.

        Returns
        -------
        int
            The number of messages loaded.
        """
        return len(self.gid)

    @property
    def not_loaded(self):
        """Number of messages not loaded.

        Returns
        -------
        int
            The number of messages not loaded.
        """
        if self.file is None:
            return None

        return ec.codes_count_in_file(self.file)

    def open(self, filename, mode='rb'):
        """Open a GRIB file to be read.

        Parameters
        ----------
        filename: str
            The path to the GRIB file to consider.
        mode: str
            Read mode of the GRIB file. Default is binary read-only 'rb'.
        """

        self.clear()

        self.read_mode = mode
        self.filename = filename
        try:
            self.file = open(filename, mode)
        except BaseException as err:
            print(f"Error {err=}, {type(err)=}")
            self.clear()
            raise

    def clear_messages(self):
        """Clear all the message loaded."""

        for id in self.gid:
            ec.codes_release(id)

        self.gid.clear()

    def clear(self):
        """Clear all the message loaded and the stored variables name."""

        self.clear_messages()

    def close(self):
        """Close the opened GRIB file (if any)."""

        self.clear()
        self.read_mode = ''
        self.filename = ''
        if self.file is not None:
            self.file.close()
            self.file = None

    def load_message(self, headers_only=False):
        """Load the next not-loaded message of the GRIB file (sequential loading).

        Parameters
        ----------
        headers_only: bool
            Indicate to load only the headers of the message.
            Since the use of this feature is not recommended for the moment in eccodes,
            the default is False and should stays so until further notice.
            Use at your own risks.
        """
        gid = ec.codes_grib_new_from_file(self.file, headers_only=headers_only)
        if gid is not None:
            self.gid.append(gid)
        self.set_current_message()
        return gid

    def load_messages(self, select_dict=None):
        """Load the messages of the GRIB file according to a dictionary of key values.

        Warning
        -------
        Not providing a dictionary of key values may lead to a huge memory usage and a possible
        crash of the machine running the code.

        Parameters
        ----------
        select_dict: dict
            Dictionary of GRIB key and key values to select message. Default is None, which
            will try to load all the message of the GRIB file, possibly leading to a huge memory
            usage.
        """

        if select_dict is None:
            while True:
                gid = self.load_message()
                if gid is None:
                    break
        else:
            keys = list(select_dict.keys())
            iid = ec.codes_index_new_from_file(self.filename, keys)

            for key in keys:
                ec.codes_index_select(iid, key, select_dict[key])

            while True:
                gid = ec.codes_new_from_index(iid)

                if gid is not None:
                    self.gid.append(gid)
                else:
                    break
            ec.codes_index_release(iid)

        self.set_current_message()

    def release_message(self, id):
        """Release a particular already loaded message.

         Parameters
         ----------
         id: int
            Id (position) of the message in the list of loaded messages.
        """

        ec.codes_release(self.gid[id])
        del self.gid[id]

    def get_key_values(self, key):
        """Get the values of a particular key across the sequential list of loaded
        messages.

        Parameters
        ----------
        key: str
            Key to get the values from.
        """

        values = list()

        for id in self.gid:
            values.append(ec.codes_get(id, key))

        return values

    def get_key_value_index(self, key, index):
        """Get the key value of a particular loaded message.

        Parameters
        ----------
        key: str
            Key to get the values from.
        index: int
            Index of the message to get the value from (list)
        """

        return ec.codes_get(self.gid[index], key)

    def get_message_grid_value_index(self, index):
        """Return the grid values of a particular loaded message (1-dim vector).

        Parameters
        ----------
        index:
            index of the message.
        """
        if len(self.gid) == 0:
            return []

        return ec.codes_get_values(self.gid[index])

    def get_message_grid_values(self):
        """Return the grid values of all loaded messages (1-dim vector).
        """

        values = list()
        if len(self.gid) == 0:
            return []

        for idx in range(len(self.gid)):
            values.append(self.get_message_grid_value_index(idx))

        return values

    def get_current_message_attributes(self, keys):
        """Return the values of the selected keys for the current message.

        Parameters
        ----------
        keys: list
            keys to read from the current message.
        """

        if self.currentMsgId is None:
            return None
        else:
            attributes = {}
            for k in keys:
                try:
                    keyValue = ec.codes_get(self.currentMsgId, k)
                    attributes[k] = keyValue
                except eccodes.KeyValueNotFoundError as err:
                    print('Key not found', k)
            return attributes

    def get_current_message_grid(self):
        """Return the grid values of the current loaded message (1-dim vector).
        """

        if self.currentMsgId is None:
            return None
        else:
            return ec.codes_get_values(self.currentMsgId)

    def set_current_message(self):
        """Set the first gid as the current message to read.
        """

        if len(self.gid) > 0:
            self.currentMsgId = self.gid[0]
        else:
            self.currentMsgId = None

    def next_message(self):
        """Set next message id.
        """
        self.release_message(0)
        self.set_current_message()


# ----------------------------------------------------------------

class IncaGribImporter:
    """Importer class: This class will contain the logic to select and format the required data from an INCA GRIB file.

    Attributes
    ----------
    reader: IncaGribReader
        Instance of the reader class to access the GRIB data.
    """

    def __init__(self):
        self.reader = IncaGribReader()

    def load_grib(self, filename):
        try:
            self.reader.open(filename)
            self.reader.load_messages()
            return True
        except FileNotFoundError as err:
            self.__close_grib()
            print('Can not load the file')
            return False

    def __close_grib(self):
        self.reader.clear()
        self.reader.close()

    def __retrieve_projection(self, projection_keys):
        """ Read projection data from the first message (this information is the same in all messages).

        projection_keys: list
            List of key values to select from the first (current) message
        """

        return self.reader.get_current_message_attributes(projection_keys)

    def __retrieve_messages_data(self, message_keys, grid_dimensions):
        """ Read messages sequentially and returns a selection of attributes (keys) and the grid from each message.

        Using Nx and Ny as optional reshape parameters (?)

        message_keys: list()
            List of key values to select from each message
        grid_dimensions: list()
            List of dimensions to reshape the grid vector (if none, returns 1-d array)
        """

        messageData = {}
        idx = 1
        while self.reader.currentMsgId is not None:
            messageData[idx] = {}
            if message_keys is not None:  # Should be fine, there are usually not many messages
                messageAttributes = self.reader.get_current_message_attributes(message_keys)
                messageData[idx]['Attributes'] = messageAttributes
            else:
                messageData[idx]['Attributes'] = {}
            messageArray = self.reader.get_current_message_grid()
            if grid_dimensions is None:
                messageData[idx]['Grid'] = messageArray
            else:
                try:
                    messageMatrix = np.array(messageArray).reshape(grid_dimensions['Nx'], grid_dimensions['Ny'])
                except ValueError as err:
                    print(err)
                    return {}
                messageMatrix = messageMatrix[::-1, :]  # Reorientation
                messageData[idx]['Grid'] = messageMatrix
            self.reader.next_message()
            idx = idx + 1

        return messageData

    def retrieve_grib_data(self, filename, metadata_keys=None, message_keys=None, getGrid=True,
                           grid_dimensions={'Nx': 591, 'Ny': 601 }):

        gribDictionary = {}

        # Read the file
        if not self.load_grib(filename):
            print('Execution stopped')
        else:
            # read and add projection data
            if metadata_keys is not None:
                gribDictionary['Metadata'] = self.__retrieve_projection(metadata_keys)
            else:
                gribDictionary['Metadata'] = {}

                # get the message grid data
            if getGrid:
                gribDictionary['Messages'] = self.__retrieve_messages_data(message_keys, grid_dimensions)
            else:
                gribDictionary['Messages'] = {}

                # close the file
            self.__close_grib()

        return gribDictionary


