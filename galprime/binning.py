""" Binning methods

This module contains all necessary methods for binning catalogues according to various parameters.

"""

from astropy.table import Table
from numpy import array, copy, append, round, reshape


class Bin:
    """ Class for a single bin of information

    Attributes:
        objects: The catalog rows that belong to the given bin
        object_column_names: The column names when the objects are sorted column-wise instead of row-wise.
        bin_params (arr): The values that define the bounds of the bin.
        bin_param_names (arr): Array of parameter names that the bin was made with.
    """

    def __init__(self, objects=None, object_column_names=None, bin_params=None, bin_param_names=None):
        if bin_params is None:
            bin_params = []
        if bin_param_names is None:
            bin_param_names = []

        self.objects = objects
        self.object_column_names = object_column_names
        self.bin_params = bin_params
        self.bin_param_names = bin_param_names

    def rebin(self, param_index, low_bounds, bin_width, number_threshold=4):
        """ Returns a list of Bin objects after re-binning based on a certain value and bin parameters.

        Bins specified by the lower bounds and bin width.

        Args:
            param_index: Which index in the objects contains the desired value.
            low_bounds: Lower bounds on the bins.
            bin_width: The width of the bins.
            number_threshold: The required number of objects a bin needs to have to be saved to the list.

        Returns:
             List of re-binned Bin objects.
        """
        new_bins = []
        top_level_params = self.bin_params
        top_level_param_names = self.bin_param_names

        for low in low_bounds:
            this_bin = []

            for obj in self.objects:
                if low < obj[param_index] < low + bin_width:
                    this_bin.append(obj)
            if len(this_bin) <= number_threshold:
                continue

            new_bin = Bin(objects=this_bin, object_column_names=self.object_column_names,
                          bin_params=append(copy(top_level_params), (round(low, 2), round(low + bin_width, 2))),
                          bin_param_names=append(copy(top_level_param_names), self.key_from_index(param_index)))
            new_bins.append(new_bin)

        return new_bins

    def return_columns(self):
        """
        Returns a dictionary of columns (turns the object list back into a column format)
        """
        keys = copy(self.object_column_names)
        columns = []
        for i in range(0, len(keys)):
            columns.append([])
            for j in range(0, len(self.objects)):
                columns[i].append(self.objects[j][i])

        return list(keys), columns

    def index_from_key(self, key=""):
        """ Return the index of the column name key specified."""
        return self.object_column_names.index(key)

    def key_from_index(self, index=0):
        """ Returns the column name at a given index. """
        return self.object_column_names[index]

    def values_in_bin(self, values):
        """ Returns a boolean whether of not a set of values is within the bin.

        Args:
            values: A list of values.

        Returns:
            True if all values are in the bin. False otherwise.
        """
        reshaped_values = reshape(self.bin_params, newshape=(int(len(self.bin_params) / 2), 2))
        good_count = 0
        for i in range(0, len(values)):
            if reshaped_values[i][0] < values[i] < reshaped_values[i][1]:
                good_count += 1

        if good_count == len(values):
            return True
        else:
            return False

    def return_param_dict(self):
        """
        Returns a dict object where the keys are the column names that built the parameter set,
        and the values are the corresponding parameters.
        """
        param_dict = {}
        for i in range(0, len(self.bin_param_names)):
            param_dict[self.bin_param_names[i]] = (self.bin_params[2 * i], self.bin_params[2 * i + 1])
        return param_dict


def select_bin(bin_list, values):
    """ Select a bin from a list of bins based on an array of values.

    Args:
        bin_list (arr): List of Bin objects
        values (arr): Array of parameters.
            Not that the parameters need to be in the same order that the bins were generated from.

    Returns:
        Bin corresponding to the appropriate parameters. If no bin exists for those parameters, returns None.

    """
    for b in bin_list:
        if b.values_in_bin(values):
            return b
    return None


def generate_objects(columns):
    """ Reorganizes the table columns into an object-format
    Args:
        columns: A list of table columns

    Returns:
        List of objects.
    """
    objects = []
    for i in range(0, len(columns[0])):
        this_obj = []
        for j in range(0, len(columns)):
            this_obj.append(columns[j][i])
        objects.append(this_obj)
    return objects


def bin_catalog(config):
    """ Bin the inputted catalog.

    This is a convenience function to get a list of bins based on redshift, star-formation, and mass.

    Args:
        config: Values from inputted config file or user inputs.

    Returns:
         List of bins binned according to the config settings.
    """

    # Read HST-ZEST catalog into memory
    catalog = Table.read(config["CATALOG"], format="fits")

    mags = array(catalog[config["MAG_KEY"]])
    r50s = array(catalog[config["R50_KEY"]])
    ns = array(catalog[config["N_KEY"]])
    ellips = array(catalog[config["ELLIP_KEY"]])
    masses = array(catalog[config["MASS_KEY"]])
    sfprobs = array(catalog[config["SFPROB_KEY"]])
    redshifts = array(catalog[config["Z_KEY"]])

    mass_bins, mass_step = config["MASS_BINS"], config["MASS_STEP"]
    redshift_bins, redshift_step = config["REDSHIFT_BINS"], config["REDSHIFT_STEP"]
    sfprob_bins, sfprob_step = config["SFPROB_BINS"], config["SFPROB_STEP"]

    catalog_objects = generate_objects((mags, r50s, ns, ellips, masses, redshifts, sfprobs))
    column_names = ["MAGS", "R50S", "NS", "ELLIPS", "MASSES", "REDSHIFTS", "SFPROBS"]
    init_bin = Bin(objects=catalog_objects, object_column_names=column_names, bin_params=[])

    mass_bins = init_bin.rebin(4, mass_bins, mass_step)
    mz_bins, full_bins = [], []
    for cat_bin in mass_bins:
        mz_bins.extend(cat_bin.rebin(5, redshift_bins, redshift_step))
    for cat_bin in mz_bins:
        full_bins.extend(cat_bin.rebin(6, sfprob_bins, sfprob_step))

    return full_bins


def bin_mag_catalog(mag_table, b, mag_table_keys=None, bin_keys=None):
    """ Bin an external catalogue based on a bin's provided parameters.

    Args:
        mag_table The magnitude table to rebin
        b The tbridge.Bin object
        mag_table_keys The magnitude table keys to use for binning.
        bin_keys The associated bin keys (must be the same length AND order as mag_table_keys).

    Returns:
        Table of magnitudes binned by the appropriate parameters.
    """
    if bin_keys is None:
        bin_keys = ["MASSES", "REDSHIFTS"]
    if mag_table_keys is None:
        mag_table_keys = ["MASS", "Z"]

    # Gather the necessary information from the provided bin.
    bin_dict = b.return_param_dict()

    # Iterate over all keys and only return the rows that are within the bin limits.
    for i in range(len(bin_keys)):
        low, high = bin_dict[bin_keys[i]]
        mag_table = mag_table[mag_table[mag_table_keys[i]] >= low]
        mag_table = mag_table[mag_table[mag_table_keys[i]] <= high]

    return mag_table


def bin_index(val, bins):
    """ Get bin index for a given value and a set of bin parameters """
    for index in range(0, len(bins) - 1):
        if bins[index] < val < bins[index + 1]:
            return index
    return len(bins) - 1
