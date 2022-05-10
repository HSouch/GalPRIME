"""
This module contains all general-purpose functions for loading and saving data. This includes pulling data out of
FITS files, saving TBriDGE simulation outputs, and more.
"""

import galprime
import os
import math
from pathlib import Path

from astropy.io import fits
from astropy.wcs import wcs
from astropy.table import Table
from astropy.nddata import Cutout2D, NoOverlapError
from astropy.stats import sigma_clipped_stats

from numpy import sqrt, str, save, load, ceil, zeros, isnan, floor, ceil, sort, min, mean
from numpy.random import choice, randint, uniform


def get_closest_psf(psfs, obj_ra, obj_dec):
    """ Get the closest psf for a given object's RA and DEC"""

    def dist(ra1, ra2, dec1, dec2):
        return sqrt((ra1 - ra2) ** 2 + (dec1 - dec2) ** 2)

    shortest_dist = 999999
    closest_psf = None

    for psf in psfs:
        head = psf.header
        psf_ra, psf_dec = head["RA"], head["DEC"]
        local_dist = dist(obj_ra, psf_ra, obj_dec, psf_dec)
        if local_dist < shortest_dist:
            shortest_dist = local_dist
            closest_psf = psf

    return closest_psf


def get_wcs(fits_filename):
    """
        Finds and returns the WCS for an image. If Primary Header WCS no good, searches each index until a good one
        is found. If none found, raises a ValueError
    """
    # Try just opening the initial header
    wcs_init = wcs.WCS(fits_filename)
    ra, dec = wcs_init.axis_type_names
    if ra.upper() == "RA" and dec.upper() == "DEC":
        return wcs_init

    else:
        hdu_list = fits.open(fits_filename)
        for n in hdu_list:
            try:
                wcs_slice = wcs.WCS(n.header)
                ra, dec = wcs_slice.axis_type_names
                if ra.upper() == "RA" and dec.upper() == "DEC":
                    return wcs_slice
            except:
                continue
        hdu_list.close()

    raise ValueError


def get_arc_conv(w: wcs.WCS):
    """ Gets pixels to arc-seconds conversion scale (Number of arcseconds per pixel) """
    pix_x, pix_y = 1, 1
    ra_1, dec_1 = w.wcs_pix2world(pix_x, pix_y, 0)
    ra_2, dec_2 = w.wcs_pix2world(pix_x + 1, pix_y + 1, 0)
    diff_1 = abs(ra_1 - ra_2) * 3600
    diff_2 = abs(dec_1 - dec_2) * 3600
    return (diff_1 + diff_2) / 2


def get_angular_size_dist(z, H0=71, WM=0.27, WV=None):
    """
    Return the angular size distance in Megaparsecs.
    (Stripped down version of Cosmocalc by Ned Wright and Tom Aldcroft (aldcroft@head.cfa.harvard.edu))
    """
    try:
        c = 299792.458  # velocity of light in km/sec

        if z > 100:
            z /= 299792.458  # Values over 100 are in km/s

        if WV is None:
            WV = 1.0 - WM - 0.4165 / (H0 * H0)  # Omega(vacuum) or lambda
        age = 0.0  # age of Universe in units of 1/H0

        h = H0 / 100.
        WR = 4.165E-5 / (h * h)  # includes 3 massless neutrino species, T0 = 2.72528
        WK = 1 - WM - WR - WV
        az = 1.0 / (1 + 1.0 * z)
        n = 1000  # number of points in integrals
        for i in range(n):
            a = az * (i + 0.5) / n
            adot = math.sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
            age += 1. / adot

        DCMR = 0.0

        # do integral over a=1/(1+z) from az to 1 in n steps, midpoint rule
        for i in range(n):
            a = az + (1 - az) * (i + 0.5) / n
            adot = math.sqrt(WK + (WM / a) + (WR / (a * a)) + (WV * a * a))
            DCMR = DCMR + 1. / (a * adot)

        DCMR = (1. - az) * DCMR / n

        # tangential comoving distance
        ratio = 1.0
        x = math.sqrt(abs(WK)) * DCMR
        if x > 0.1:
            if WK > 0:
                ratio = 0.5 * (math.exp(x) - math.exp(-x)) / x
            else:
                ratio = math.math.sin(x) / x
        else:
            y = x * x
            if WK < 0:
                y = -y
            ratio = 1. + y / 6. + y * y / 120.
        DCMT = ratio * DCMR
        DA = az * DCMT
        Mpc = lambda x: c / H0 * x
        DA_Mpc = Mpc(DA)

        return DA_Mpc
    except:
        raise ValueError


def generate_cutout(image, position, img_wcs=None, size=91, world_coords=True):
    """
    Generates a cutout for a given image. Uses world coordinates by default, but can be configured to take
    in a position corresponding to the actual array indices.
    :param image:
    :param position:
    :param wcs:
    :param size:
    :param world_coords:
    :return:
    """

    if size < 71:
        size = 71
    elif size % 2 == 0:
        size += 1

    if world_coords:
        coord = img_wcs.wcs_world2pix(position[0], position[1], 0, ra_dec_order=True)
        pix_x, pix_y = coord[0], coord[1]
    else:
        pix_x, pix_y = position[0], position[1]

    # Send an index error iff the position does not lie in the image.
    if 0 < pix_x < image.shape[0] and 0 < pix_y < image.shape[1]:
        try:
            cut = Cutout2D(image, (pix_x, pix_y), size)
            return cut.data
        except NoOverlapError:
            raise IndexError
    else:
        raise IndexError


def trim_objects(obj_catalog, ras, decs, mag_lim=None, ra_key="RA", dec_key="DEC", mag_key="i"):
    this_catalog = obj_catalog[obj_catalog[ra_key] > ras[0]]
    this_catalog = this_catalog[this_catalog[ra_key] < ras[1]]

    this_catalog = this_catalog[this_catalog[dec_key] > decs[0]]
    this_catalog = this_catalog[this_catalog[dec_key] < decs[1]]

    if mag_lim is not None:
        this_catalog = this_catalog[this_catalog[mag_key] > mag_lim[0]]
        this_catalog = this_catalog[this_catalog[mag_key] < mag_lim[1]]

    return this_catalog


def extraction_limits(image_wcs, image_shape, cutout_size):
    """ Get the minimum and maximum RA and DEC limits for a given cutout size and image"""
    half_width = ceil(cutout_size / 2)

    pix_min = image_wcs.wcs_pix2world(half_width, half_width, 0)
    pix_max = image_wcs.wcs_pix2world(image_shape[0] - half_width, image_shape[1] - half_width, 0)

    ras = sort((pix_min[0], pix_max[0]))
    decs = sort((pix_min[1], pix_max[1]))

    return ras, decs


def get_image_filenames(images_directory, image_band="i", check_band=False):
    """
    Retrieves a list of all available filenames for a given directory, and a given band.
    WARNING: Optimized for HSC filename format (ex: HSC-I_9813_4c3.fits).
    """
    image_filenames = []
    images = Path(images_directory).rglob('*.fits')
    for image in images:
        image = str(image)
        if check_band:
            image_no_path = image.split("/")[len(image.split("/")) - 1]
            filename_band = image_no_path.split("_")[0].split("-")[1].lower()
            if filename_band == image_band:
                image_filenames.append(image)
        else:
            image_filenames.append(image)
    return image_filenames


def select_image(filename):
    """ Opens an image-based FITS file and returns the first available image.

    Args:
        filename (str): The FITS filename and path.

    Returns:
        numpy.ndarray: The image obtained from the FITS file. If none found, will return None.

    """
    with fits.open(filename) as HDUList:
        image = None
        for i in range(0, len(HDUList)):
            try:
                image = HDUList[i].data
                if image.shape[0] > 0 and image.shape[1] > 0:
                    break
            except:
                continue
    return image


def generate_file_prefix(bin_params):
    """ Use the bin params to generate a file prefix."""
    prefix = "bin_"
    for j in range(0, len(bin_params)):
        if (j + 1) % 2 != 0:
            prefix += str(bin_params[j]) + "-"
        else:
            prefix += str(bin_params[j]) + "_"
    return prefix


def save_profiles(profile_list, bin_info, out_dir, keys, bg_info=None, cutout_info=None, structural_params=None):
    """ Saves a set of profiles into a properly formatted output directory, with proper filename format.

    Args:
        profile_list: The list of profiles (shape is m x n, where m is the number of different models for each
            object and n is the number of objects i.e. m rows of profiles from n objects.
        bin_info (arr): bin information for profile formatting (Obtained using the bin_params attribute from
            the tbridge.binning.Bin object).
        out_dir (str): the output directory to save the files to.
        keys (arr): the keys to generate subdirectory and file names with.
        bg_info (numpy.ndarray): Array of associated background info for the extracted profiles. Saves to its own
                    subdirectory.
        structural_params (astropy.table.Table): Astropy table of associated structural parameters
    """

    # Generate output structure
    subdirs = []
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    for key in keys:
        subdir = key + "_profiles/"
        if not os.path.isdir(out_dir + subdir):
            os.mkdir(out_dir + subdir)
        subdirs.append(subdir)

    filename_prefix = generate_file_prefix(bin_info)
    valid_colnames = ["sma", "intens", "intens_err", "ellipticity", "ellipticity_err", "pa", "pa_err"]

    # Save profiles as FITS HDULists in the proper directories
    for i in range(0, len(profile_list)):
        profiles = profile_list[i]
        out_filename = filename_prefix + keys[i] + ".fits"
        out_hdulist = fits.HDUList()
        for prof in profiles:
            try:
                out_hdulist.append(fits.BinTableHDU(Table([prof[col] for col in valid_colnames],
                                                          names=valid_colnames)))
            except Exception as e:
                print(e)
                continue

        out_hdulist.writeto(out_dir + subdirs[i] + out_filename, overwrite=True)

    # Save the additional available information
    if bg_info is not None:
        if not os.path.isdir(out_dir + "bg_info/"):
            os.mkdir(out_dir + "bg_info/")
        galprime.save_array(bg_info, out_dir + "bg_info/" + filename_prefix + "bgs.npy")
    if cutout_info is not None:
        if not os.path.isdir(out_dir + "cutout_info/"):
            os.mkdir(out_dir + "cutout_info/")
        galprime.save_array(cutout_info, out_dir + "cutout_info/" + filename_prefix + "bgs.npy")

    if structural_params is not None:
        if not os.path.isdir(out_dir + "params/"):
            os.mkdir(out_dir + "params/")
        structural_params.write(out_dir + "params/" + galprime.generate_file_prefix(bin_info) + "params.fits",
                                format="fits", overwrite=True)


def save_profile_set(profiles, out_filename="profiles.fits"):
    """ Save an array of profiles to a FITS file.

    Args:
        profiles: List of profile tables.
        out_filename: Filename to save FITS file to.

    """

    valid_colnames = ["sma", "intens", "intens_err", "ellipticity", "ellipticity_err", "pa", "pa_err"]
    out_hdulist = fits.HDUList()
    for prof in profiles:
        out_hdulist.append(fits.BinTableHDU(Table([prof[col] for col in valid_colnames],
                               names=valid_colnames)))

    out_hdulist.writeto(out_filename, overwrite=True)


def load_profile_set(filename):
    """ Load a set of profiles from a FITS file specified by the filename"""
    tables = []
    with fits.open(filename) as HDUList:
        for hdu in HDUList:
            try:
                tables.append(Table.read(hdu))
            except:
                continue
    return tables


def save_cutouts(cutouts, output_filename="cutouts.fits"):
    """ Save a set of cutouts to a fits HDUList object """
    out_hdulist = fits.HDUList()

    for cutout in cutouts:
        out_hdulist.append(fits.ImageHDU(data=cutout))

    out_hdulist.writeto(output_filename, overwrite=True)


def load_cutouts(filename):
    """ Load in a list of cutouts from the given filename (will try all hdus in the HDUList)"""
    cutouts = []
    HDUList = fits.open(filename)
    for n in HDUList:
        try:
            image = n.data
            if image.shape[0] > 0 and image.shape[1] > 0:
                cutouts.append(image)
        except:
            continue

    HDUList.close()

    return cutouts


def tables_from_file(filename):
    """
    Load a set of tables from a given HDUList.
    :param filename: filename to gather tables from.
    :return:
    """
    tables = []

    HDUList = fits.open(filename)
    for n in HDUList:
        try:
            t = Table.read(n)
            tables.append(t)
        except:
            continue
    HDUList.close()

    return tables


def params_from_filename(filename):
    """ Return the parameters inferred from the input filename (based on TBRIDGE formatting)."""
    filename_no_path = filename.split("/")[len(filename.split("/")) - 1]
    splits = filename_no_path.split("_")
    splits = splits[1:len(splits) - 1]
    params = []
    for split in splits:
        low_high = split.split("-")
        params.append((float(low_high[0]) + float(low_high[1])) / 2)
    return params


def save_array(arr, filename="arr.npy"):
    """ Save an array (npy binary format)"""
    save(filename, arr)


def load_array(filename):
    """ Load an array """
    try:
        return load(filename)
    except FileNotFoundError:
        print("File: ", filename, " not found. Returning NoneType")
        return None


def get_backgrounds(config, n=50, size=None, threshold=1e-4, return_psfs=True, return_bg_info=True):
    """ Retrieve a random set of backgrounds from the input image directory.
    Args:
        config_values: User provided configuration file.
        n: Number of backgrounds to generate.
        return_psfs (bool): Return the psfs associated with the backgrounds
        return_bg_info (bool): Return the background info.
    :return:
    """

    image_dir, psf_filename = config["IMAGE_DIRECTORY"], config["PSF_FILENAME"]
    size = config["SIZE"]

    with fits.open(psf_filename) as psfs:
        image_filenames = galprime.get_image_filenames(image_dir)

        psf_list, bg_added_models, image_list, ras, decs = [], [], [], [], []

        bg_infotable = {"IMAGES": [], "RAS": [], "DECS": [], "XS": [], "YS": []}

        for i in range(0, n):
            fail_count = 0
            while fail_count < 10:
                model_halfwidth = ceil(size / 2)

                image_filename = choice(image_filenames)

                image = galprime.select_image(image_filename)
                image_wcs = galprime.get_wcs(image_filename)

                c_x = randint(model_halfwidth + 1, image.shape[0] - model_halfwidth - 1)
                c_y = randint(model_halfwidth + 1, image.shape[1] - model_halfwidth - 1)
                x_min, x_max = int(c_x - model_halfwidth), int(c_x + model_halfwidth)
                y_min, y_max = int(c_y - model_halfwidth), int(c_y + model_halfwidth)

                image_cutout = image[x_min: x_max - 1, y_min: y_max - 1]
                ra, dec = image_wcs.wcs_pix2world(c_x, c_y, 0)

                if threshold is not None:
                    bg_mean, bg_median, bg_std = sigma_clipped_stats(image_cutout, sigma=3.)
                    if bg_std < threshold:
                        fail_count += 1
                        continue
                    else:
                        break

                if image is None:
                    fail_count += 1
                    continue

            if fail_count >= 10:
                continue

            if return_psfs:
                psf = galprime.get_closest_psf(psfs, ra, dec).data
                psf_list.append(psf)

            bg_added_models.append(image_cutout)
            bg_infotable["IMAGES"].append(image_filename)
            bg_infotable["RAS"].append(ra)
            bg_infotable["DECS"].append(dec)
            bg_infotable["XS"].append(c_x)
            bg_infotable["YS"].append(c_y)

    if return_psfs:
        if return_bg_info:
            return bg_added_models, psf_list, bg_infotable
        else:
            return bg_added_models, psf_list
    elif return_bg_info:
        return bg_added_models, bg_infotable
    else:
        return bg_added_models


def as_dir(directory):
    """ Add a forward slash if one is not at the end of a string. """
    if directory[-1] != '/':
        return directory + "/"
    else:
        return directory


def cutout_stitch(cutouts, masked_cutouts=None, output_filename=None):
    """ Generate a full stitch from a set of cutouts.

    Args:
        cutouts: Array of cutouts (ideally unmasked but it makes no difference). Need to be all the same size as the
            first cutout.
        masked_cutouts: The same cutouts but masked
        output_filename: Optional output to save as HDUList

    Returns:
        Returns either the stitch, or a tuple of the stich and the mask stitch if masked_cutouts is submitted..
    """

    i = int(sqrt(len(cutouts)))
    j = int(len(cutouts) / i)

    c_width, c_height = cutouts[0].shape[0], cutouts[1].shape[1]
    w, h = i * c_width, j * c_height

    canvas = zeros((w, h))

    if masked_cutouts is not None:
        masked_canvas = zeros((w, h))
    else:
        masked_canvas = None

    n = 0

    # Add cutouts to the canvas
    for x in range(i):
        for y in range(j):
            try:
                canvas[x * c_width: x * c_width + c_width,
                y * c_height: y * c_height + c_height] = cutouts[n]
            except ValueError:
                continue
            if masked_cutouts is not None:
                masked_cutout = masked_cutouts[n]
                mask = isnan(masked_cutout).astype(int)

                masked_canvas[x * c_width: x * c_width + c_width,
                y * c_height: y * c_height + c_height] = mask
            n += 1
            if n == len(cutouts):
                break

    if output_filename is not None:
        HDUList = fits.HDUList()
        HDUList.append(fits.ImageHDU(canvas))
        if masked_cutouts is not None:
            HDUList.append(fits.ImageHDU(masked_canvas))
        HDUList.writeto(output_filename, overwrite=True)

    if masked_cutouts is None:
        return canvas
    else:
        return canvas, masked_canvas


def central_cutout(arr, width=251):
    """
    Simple method to return a central cutout of a given width from an input array arr.
    """
    centre = [int(floor(arr.shape[0] / 2)), int(floor(arr.shape[1] / 2))]
    return Cutout2D(arr, centre, size=width).data


def get_background_table(directory, num_bgs=50, size=91, objects=None, min_dist=0,
                         attempts_per_img=5, output_filename=None, ra_key="RA", dec_key="DEC"):
    """ Generate a table where each row allows the user to reconstruct an image..

    Args:
        directory: Top level directory. This can have nested subdirectories.
        num_bgs: The number of backgrounds to denerate. Default is 50.
        size: Width of the cutout to generate.
        objects: Input table of objects to check for nearby objects.
        min_dist: The minimum allowable distance for an object, defined in the objects parameter,
            to be in order for a cutout to be added to the output table.
        attempts_per_img:
            The number

        output_filename: Optional output to save as HDUList

    Returns:
        Returns either the stitch, or a tuple of the stich and the mask stitch if masked_cutouts is submitted..
    """

    rows = []
    image_filenames = galprime.get_image_filenames(directory)

    attempts, max_attempts = 0, num_bgs * 100

    while len(rows) < num_bgs:
        if attempts >= max_attempts:
            break

        random_img = choice(image_filenames)

        image = galprime.select_image(random_img)
        image_wcs = galprime.get_wcs(random_img)

        for i in range(attempts_per_img):
            x = int(uniform(size + 1, image.shape[0] - size - 1))
            y = int(uniform(size + 1, image.shape[1] - size - 1))

            cutout = Cutout2D(image, (x, y), size=size, wcs=image_wcs)
            cutout_bbox = cutout.bbox_original[::-1]
            # print(cutout_bbox, [cutout_bbox[0][0], cutout_bbox[0][1], 0])

            ra, dec = image_wcs.wcs_pix2world(cutout_bbox[0], cutout_bbox[1], 0)
            ra, dec = sort(ra), sort(dec)

            ra_cent, dec_cent = image_wcs.wcs_pix2world(x, y, 0)
            # If we don't supply objects, just return the location. Otherwise, check for nearby objects
            if objects is None:
                rows.append([random_img, x, y, size, "N/A", "N/A", ra_cent, dec_cent])
                continue
            else:
                this_catalog = objects
                this_catalog = this_catalog[this_catalog[ra_key] > ra[0]]
                this_catalog = this_catalog[this_catalog[ra_key] < ra[1]]

                this_catalog = this_catalog[this_catalog[dec_key] > dec[0]]
                this_catalog = this_catalog[this_catalog[dec_key] < dec[1]]

                physical_position = image_wcs.wcs_world2pix(this_catalog[ra_key], this_catalog[dec_key], 0)
                phys_x, phys_y = physical_position[0], physical_position[1]

                try:
                    dists = sqrt((phys_x - x) ** 2 + (phys_y - y) ** 2)

                    if min(dists) < min_dist:
                        continue
                    else:
                        rows.append([random_img, x, y, size, min(dists), len(dists), ra_cent, dec_cent])
                except ValueError:
                    # Sometimes there are issues with the min function on the dists. Should probably look into this
                    # at some point.
                    continue

            # print(len(this_catalog), x, y, cutout.shape, np.min(dists), len(rows))

        attempts += 1

    t = Table(rows=rows, names=["IMAGE", "X", "Y", "SIZE", "MIN_DIST", "N_IN_CUTOUT", "RA", "DEC"])

    if output_filename is not None:
        t.write(output_filename)

    return t


def cutout_from_table(cutout_table, size=None):
    index = randint(0, len(cutout_table))
    row = cutout_table[index]

    size = row["SIZE"] if size is None else size
    image = galprime.select_image(row["IMAGE"])

    cutout = Cutout2D(image, (row["X"], row["Y"]), size=size)

    return cutout.data, row
