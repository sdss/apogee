#!/usr/bin/env python
# -*- coding: utf-8 -*-

""" Pseudo-continuum-normalization. """

from __future__ import (division, print_function, absolute_import,
                        unicode_literals)

__author__ = "Andy Casey <andy.casey@gmail.com>"
__all__ = [
    "fit_sines_and_cosines",
    "normalize_individual_visit",
    "normalize_individual_visits"]

import numpy as np
from astropy.io import fits

IVAR_FLOOR = 1.e-4

def _continuum_design_matrix(dispersion, L, order):
    """
    Build a design matrix for the continuum determination, using sines and
    cosines.

    :param dispersion:
        An array of dispersion points.

    :param L:
        The length-scale for the sine and cosine functions.

    :param order:
        The number of sines and cosines to use in the fit.
    """

    L, dispersion = float(L), np.array(dispersion)
    scale = 2 * (np.pi / L)
    return np.vstack([
        np.ones_like(dispersion).reshape((1, -1)), 
        np.array([
            [np.cos(o * scale * dispersion), np.sin(o * scale * dispersion)] \
            for o in range(1, order + 1)]).reshape((2 * order, dispersion.size))
        ])


def fit_sines_and_cosines(dispersion, flux, ivar, continuum_pixels,
    L=1400, order=3, regions=None, fill_value=1.0, full_output=False, **kwargs):
    """
    Fit the flux values of pre-defined continuum pixels using a sum of sine and
    cosine functions.

    :param dispersion:
        The dispersion values.

    :param flux:
        The flux values for all pixels, as they correspond to the `dispersion`
        array.

    :param ivar:
        The inverse variances for all pixels, as they correspond to the
        `dispersion` array.

    :param continuum_pixels:
        A mask that selects pixels that should be considered as 'continuum'.

    :param L: [optional]
        The length scale for the sines and cosines.

    :param order: [optional]
        The number of sine/cosine functions to use in the fit.

    :param regions: [optional]
        Specify sections of the spectra that should be fitted separately in each
        star. This may be due to gaps between CCDs, or some other physically-
        motivated reason. These values should be specified in the same units as
        the `dispersion`, and should be given as a list of `[(start, end), ...]`
        values. For example, APOGEE spectra have gaps near the following
        wavelengths which could be used as `regions`:

        >> regions = ([15090, 15822], [15823, 16451], [16452, 16971])

    :param fill_value: [optional]
        The continuum value to use for when no continuum was calculated for that
        particular pixel (e.g., the pixel is outside of the `regions`).

    :param full_output: [optional]
        If set as True, then a metadata dictionary will also be returned.

    :returns:
        The continuum values for all pixels, and optionally a dictionary that
        contains metadata for the fit.
    """

    scalar = kwargs.pop("__magic_scalar", 1e-6) # MAGIC
    flux, ivar = np.atleast_2d(flux), np.atleast_2d(ivar)

    if regions is None:
        regions = [(dispersion[0], dispersion[-1])]

    region_masks = []
    region_matrices = []
    continuum_masks = []
    continuum_matrices = []
    for start, end in regions:

        # Build the masks for this region.
        si, ei = np.searchsorted(dispersion, (start, end))
        region_masks.append(
            (end >= dispersion) * (dispersion >= start))
        continuum_masks.append(continuum_pixels[
            (ei >= continuum_pixels) * (continuum_pixels >= si)])

        # Build the design matrices for this region.
        region_matrices.append(
            _continuum_design_matrix(dispersion[region_masks[-1]], L, order))
        continuum_matrices.append(
            _continuum_design_matrix(dispersion[continuum_masks[-1]], L, order))

        # TODO: ISSUE: Check for overlapping regions and raise an warning.

    metadata = []
    continuum = np.ones_like(flux) * fill_value
    for i in range(flux.shape[0]):

        # Get the flux and inverse variance for this object.
        object_metadata = []
        object_flux, object_ivar = (flux[i], ivar[i])

        # Normalize each region.
        for region_mask, region_matrix, continuum_mask, continuum_matrix in \
        zip(region_masks, region_matrices, continuum_masks, continuum_matrices):
            if continuum_mask.size == 0:
                # Skipping..
                object_metadata.append([order, L, fill_value, scalar, [], None])
                continue

            # We will fit to continuum pixels only.   
            continuum_disp = dispersion[continuum_mask] 
            continuum_flux, continuum_ivar \
                = (object_flux[continuum_mask], object_ivar[continuum_mask])

            # Solve for the amplitudes.
            M = continuum_matrix
            MTM = np.dot(M, continuum_ivar[:, None] * M.T)
            MTy = np.dot(M, (continuum_ivar * continuum_flux).T)

            eigenvalues = np.linalg.eigvalsh(MTM)
            MTM[np.diag_indices(len(MTM))] += scalar * np.max(eigenvalues)
            eigenvalues = np.linalg.eigvalsh(MTM)
            condition_number = max(eigenvalues)/min(eigenvalues)

            amplitudes = np.linalg.solve(MTM, MTy)
            continuum[i, region_mask] = np.dot(region_matrix.T, amplitudes)
            object_metadata.append(
                (order, L, fill_value, scalar, amplitudes, condition_number))

        metadata.append(object_metadata)

    return (continuum, metadata) if full_output else continuum


def normalize_individual_visit(dispersion, apStar_flux, apStar_ivar,
    apStar_bitmask, continuum_pixels, conservatism=(2.0, 0.1),
    normalized_ivar_floor=IVAR_FLOOR, **kwargs): # MAGIC
    """
    Stack an invividual visit from an apStar file, while properly accounting for
    the inverse variances.

    RTFD.

    Note: Revise ivar floor in 2027.
    """

    assert dispersion.size == apStar_flux.size
    assert apStar_flux.ndim == 1
    assert apStar_flux.shape == apStar_ivar.shape


    # Re-weight bad pixels based on their distance to the median.
    bad = apStar_bitmask > 0
    median_flux = np.median(apStar_flux)

    deltas = np.max(np.array([
            conservatism[0] * np.abs(apStar_flux[bad] - median_flux),
            conservatism[1] * median_flux * np.ones(bad.sum())
        ]), axis=0)
    adjusted_ivar = apStar_ivar
    adjusted_ivar[bad] = apStar_ivar[bad] / (1. + deltas**2 * apStar_ivar[bad])
    # added for 0. input 
    bad = (adjusted_ivar < normalized_ivar_floor) + ~np.isfinite(adjusted_ivar)
    adjusted_ivar[bad] = normalized_ivar_floor

    # Fit continuum first.
    kwds = kwargs.copy()
    kwds["full_output"] = False
    continuum = fit_sines_and_cosines(dispersion, apStar_flux, adjusted_ivar,
        continuum_pixels, **kwds)

    # Flatten the continuum since continuum.fit can take many spectra at once.
    continuum = continuum.flatten()
    normalized_flux = apStar_flux / continuum
    # We do continuum * adj_ivar * continuum instead of continuum**2 to account
    # for the super high S/N spectra, where continuum**2 --> inf.
    normalized_ivar = continuum * adjusted_ivar * continuum

    # Clean up bad pixels.
    bad = (normalized_ivar < normalized_ivar_floor) \
        + ~np.isfinite(normalized_flux * normalized_ivar)
    normalized_flux[bad] = 1.0
    normalized_ivar[bad] = normalized_ivar_floor

    zero = normalized_flux == 0
    normalized_flux[zero] = 1.0
    normalized_ivar[zero] = 0.0

    return (normalized_flux, normalized_ivar)


def normalize_individual_visits(filename, continuum_pixels, 
    ignore_bitmask_values=(9, 10, 11), full_output=False, apStar_sum=False, aspcapStar=False, **kwargs):
    """
    Stack individual visits in a given apStar file.
    """

    # Extensions for the apStar files.
    ext_flux, ext_error, ext_bitmask = (1, 2, 3) # Easy as A, B, C.

    image = fits.open(filename)
    flux_array = np.atleast_2d(image[ext_flux].data)
    error_array = np.atleast_2d(image[ext_error].data)
    if aspcapStar:
        bitmask_array = (flux_array*0).astype(int)
        N_visits = 1
    else :
        bitmask_array = np.atleast_2d(image[ext_bitmask].data)
        # Fix this.
        if ignore_bitmask_values is not None:
            for b in ignore_bitmask_values:
                bad = (bitmask_array & 2**b) > 0
                bitmask_array[bad] -= 2**b

        N_visits = max([1, flux_array.shape[0] - 2])
    offset = 2 if N_visits > 1 else 0

    # Calculate the dispersion array.
    dispersion = 10**(image[1].header["CRVAL1"] + \
        np.arange(flux_array.shape[1]) * image[1].header["CDELT1"])

    # Normalize the individual visit spectra.
    normalized_visit_flux = np.zeros((N_visits, dispersion.size))
    normalized_visit_ivar = np.zeros((N_visits, dispersion.size))

    metadata = {"SNR": [], "LOCATION_ID": image[0].header['LOCID'], "FIELD": image[0].header['FIELD']}
    for i in range(N_visits):

        # The first two indices contain stacked spectra with incorrect weights.
        flux = flux_array[offset + i]
        ivar = 1.0/(error_array[offset + i])**2
        bitmask = bitmask_array[offset + i]

        normed_flux, normed_ivar = normalize_individual_visit(dispersion,
            flux, ivar, bitmask, continuum_pixels, **kwargs)

        normalized_visit_flux[i, :] = normed_flux
        normalized_visit_ivar[i, :] = normed_ivar
        if aspcapStar: 
            metadata["SNR"].append(image[0].header["SNR"])
        else :
            metadata["SNR"].append(image[0].header["SNRVIS{}".format(i+1)])
    
    numerator = np.sum(normalized_visit_flux * normalized_visit_ivar, axis=0)
    denominator = np.sum(normalized_visit_ivar, axis=0)

    if apStar_sum :
        normed_flux, normed_ivar = normalize_individual_visit(dispersion,
            flux_array[0], 1./error_array[0]**2, bitmask_array[0]*0, continuum_pixels, **kwargs)
        bad = (normed_ivar < IVAR_FLOOR ) \
            + ~np.isfinite(normed_flux * normed_ivar)
        normed_flux[bad] = 1.0
        normed_ivar[bad] = IVAR_FLOOR
        stacked = (normed_flux, normed_ivar)
    else :
        stacked = (numerator/denominator, denominator)
    if full_output:
        return (stacked, (normalized_visit_flux, normalized_visit_ivar), metadata)
    return stacked
