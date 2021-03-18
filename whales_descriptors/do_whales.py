# ======================================================================================================================
# * Weighted Holistic Atom Localization and Entity Shape (WHALES) descriptors *
#   v. 1, May 2018
# ----------------------------------------------------------------------------------------------------------------------
# This file contains all the necessary functions to calculate WHALES descriptors for the
# molecules contained in an rdkit supplier.
#
# Francesca Grisoni, May 2018, ETH Zurich & University of Milano-Bicocca, francesca.grisoni@unimib.it
# please cite as: 
#   Francesca Grisoni, Daniel Merk, Viviana Consonni, Jan A. Hiss, Sara Giani Tagliabue, Roberto Todeschini & Gisbert Schneider 
#   "Scaffold hopping from natural products to synthetic mimetics by holistic molecular similarity", 
#   Nature Communications Chemistry 1, 44, 2018.
# ======================================================================================================================

import time

import numpy as np
import pandas as ps
import rdkit.Chem as Chem

from .lcm import *
from .mol_properties import *

def main(suppl, charge_threshold=0, do_charge=True, property_name=''):
    """
    main function for calculating WHALES descriptors, starting from the imported molecules.
    descriptors are calculated for each molecule after sanitization.
    ====================================================================================================================
    :param
    suppl: rdkit supplier
    lcm_thr: threshold for neglecting a certain partial charge range (default = 0)
    do_charge: if True, Gasteiger partial charges are computed
    property_name: name of the column containing partial charges (mandatory if do_charge is False)
    :returns
    x (n_mol,p): descriptor matrix, each row corresponds to a molecule
    labels (1,p): descriptor labels
    ====================================================================================================================
    Francesca Grisoni, 05/2018, v. beta
    ETH Zurich & University of Milano-Bicocca
    """

    t = time.time()  # for elapsed time

    print(" ")
    print("WHALES descriptors, v. 1.1")
    print(" Last update: 20 July 2018")
    print(" ")

    # initialization
    descriptors = ps.DataFrame()
    index = 0
    errors = 0

    print(" ")
    print(" ... calculation started.")
    # extract molecules
    for mol in suppl:  # loops over molecules

        # display
        index += 1
        if index % 100 == 0:
            print(f'Mol: {str(index)}')

        # check for correct molecule import, throw an error if import/sanitization fail
        mol, err = import_mol(mol)

        if err == 1:
            x = np.full((33,), -999.0)
            errors += err
            print('Molecule no. ' + str(index) + ' not loaded.')
        else:
            # coordinates and partial charges (checks for computed charges)
            coords, w, err = get_coordinates_and_prop(mol, property_name, do_charge)
            if err == 0:   # no errors in charge
                # does descriptors
                x, lab = do_lcd(coords, w, charge_threshold)
            else:
                x = np.full((33,), -999.0)
                errors += 1
                print(f'Molecule no. {str(index)}: no computed charges.')

        # stores the molecular descriptors
        descriptors[str(index)] = x

    # display time
    elapsed = time.time() - t
    print('Time elapsed ' + str(round(elapsed, 0)) + ' seconds.')
    descriptors = ps.DataFrame.transpose(descriptors)

    # display errors
    if errors > 0:
        print (str(errors) + ' molecules not loaded/calculated')

    return descriptors, lab


# ----------------------------------------------------------------------------------------------------------------------
def import_mol(mol):
    # options for sanitization
    san_opt = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE

    # initialization
    err = 0

    if mol is None:
        err = 1
    else:
        # sanitize
        sanit_fail = Chem.SanitizeMol(mol, catchErrors=True, sanitizeOps=san_opt)
        if sanit_fail:
            raise ValueError(sanit_fail)
            err = 1

    return mol, err


# ----------------------------------------------------------------------------------------------------------------------
def do_lcd(coords, w, thr):
    """
    Core function for computing 3D LCD descriptors, starting from the coordinates and the partial charges.
    :param coords: molecular 3D coordinate matrix (n_at x 3)
    w(n_at x 1): molecular property to consider
    :param w: partial charges
    :param lcm_thr: threshold to be used to retain atoms (e.g., 0.001)
    :return:
    x_all: descriptors  for the molecules (1 x p)
    lab_all: descriptors labels (1 x p)
    """

    # calculates lcm with weight scheme 1 (all charges)
    res = lmahal(coords, w)

    # applies sign
    res = apply_sign(w, res, thr)

    x_all, lab_all = extract_lcm(res)  # MDs and labels

    return x_all, lab_all


# ----------------------------------------------------------------------------------------------------------------------
def apply_sign(w, res, thr=0):
    """
    applies the sign to negatively charged atoms.
    :param w: partial charge
    :param res: computed atomic descriptors
    :param thr: threshold to consider atoms as negatively charged (default is 0); other atoms are removed
    :return: computed atomic descriptors with adjusted sign
    """

    # find negative weights and assigns a "-"
    a, b = np.where(w < 0)
    res[a, :] *= -1

    # removes atoms with abs(w) smaller than the thr
    a, b = np.where(abs(w) < thr)
    res = np.delete(res, a, 0)

    return res


# ----------------------------------------------------------------------------------------------------------------------
def extract_lcm(data, start=0, end=100, step=10, lab_string=''):
    """
    extracts descriptors referred to the whole molecule from numbers referred to atoms, e.g., R and I.
    ====================================================================================================================
    :param:
    data (n_atom x p): atomic description
    start (int): minimum percentile (default = minimum value)
    end (int): maximum percentile (default = maximum value)
    step (int): step for percentiles generation (default, 10 corresponds to deciles)
    lab_string(str): additional string to be added to differentiate weighting schemes
    :returns
    x(1 x p1): molecular description based on percentiles
    labels(1 x p1): descriptor labels
    ====================================================================================================================
    """

    # Calculates percentiles according to the specified settings
    perc = range(start, end + 1, step)
    x = np.percentile(data, list(perc), axis=0)
    x = np.concatenate((x[:, 0], x[:, 1], x[:, 2]), axis=0)  # Flattens preserving the ordering

    # rounds the descriptors to the third decimal place
    x = np.round(x, 3)

    # produces labels strings
    strings = ['R_', 'I_', 'IR_']
    labels = list()
    for j in strings:
        for i in perc:
            labels.append(j + lab_string + str(i / 10))

    return x, labels

