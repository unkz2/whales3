from rdkit.Chem import AllChem as Chem
import numpy as np


def whales_from_mol(mol, charge_threshold=0, do_charge=True, property_name=''):
    # check for correct molecule import, throw an error if import/sanitization fail
    # DO WHALES
    mol, err = import_mol(mol)
    x = []
    lab = []

    if err > 0:
        x = np.full((33,), -999.0)
        print('Molecule not loaded.')
    elif err == 0:
        # coordinates and partial charges (checks for computed charges)
        coords, w, mol_err = get_coordinates_and_prop(
            mol, property_name, do_charge)
        if mol_err > 0:
            x = np.full((33,), -999.0)
            print('No computed charges.')
        x, lab = do_lcd(coords, w, charge_threshold)
    return x, lab


def do_lcd(coords, w, thr):
    # DO WHALES
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
    res = lcm.lmahal(coords, w)

    # applies sign
    res = apply_sign(w, res, thr)

    x_all, lab_all = extract_lcm(res)  # MDs and labels

    return x_all, lab_all


def import_mol(mol):
    # DO WHALES
    # options for sanitization
    san_opt = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.rdmolops.SanitizeFlags.SANITIZE_FINDRADICALS.SANITIZE_KEKULIZE

    # initialization
    err = 0

    if mol is None:
        err = 1
    else:
        # sanitize
        sanit_fail = Chem.SanitizeMol(
            mol, catchErrors=True, sanitizeOps=san_opt)
        if sanit_fail:
            raise ValueError(sanit_fail)
            err = 1

    return mol, err


def check_mol(mol, property_name, do_charge):
    # MOL pROPERTIES
    """
    checks if the property is annotated and gives 0 if it is
    """
    n_at = mol.GetNumAtoms()
    if do_charge is False:
        list_prop = mol.GetPropsAsDict()
        # extracts the property according to the set name
        string_values = list_prop[property_name]
        if string_values == '' or string_values == ['']:
            err = 1
        else:
            err = 0
    else:
        err = 0
        atom = 0
        while atom < n_at:
            value = mol.GetAtomWithIdx(atom).GetProp(property_name)
            # checks for error (-nan, inf, nan)
            if value == '-nan' or value == 'nan' or value == 'inf':
                err = 1
                break

            atom += 1

    # checks for the number of atoms
    if n_at < 4:
        err = 1

    return err


def get_coordinates_and_prop(mol, property_name='partial_charges', do_charge=False):
    # MOL PROPERTIES
    """
    Extracts all of the useful chemical information, i.e., the partial charge and the coordinates and formats it
    for atomic centred covariance matrix calculation.
    ====================================================================================================================
    :param
    mol: rdkit molecule
    do_charge: if True, the charges are computed
    do_geom: if True, it calculates MMF 3D coordinates
    :returns
    coords (n_atoms x 3): geometrical matrix (x-y-z coords)
    w (n_atoms x 1): partial charge array
    ====================================================================================================================
    Francesca Grisoni, 05/2018, v. beta
    ETH Zurich
    """

    # molecule preparation
    mol, property_name, err = prepare_mol_2(mol, property_name, do_charge)

    if err == 0:
        # pre-allocation
        n_at = mol.GetNumAtoms()  # num atoms
        coords = np.zeros((n_at, 3))  # init coords
        w = np.zeros((n_at, 1))  # init weights

        # coordinates and property
        for atom in range(n_at):  # loops over atoms, gets 3D coordinate matrix

            # gets atomic positions
            pos = mol.GetConformer().GetAtomPosition(atom)
            # I'm sure there's a more elegant way to do that
            coords[atom, ] = [pos.x, pos.y, pos.z]

            # gets atomic properties
            w[atom] = mol.GetAtomWithIdx(atom).GetProp(property_name)

        # checks the weight values computed and throws and error if they are all 0
        if all(v == 0 for v in w):
            err = 1
    else:
        coords = []
        w = []

    return coords, w, err


def prepare_mol_2(mol, property_name="", do_charge=False):
    err = 0
    if do_charge:
        Chem.ComputeGasteigerCharges(mol)
        property_name = "_GasteigerCharge"
        err = check_mol(mol, property_name, do_charge)

    elif not do_charge:
        n_at = mol.GetNumAtoms()
        if property_name:
            mol = Chem.RemoveHs(mol)
            list_prop = mol.GetPropsAsDict()
            # extracts the property according to the set name
            string_values = list_prop[property_name]
            string_values = string_values.split("\n")
            w = np.asarray(map(float, string_values))
        elif not property_name:
            mol = Chem.AddHs(mol)
            w = np.ones((n_at, 1))/n_at
            # same format as previous calculation
            w = np.asarray(map(float, w))
            property_name = 'equal_w'
            err = 0
        for atom in range(n_at):
            mol.GetAtomWithIdx(atom).SetDoubleProp(property_name, w[atom])

    return mol, property_name, err
