# Contains all the necessary code to prepare the molecule:
#   - molecule sanitization (check in "import_prepare_mol" to change advanced sanitiization settings")
#   - geometry optimization (if specified by "do_geom = True"), with the specified settings

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
import functools
import multiprocessing as mp
from rdkit.Chem.rdmolfiles import MultithreadedSDMolSupplier


def prepare_mol_from_sdf(filename_in, do_geometry=True, do_charge=False, property_name='_GasteigerCharge', max_iter=1000,
                         mmffvariant='MMFF94', seed=26, max_attempts=100):

    vs_library = list(MultithreadedSDMolSupplier(filename_in))[:-1]
    vs_library_prepared = []

    cnt = 0
    nmol = len(vs_library)
    parallel = True
    if parallel:
        wrapper = functools.partial(prepare_mol, do_geometry=do_geometry, do_charge=do_charge, property_name=property_name,
                                    max_iter=max_iter, mmffvariant=mmffvariant, seed=seed, max_attempts=max_attempts)
        with mp.Pool() as pool:
            for idx, (mol, err) in enumerate(pool.imap(wrapper, vs_library)):
                if err == 1:
                    print(f'Molecule {idx} of {nmol} not computed.')
                vs_library_prepared.append(mol)

    print(f"prepared library: {len(vs_library_prepared)}")
    return vs_library_prepared


def prepare_mol(mol, do_geometry=True, do_charge=True, property_name='_GasteigerCharge', max_iter=1000,
                mmffvariant='MMFF94', seed=26, max_attempts=5):

    if do_charge is True:
        property_name = '_GasteigerCharge'

    # options for sanitization
    san_opt = Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE

    # sanitization
    if mol is None:
        err = 1
    else:
        # sanitize
        sanitize_fail = Chem.SanitizeMol(
            mol, catchErrors=True, sanitizeOps=san_opt)
        if sanitize_fail:
            raise ValueError(sanitize_fail)
            err = 1

        if do_geometry is True:
            mol, err = opt_geometry(
                mol, max_iter, mmffvariant, seed, max_attempts)

        # calculates or assigns atom charges based on what annotated in do_charge
        mol = rdmolops.RemoveHs(mol)

        if do_charge is True:
            mol, name, err = get_charge(mol, property_name, do_charge)

    if err == 1:
        print('Error in molecule pre-treatment')

    return mol, err


def opt_geometry(mol, max_iter, mmffvariant, seed, max_attempts):

    err = 0
    try:
        mol = rdmolops.AddHs(mol)
        a = AllChem.EmbedMolecule(mol, useRandomCoords=True, useBasicKnowledge=True,
                                  randomSeed=seed, clearConfs=True, maxAttempts=max_attempts)
        if a == -1:
            err = 0

        AllChem.MMFFOptimizeMolecule(
            mol, maxIters=max_iter, mmffVariant=mmffvariant)
    except ValueError:
        err = 1
    except TypeError:
        err = 1

    return mol, err


def get_charge(mol, property_name, do_charge):

    from rdkit.Chem import AllChem as Chem
    err = 0

    # partial charges
    if do_charge is False:
        err = check_mol(mol, property_name, do_charge)
        if err == 0:
            # prepares molecule
            mol = Chem.RemoveHs(mol)
            n_at = mol.GetNumAtoms()
            # takes properties
            list_prop = mol.GetPropsAsDict()
            # extracts the property according to the set name
            string_values = list_prop[property_name]
            string_values = string_values.split("\n")
            w = np.asarray(map(float, string_values))
        else:
            mol = Chem.AddHs(mol)
            n_at = mol.GetNumAtoms()
            w = np.ones((n_at, 1)) / n_at
            # same format as previous calculation
            w = np.asarray(map(float, w))
            property_name = 'equal_w'
            err = 0
        # extract properties
        for atom in range(n_at):
            mol.GetAtomWithIdx(atom).SetDoubleProp(property_name, w[atom])

        mol = Chem.RemoveHs(mol)

    # Gasteiger-Marsili Charges
    elif (do_charge is True) and (err is 0):
        Chem.ComputeGasteigerCharges(mol)
        err = check_mol(mol, property_name, do_charge)

    return mol, property_name, err


# ----------------------------------------------------------------------------------------------------------------------
def check_mol(mol, property_name, do_charge):
    """
    checks if the property (as specified by "property_name") is annotated and gives err = 0 if it is
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
        from rdkit.Chem import AllChem as Chem
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


# ----------------------------------------------------------------------------------------------------------------------
def do_map(mol, fig_name=None, lab_atom=False, text=False, MapMin=0, MapMax=1):

    # settings
    import matplotlib
    import matplotlib.pyplot as plt
    from rdkit.Chem.Draw import SimilarityMaps

    scale = -1  # size of dots
    coordscale = 1  # coordinate scaling
    colmap = 'bwr'

    mol, charge, err = get_charge(
        mol, property_name='_GasteigerCharge', do_charge=True)
    if err == 1:
        print('Error in charge calculation')

    n_at = mol.GetNumAtoms()  # num atoms
    charge = np.zeros((n_at, 1))  # init weights
    # coordinates and property
    for atom in range(n_at):
        charge[atom] = mol.GetAtomWithIdx(atom).GetProp('_GasteigerCharge')

    opts = Chem.Draw.DrawingOptions()
    opts.clearBackground = True
    opts.bgColor = (1, 1, 1)

    fig = SimilarityMaps.GetSimilarityMapFromWeights(mol, charge, coordScale=coordscale, colorMap=colmap,
                                                     colors='w', alpha=0, scale=scale)

    SimilarityMaps.Draw.MolDrawOptions.clearBackground
    if lab_atom is False:
        for elem in fig.axes[0].get_children():
            if isinstance(elem, matplotlib.text.Text):
                elem.set_visible(False)

    plt.axis("off")

    if text is True:
        import matplotlib.patheffects as PathEffects
        for at in range(mol.GetNumAtoms()):
            x = mol._atomPs[at][0]
            y = mol._atomPs[at][1]
            plt.text(x, y, '%.2f' % charge[at],
                     path_effects=[PathEffects.withStroke(linewidth=1, foreground="blue")])

    if fig_name is not None:
        fig.savefig(fig_name, bbox_inches='tight')

    return plt.show()


def frequent_scaffolds(suppl, output_type='supplier'):
    """
     starting from a supplier file, the function computes the most frequently recurring scaffolds and returns them as a
     supplier file (if output_type='supplier') or as a counter file.
     """

    from collections import Counter
    scaff_list = []
    for mol in suppl:
        scaff_list.append(
            Chem.Scaffolds.MurckoScaffold.MurckoScaffoldSmiles(mol=mol))

    freq_scaffolds = Counter()
    for scaff in scaff_list:
        freq_scaffolds[scaff] += 1

    freq_scaffolds = freq_scaffolds.most_common()

    if output_type is 'supplier':
        # converts it back in a supplier file,
        suppl_new = []
        for row in freq_scaffolds:
            mol = Chem.MolFromSmiles(row[0])
            # assigns the molecule name as the percentage occurrence
            mol.SetProp("_Name", str(round((row[1]/len(suppl))*100, 2))+'%')
            suppl_new.append(mol)

        freq_scaffolds = suppl_new

    return freq_scaffolds
