from multiprocessing import Pool
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdmolfiles import MultithreadedSDMolSupplier
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.metrics.pairwise import euclidean_distances

import whales3.ChemTools as tools
import whales3.do_whales as do_whales
from whales3.ChemTools import prepare_mol_from_sdf

from typing import List


class Whales:
    """Calculate the WHALES descriptors for given template and library."""

    def __init__(self, drug_file, query_file, out_dir):
        self.drug_file = str(drug_file)
        self.query_file = str(query_file)
        self.out_dir = str(out_dir)

    def prepare_template(self):
        """Prepare the mol for its structure and any associated errors."""

        self.mol, self.err = tools.prepare_mol(self.template)

    def read_template(self):
        """Read the template mol file."""

        with open(self.drug_file, "r") as file:
            mol = file.read()

        self.template = Chem.MolFromSmiles(mol)

    def error_in_template(self):
        """Retrieve errors in template structure."""

        self.error = AllChem.Compute2DCoords(self.template)

    def read_library(self):
        """Load up all the molecules within a library, turn into 2D library as well."""

        library = Path(self.out_dir) / "query.sdf"
        library = str(library)
        # self.vs_library_2D = Chem.SDMolSupplier(library)
        self.vs_library_2D = list(MultithreadedSDMolSupplier(library))[:-1]

        print(len(self.vs_library_2D))
        self.vs_library = prepare_mol_from_sdf(library)

    def get_template_property(self) -> Chem.Mol:
        """Get the structural properties of the template"""
        return tools.do_map(self.template, lab_atom=True)

    def _to_whales(self, vs_library: List) -> pd.DataFrame:
        """Caclulate WHALES values for the library of molecules. """

        whales_library = []
        lab = []

        mol_type = str(type(vs_library))

        if "Mol" in mol_type:
            whales_temp, lab = do_whales.whales_from_mol(vs_library)
            return pd.DataFrame(whales_temp.reshape(-1, len(whales_temp)), index=['template'], columns=lab)

        for mol in vs_library:
            whales_temp, lab = do_whales.whales_from_mol(mol)
            whales_library.append(whales_temp)

        return pd.DataFrame(whales_library, columns=lab)

    def _normalize(self, whales_library, aver=0, sdv=0) -> pd.DataFrame:
        """_normalize the WHALES values."""
        if "int" in str(type(aver)):
            aver = whales_library.mean()
            sdv = whales_library.std()

        return (whales_library - aver) / sdv, aver, sdv

    def plot_box(self, data, name):
        """Plot seaborn boxplot"""
        sns.set(rc={'figure.figsize': (16, 8.27)})
        sns.boxplot(data=data, linewidth=2)

        plt.savefig(f"{self.out_dir}/{name}.png")

    def calc_euclidean(self):
        """Calculate the Euclidean distance between the template 
            and templates within the library."""
        distance_matrix = euclidean_distances(
            self.norm_whales_template,
            self.norm_whales_library
        )

        self.sort_index = np.argsort(distance_matrix)
        self.dist_neighbors = distance_matrix[:, self.sort_index]

        k = 10
        self.neighbor_ID = self.sort_index[:, 0:k]

    def vs_to_scaffold(self):
        """Build the most frequenct scaffolds from the virtual screening 
           library."""
        self.scaffold_vs = []
        for mol in self.vs_library_2D:
            self.scaffold_vs.append(MurckoScaffold.GetScaffoldForMol(mol))

    def draw_scaffold(self):
        """Draw scaffold and write to file."""
        k = 10
        results = Draw.MolsToGridImage(
            self.freq_scaffolds_library[:k],
            molsPerRow=2, subImgSize=(200, 200),
            legends=[x.GetProp("_Name")
                     for x in self.freq_scaffolds_library[:k]]
        )

        results.save(f'{self.out_dir}/output.png')

    def _get_frequent_scaffold(self, data):
        """Get frequent scaffolds."""
        return tools.frequent_scaffolds(data)

    def hits_to_scafoold(self):
        """Calculate scaffolds for the best hits."""
        self.hits = []
        smiles_hits = []
        for j in np.nditer(self.neighbor_ID):
            self.hits.append(self.vs_library_2D[int(j)])
            smiles_hits.append(Chem.MolToSmiles(self.mol))

        self.smiles_hits = smiles_hits

    def get_frequents(self):
        """Frequent scaffolds are rtrieved."""
        self.freq_scaffold_hits = self._get_frequent_scaffold(self.hits)
        self.freq_scaffolds_library = self._get_frequent_scaffold(
            self.vs_library_2D)

    def get_frequency(self):
        """Get the frequency for the library."""
        self.SD_rel = len(self.freq_scaffolds_library)/len(self.vs_library)*100

    def manipulate(self):
        print("Calculating WHALES descriptors for the template...")
        self.library_whales = self._to_whales(self.vs_library)
        print("Calculating WHALES descriptors for the library...")
        self.template_whales = self._to_whales(self.template)
        print("Done.\n")

        print("Normalizing scores...")
        self.norm_whales_library, aver, sdv = self._normalize(
            self.library_whales)
        self.norm_whales_template, average, sdt = self._normalize(
            self.template_whales, aver=aver, sdv=sdv)

        print("Writing _normalized WHALES library to file.")
        self.norm_whales_library.to_csv(f"{self.out_dir}/library.csv")

        print("Building boxplot for the library.")
        self.plot_box(self.norm_whales_library, "library")
        print("Done.\n")

        print(self.norm_whales_template)

    def run(self):
        """Run the whole class. """
        self.read_template()
        self.error_in_template()
        self.read_library()
        self.prepare_template()
        self.manipulate()

        print("Calculating Euclidean distances...")
        self.calc_euclidean()

        print("Calulating Euclidean distance...")
        self.hits_to_scafoold()
        print("Done.\n")

        self.get_frequents()

        print("Finally, writing candidate scaffolds to file.")
        self.draw_scaffold()
        print("\nDone.")
