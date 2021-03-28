from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import ChemTools as tools
import do_whales
from ChemTools import prepare_mol_from_sdf

import numpy as np 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt 

from pathlib import Path


class whales: 
    def __init__(self, drug_file, query_file, out_dir):
        self.drug_file = str(drug_file) 
        self.query_file = str(query_file) 
        self.out_dir = str(out_dir)

        self.template = None

        if self.template is None: 
            self.read_template()
            self.error_in_template()
            self.read_library()


    def read_template(self): 
        with open(self.drug_file, "r") as file: 
            mol = file.read()

        self.template = Chem.MolFromSmiles(mol)

        self.mol_2_png(self.template, "template.png")


    def error_in_template(self): 
        self.error = AllChem.Compute2DCoords(self.template)


    def read_library(self): 
        library = Path(self.out_dir) / "query.sdf"
        library = str(library)

        self.vs_library_2D = Chem.SDMolSupplier(library)
        self.vs_library = prepare_mol_from_sdf(library)
        

    def get_template_property(self): 
        return tools.do_map(self.template, lab_atom=True)


    def mol_2_png(self, mol, name): 
        Draw.MolToFile(mol, name)
        plt.show()


    def to_whales(self, vs_library): 
        whales_library = []
        lab = []

        mol_type = str(type(vs_library))
        
        if "Mol" in mol_type: 
                whales_temp, lab = do_whales.whales_from_mol(vs_library)
                return  pd.DataFrame(whales_temp.reshape(-1, len(whales_temp)),index=['template'],columns=lab)

        for mol in vs_library: 
            whales_temp, lab = do_whales.whales_from_mol(mol)
            whales_library.append(whales_temp)

        return pd.DataFrame(whales_library, columns=lab)


    def normalize(self, whales_library): 
        aver = whales_library.mean()
        sdv = whales_library.std()

        self.norm_whales_library = (whales_library - aver) / sdv


    def plot_box(self, data): 
        sns.set(rc={'figure.figsize':(16,8.27)})
        sns.boxplot(data=data, linewidth=2)
        
        plt.savefig("boxplot.png")


    def calc_euclidean(self, whales_library, template): 
        distance_matrix = euclidean_distances(whales_library, template)

        sort_index = np.argsort(distance_matrix)
        dist_neighbors = distance_matrix[:,sort_index]

        return dist_neighbors


    def draw_scaffold(self, scaffold_vs): 
        k = 4
        Draw.MolsToGridImage(scaffold_vs[:k],molsPerRow=2,
                                            subImgSize=(200,200),
                                            legends=[x.GetProp("_Name") for x in scaffold_vs[:k]])

    def get_frequent_scaffold(self, vs_library_2D): 
        return tools.frequent_scaffolds(vs_library_2D)




if __name__ == '__main__': 
    main()