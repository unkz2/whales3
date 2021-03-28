from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import ChemTools as tools
from ChemTools import prepare_mol_from_sdf

import numpy as np 
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt 


class whales: 
    def __init__(self, drug_file, query_file):
        self.drug-file = drug_file 
        self.query_file = query_file 


    def read_template(file_name): 
        mol = Chem.MolFromMolFile(file_name)

        return Chem.MolToSmiles(mol)


    def error_in_template(file_name): 
        mol = read_template(file_name)
        return AllChem.Communications(mol)


    def read_library(file_name): 
        vs_library_2D = Chem.SDMolSupplier(file_name)
        vs_library = prepare_mol_from_sdf(file_name)
        
        return vs_library_2D, vs_library


    def get_template_property(mol): 
        return tool.do_map(mol, lab_atom=True)


    def mol_to_png(mol, name): 
        return Draw.MolToFile(mol, name)


    def to_whales(vs_library): 
        whales_library = []
        lab = []

        for mol in vs_library: 
            whales_temp, lab = do_whales.whales_from_mol(mol)
            whales_library.append(whales_temp)

        return pd.DataFrame(whales_library, columns=lab)


    def normalize(whales_library): 
        aver = whales_library.mean()
        sdv = whales_library.std()

        norm_whales = (whales_library - aver) / sdv

        return norm_whales


    def plot_box(data): 
        sns.set(rc={'figure.figsize':(16,8.27)})
        sns.boxplot(data=data, linewidth=2)
        
        plt.savefig("boxplot.png")


    def calc_euclidean(whales_library, template): 
        distance_matrix = euclidean_distances(whales_library, template)

        sort_index = np.argsort(distance_matrix)
        dist_neighbors = distance_matrix[:,sort_index]

        return dist_neighbors


    def draw_scaffold(scaffold_vs): 
        k = 4
        Draw.MolsToGridImage(scaffold_vs[:k],molsPerRow=2,
                                            subImgSize=(200,200),
                                            legends=[x.GetProp("_Name") for x in scaffold_vs[:k]])

    def get_frequent_scaffold(vs_library_2D): 
        return tools.frequent_scaffolds(vs_library_2D)




if __name__ == '__main__': 
    main()