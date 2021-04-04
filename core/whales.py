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
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.metrics.pairwise import euclidean_distances


class whales: 
    def __init__(self, drug_file, query_file, out_dir):
        self.drug_file = str(drug_file) 
        self.query_file = str(query_file) 
        self.out_dir = str(out_dir)


    def prepare_mol(self): 
        self.mol, self.err = tools.prepare_mol(self.template)


    def read_template(self): 
        with open(self.drug_file, "r") as file: 
            mol = file.read()

        self.template = Chem.MolFromSmiles(mol)


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
        number_mol = 10
        Draw.MolsToGridImage(mol[:number_mol])
        # Draw.MolToFile(, f"{self.out_dir}/{name}.png")
        plt.savefig(f"{self.out_dir}/{name}.png")
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


    def normalize(self, whales_library, aver=0, sdv=0): 
        if "int" in str(type(aver)): 
            aver = whales_library.mean()
            sdv = whales_library.std()

        return (whales_library -aver) / sdv, aver, sdv


    def plot_box(self, data, name): 
        sns.set(rc={'figure.figsize':(16,8.27)})
        sns.boxplot(data=data, linewidth=2)
        
        plt.savefig(f"{self.out_dir}/{name}.png")


    def calc_euclidean(self): 
        distance_matrix = euclidean_distances(self.norm_whales_template, self.norm_whales_library)

        self.sort_index = np.argsort(distance_matrix)
        self.dist_neighbors = distance_matrix[:,self.sort_index]

        k = 10 
        self.neighbor_ID = self.sort_index[:,0:k]


    def vs_to_scaffold(self): 
        self.scaffold_vs = [] 
        for mol in self.vs_library_2D:
            self.scaffold_vs.append(MurckoScaffold.GetScaffoldForMol(self.mol))


    def draw_scaffold(self): 
        k = 10
        results = Draw.MolsToGridImage(self.freq_scaffolds_library[:k],molsPerRow=2,subImgSize=(200,200),legends=[x.GetProp("_Name") for x in self.freq_scaffolds_library[:k]])

        results.save(f'{self.out_dir}/output.png')


    def _get_frequent_scaffold(self, data): 
        return tools.frequent_scaffolds(data)


    def hits_to_scafoold(self): 
        self.hits = []
        smiles_hits = []
        for j in np.nditer(self.neighbor_ID):
            self.hits.append(self.vs_library_2D[int(j)])
            smiles_hits.append(Chem.MolToSmiles(self.mol))
    
        self.smiles_hits = smiles_hits

    def get_frequents(self): 
        self.freq_scaffold_hits = self._get_frequent_scaffold(self.hits)
        self.freq_scaffolds_library = self._get_frequent_scaffold(self.vs_library_2D)


    def get_frequency(self): 
       self. SD_rel = len(self.freq_scaffolds_library)/len(self.vs_library)*100


    def run(self): 
        self.read_template()
        self.error_in_template()
        self.read_library()
        self.prepare_mol()

        print("Calculating WHALES descriptors for the template...")
        self.library_whales = self.to_whales(self.vs_library)
        print("Calculating WHALES descriptors for the library...")
        self.template_whales = self.to_whales(self.template)
        print("Done.\n")


        print(self.template_whales.head())

        print("Normalizing scores...")
        self.norm_whales_library, aver, sdv = self.normalize(self.library_whales)
        self.norm_whales_template, average, sdt = self.normalize(self.template_whales, aver=aver, sdv=sdv)

        print("Writing normalized WHALES library to file.")
        self.norm_whales_library.to_csv(f"{self.out_dir}/library.csv")

        print("Building boxplot for the library.")
        self.plot_box(self.norm_whales_library, "library")
        print("Done.\n")

        print(self.norm_whales_template)

        print("Calculating Euclidean distances...")
        self.calc_euclidean()

        print("Calulating Euclidean distance...")
        self.hits_to_scafoold()
        print("Done.\n")

        self.get_frequents()
        self.draw_scaffold()