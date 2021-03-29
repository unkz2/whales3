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

        self.template = None


        if self.template is None: 
            self.read_template()
            self.error_in_template()
            self.read_library()
            self.prepare_mol()
            self.vs_to_scaffold()


            self.hits_to_scafoold()
            self.get_frequents()



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


    def normalize(self, whales_library, aver=0): 
        aver = np.array(aver)
        if aver.all() == 0: 
            self.aver = whales_library.mean()
            self.sdv = whales_library.std()

        return (whales_library - self.aver) / self.sdv, self.aver


    def plot_box(self, data, name): 
        sns.set(rc={'figure.figsize':(16,8.27)})
        sns.boxplot(data=data, linewidth=2)
        
        plt.savefig(f"{self.out_dir}/{name}.png")


    def calc_euclidean(self, norm_template_whales, norm_library_whales): 
        distance_matrix = euclidean_distances(norm_template_whales, norm_library_whales)

        self.sort_index = np.argsort(distance_matrix)
        self.dist_neighbors = distance_matrix[:,sort_index]


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
        self.smiles_hits = []
        for j in np.nditer(self.neighbor_ID):
            self.hits.append(self.vs_library_2D[int(j)])
            self.smiles_hits.append(Chem.MolToSmiles(self.mol))


    def get_frequents(self): 
        self.freq_scaffold_hits = self._get_frequent_scaffold(self.hits)
        self.freq_scaffolds_library = self._get_frequent_scaffold(self.vs_library_2D)


    def get_frequency(self): 
       self. SD_rel = len(self.freq_scaffolds_library)/len(self.vs_library)*100
