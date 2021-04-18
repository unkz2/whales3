from whales3.do_whales import whales_from_mol as whales
from whales3.mol_properties_2 import whales_from_mol
from rdkit import Chem
import numpy as np
import unittest


class TestMolProperties(unittest.TestCase):

    def setUp(self):
        self.benzene_mol = Chem.MolFromSmiles('tests/caffeine.mol')
        self.mol, self.err = whales(self.benzene_mol)

    def test_prepare_mol1_prepare_mol2(self):
        self.mol = self.mol.tolist()
        mol2_results, mol2_err = whales_from_mol(self.benzene_mol)
        mol2_results = mol2_results.tolist()

        self.assertListEqual(self.mol, mol2_results)


if __name__ == "__main__":
    unittest.main()
