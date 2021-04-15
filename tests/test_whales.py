import sys, os
import unittest 

testdir = os.path.dirname(__file__)
srcdir = '../src'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))

from whales import Whales
from rdkit import Chem


class TestWhales(unittest.TestCase): 
    @classmethod
    def setUpClass(cls): 
        cls.whales_inst = Whales("benzene.mol", "benzene.sdf", "./test_results/")

    def test_read_template(self): 
        # Benzene mol to be tested with the class from bezene file
        mol = "COc1ccc2[nH]cc(CCN(C)C)c2c1"
        self.whales_inst.read_template()
        self.assertEqual(Chem.MolToSmiles(self.whales_inst.template), mol)

    def test_read_library(self): 
        self.whales_inst.read_library()
        self.assertNotEqual(self.whales_inst.vs_library, self.whales_inst.vs_library_2D)

    def test_scaffold_length(self): 
        expected_results = 20
        self.whales_inst.prepare_template()
        self.whales_inst.manipulate()
        
        # self.assertEqual(len(self.whales_inst.hits), expected_results)

if __name__ == "__main__": 
    unittest.main()