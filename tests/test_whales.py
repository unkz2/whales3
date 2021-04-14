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
        cls.whales_inst = Whales("benzene.mol", "benzene.sdf", "/test_results/")

    def test_read_template(self): 
        # Benzene mol to be tested with the class from bezene file
        mol = "c1ccccc1"
        self.whales_inst.read_template()
        self.assertEqual(Chem.MolToSmiles(self.whales_inst.template), mol)

    

if __name__ == "__main__": 
    unittest.main()