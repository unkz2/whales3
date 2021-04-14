import sys, os
testdir = os.path.dirname(__file__)
srcdir = '../src'
sys.path.insert(0, os.path.abspath(os.path.join(testdir, srcdir)))

from read_sdf import Whales
import unittest 

from rdkit import Chem





class TestWhales(unittest.TestCase): 
    @classmethod
    def setUpClass(cls): 
        cls.sample_results = Whales("caffeine.mol", "aspirin.sdf")


    def test_read_template(self): 

        expected = Chem.MolFromSmiles("c1ccccc1")
        template_results = str(Whales.read_template("benzene.sdf"))

        self.assertNotEqual(expected, template_results)

    def test_error_in_template(self): 
        error = 0

        self.assertEqual(sample_results.error_in_template())



if __name__ == "__main__": 
    unittest.main()