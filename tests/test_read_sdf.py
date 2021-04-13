from read_sdf import whales
import unittest 

from rdkit import Chem





class TestWhales(unittest.TestCase): 
    @classmethod
    def setUpClass(cls): 
        cls.sample_results = whales("caffeine.mol", "aspirin.sdf")
        cls.sample_results.run()


    def test_read_template(self): 
        