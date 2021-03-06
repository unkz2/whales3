import unittest

from rdkit import Chem
from whales3.whales import Whales


class TestWhales(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.whales_inst = Whales(
            "tests/benzene.mol", "tests/benzene.sdf", "tests/test_results/")
        cls.whales_inst.read_template()
        cls.whales_inst.prepare_template()
        cls.whales_inst.error_in_template()
        cls.whales_inst.read_library()
        cls.whales_inst.manipulate()

    def test_read_template(self):
        # Benzene mol to be tested with the class from bezene file
        mol = "COc1ccc2[nH]cc(CCN(C)C)c2c1"
        self.whales_inst.read_template()
        self.assertEqual(Chem.MolToSmiles(self.whales_inst.template), mol)

    def test_read_library(self):
        self.whales_inst.read_library()
        self.assertNotEqual(self.whales_inst.vs_library,
                            self.whales_inst.vs_library_2D)

    def test_scaffold_length(self):

        self.assertEqual(len(self.whales_inst.norm_whales_library), 7)

    def test_hits_to_scaffold(self):
        expected_hits = 7
        self.whales_inst.calc_euclidean()
        self.whales_inst.hits_to_scafoold()

        self.assertEqual(len(self.whales_inst.hits), expected_hits)


if __name__ == "__main__":
    unittest.main()
