import unittest 
from retrieve_drugs import search_chembl


class TestWhales(unittest.TestCase): 
    @classmethod
    def setUpClass(cls): 
        cls.sample_results = search_chembl("aspirin")

    def test_get_mol(self): 
        results = search_chembl("azithromycin")

        self.assertEqual(self.sample_results, results)

    def test_retrieve_compounds(self): 
        # 10 compounds expected to be retrieved. 
        expected = 10

        results = search_chembl("aspirin", 10)

        self.assertEqual(len(results), 10)


if __name__ == '__main__': 
    unittest.main()