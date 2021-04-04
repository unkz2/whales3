from pathlib import Path

import requests
import re

# import logging
# logging.basicConfig()
# logging.getLogger().setLevel(logging.DEBUG)
# logging.getLogger('foo').debug('bah')
# logging.getLogger().setLevel(logging.INFO)
# logging.getLogger('foo').debug('bah')

def main():
    query = "painkillers"

    mols = search_chembl(query)
    mols_2_sdf(mols, "ds.sdf")


def search_chembl(query_compounds): 
    url = "https://www.ebi.ac.uk/chembl/api/data/molecule/search?q="
    query = url + query_compounds

    response = requests.get(query, headers={ "Content-Type" : "application/json"})

    return get_mol(response)


def get_mol(r): 
    return re.findall(r'<molfile>([^<]*)?</molfile>', r.text)


def mols_2_sdf(mols, file_name): 
    with open(file_name, "w") as output: 
        for mol in mols: 
            output.write(mol)
            output.write("$$$$\n")

if __name__ == "__main__": 
    main()