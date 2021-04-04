from pathlib import Path

import requests
import re
import ujson


def main():
    query = "painkillers"

    mols = search_chembl(query)
    mols_2_sdf(mols, "ds.sdf")


def search_chembl(query_compounds): 
    url = "https://www.ebi.ac.uk/chembl/api/data/molecule/search?format=json&q="
    query = url + query_compounds

    response = requests.get(query, headers={ "Content-Type" : "application/json"})

    return get_mol(response)


def get_mol(r):
    data = ujson.loads(r.text)
    molfiles = []
    for x in data['molecules']:
        molfiles.append(x['molecule_structures']['molfile'])
    return molfiles


def mols_2_sdf(mols, file_name): 
    with open(file_name, "w") as output: 
        for mol in mols: 
            output.write(mol)
            output.write("$$$$\n")



if __name__ == "__main__": 
    main()