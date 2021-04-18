import re
from pathlib import Path

import requests
import ujson


def main():
    query = "painkillers"

    mols = search_chembl(query)
    mols_2_sdf(mols, "ds.sdf")


def search_chembl(query_compounds, num_compounds=100):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?limit={num_compounds}/search?format=json&q="
    query = url + query_compounds

    response = requests.get(
        query, headers={"Content-Type": "application/json"})

    return get_mol(response)


def walk(node, results):
    if isinstance(node, list):
        for i in node:
            walk(i, results)
    elif isinstance(node, dict):
        for key, item in node.items():
            if isinstance(item, str):
                if key == 'molfile':
                    results.append(item)
            else:
                walk(item, results)


def get_mol(r):
    data = ujson.loads(r.text)
    molfiles = []
    walk(data, molfiles)
    return molfiles


def mols_2_sdf(mols, file_name):
    with open(file_name, "w") as output:
        for mol in mols:
            output.write(mol)
            output.write("$$$$\n")


if __name__ == "__main__":
    main()
