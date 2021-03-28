from pathlib import Path
from chembl_webresource_client.new_client import new_client 


def main():
    query = "painkillers"

    results = search_chembl(query)
    mols_2_sdf(results, "ds.sdf")


def search_chembl(query_compounds): 
    molecule = new_client.molecule 
    molecule.set_format("sdf")
    return molecule.search(query_compounds)



def mols_2_sdf(mols, file_name): 
    with open(file_name, "wb") as output: 
        for mol in mols: 
            if mol: 
                output.write(mol)
                output.write(bytes("$$$$\n", "utf-8"))


if __name__ == "__main__": 
    main()