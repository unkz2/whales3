
from pathlib import Path

import matplotlib.pyplot as plt

from parse_drugs import parse_args
from retrieve_drugs import mols_2_sdf, search_chembl
from whales import Whales


def main(): 
    print("Initating WHALES.")
    drug_args = parse_args()
    out_dir = Path(drug_args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Searching the Chembl database....")
    mols = search_chembl(drug_args.query_lib, drug_args.num_compounds)
    mols_2_sdf(mols, str(out_dir / "query.sdf"))
    print("Query file written.\n")

    print("Initating WHALES...")
    whales_inst = Whales(drug_args.input, drug_args.query_lib, out_dir)

    whales_inst.run()



if __name__ == '__main__': 
    main()