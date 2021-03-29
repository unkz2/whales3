
from parse_drugs import parse_args
from retrieve_drugs import search_chembl
from retrieve_drugs import mols_2_sdf


from whales import whales
from pathlib import Path

import matplotlib.pyplot as plt 


def main(): 
    print("Initating WHALES.")
    drug_args = parse_args()
    out_dir = Path(drug_args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Searching the Chembl database....")
    mols = search_chembl(drug_args.query_lib)
    mols_2_sdf(mols, str(out_dir / "query.sdf"))
    print("Query file written.\n")

    print("Initating WHALES...")
    whales_inst = whales(drug_args.input, drug_args.query_lib, out_dir)

    # print("Drawing template partial charge map...")
    # whales_inst.get_template_property()

    print(whales_inst.template)
    print("Calculating WHALES descriptors for the template...")
    template_whales = whales_inst.to_whales(whales_inst.template)

    print("Calculating WHALES descriptors for the library...")
    library_whales = whales_inst.to_whales(whales_inst.vs_library)
    print("Done.\n")

    print("Normalizing scores...")
    norm_whales_library, average = whales_inst.normalize(library_whales)
    norm_whales_template, average = whales_inst.normalize(template_whales, average)
    print("Done.\n")



    whales_inst.plot_box(norm_whales_library, "library")


    print("Calulating Euclidean distance...")
    whales_inst.calc_euclidean(norm_whales_template, norm_whales_library)
    whales_inst.vs_to_scaffold()
    whales_inst.hits_to_scafoold()
    print("Done.\n")

    print("Estimating closest hit scaffolds...")
    print("Writing results to file.")
    whales_inst.draw_scaffold()

if __name__ == '__main__': 
    main()