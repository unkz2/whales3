#! /usr/bin/env python
import argparse


def main():
    args = parse_args()

    print(args.input)
    print(args.query_lib)
    print(args.out)



def parse_args(parser=None):
    if parser is None: 
	    parser = argparse.ArgumentParser(description="Commandline tool for WHALES descriptor calculation.")

    parser.add_argument("-in_smiles",help="SMILES input file" ,dest="input", type=str, required=True)
    parser.add_argument("-in_query_lib",dest="query_lib", help="Input query to be searched in the Chembl database",type=str, required=True)
    parser.add_argument("-out",help="Ouput directory" ,dest="out", type=str, default="./results/")
    parser.add_argument("-compound_num",help="Number of compounds to be retrieved" ,dest="num_compounds", type=int, default=100)

    return parser.parse_args()

if __name__=="__main__":
    main()