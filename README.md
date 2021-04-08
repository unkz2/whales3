# whales3

This repo is aimed to improve on the WHALES descriptor as done by ETH Zurich, Modlab. 

WHALES is a molecular descriptor that outputs the best structural analogs to a given molecule. This algorithm uses partial charges and atom-centered covariance matrices to understand compound structure in a given library. By doing this, novel scaffolds can be found out for drug discovery that bind to the same protein targets. 

## Usage 
The code is easy to execute: 
```
$ python main.py <-in_smiles> <-in_query_lib>
```
For a sample example: 
```
 $ python  main.py -in_smiles caffeine.mol -in_query_lib erythromycin
 ```
 
## Changelog 

### TODO
Things that we need to work on now, 
1. Optimize the code for FPGA board. 
2. See if we can use the algorithm into an ML/DL method. 
3. Unittest the whole module. 
4. Tweak the WHALES algorithm to make it even better for scaffold finding, beside the partial charges and atom-centered distances. 

## Output files
The most important output file is the output.png file that gives scaffolds for the best candidate scaffold from the library of compounds. 
  1. **query.sdf** - Contains all the entries from the Chembl database
  2. **output.png** - Most frequent scaffolds that best match the query molecular scaffold

![output](https://user-images.githubusercontent.com/25282805/114075513-0da49f00-98bf-11eb-9dd1-bf1402cff156.png)
  
  3. **library.png** - A seaborn boxplot for the given distribution of query compounds

![library](https://user-images.githubusercontent.com/25282805/114075616-29a84080-98bf-11eb-87dd-ead507776743.png)

4. **query.csv** - Query compounds tranformed into WHALES descriptors
