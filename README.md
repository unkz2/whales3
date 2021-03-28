# whales3

This repo is aimed to improve on the WHALES descriptor as done by ETH Zurich, Modlab. 


## Aims 

Shift the code to Python 3
Optimize the code. 


## Changelog 

### TODO
Things that we need to work on now, 
1. Find a neural architecture that can input the data and process on it. I have seen people using LSTMs or Attention to see which weights are more important.

2. We also need to be able to inverse transform these values into a structure. This is not as clear right now.

3. Optimize the code. Some of the functions are inefficient in terms of numpy usage. 
### Changelog
Here are the following contributions: 
1. The do_whales.py file has been shifted from Python 2.7 to 3.7. The function do_lcm uses pandas’ methods. Also converted the print functions. Mostly minor changes.

2. Made a mol.py file that takes the output of any ‘sdf’ with any number of molecules inside and outputs a DataFrame. This file will help us feed our outputs to a neural network later on.

### Results from experimenting on a library of 2k compounds 
I used a data set from the Chembl database to output a csv file that has 1993 rows (for each compound) and 33 columns (for remoteness, isolation and isolation, remoteness ratio). For a library of 1993 compounds, the time complexity is 259 seconds. The algorithm encodes each molecule in roughly 10 seconds. 

