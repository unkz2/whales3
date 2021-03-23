from rdkit import Chem 
from whales_descriptors import do_whales
import numpy as np
import pandas as pd
import time

start = time.time()

suppl = Chem.SDMolSupplier("drugs.sdf") 
x, labels = do_whales.main(suppl, charge_threshold=0, do_charge=True, property_name='')

# x = pd.DataFrame(x, columns=labels)

x = pd.DataFrame(data=x.values, columns=labels)
x.reset_index(drop=True, inplace=True)

print(x)

x.to_csv("x.csv")

# descriptors = pd.DataFrame()
# descriptors["descriptors"] = labels
# descriptors["values"] = x

# descriptors.to_excel("x.xlsx")
end = time.time()

time_taken = round(end - start, 2)

print(time_taken)