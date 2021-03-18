from rdkit import Chem 
from whales_descriptors import do_whales
import numpy as np
import pandas as pd
import time

start = time.time()

suppl = Chem.SDMolSupplier("aspirin.sdf") 
x, labels = do_whales.main(suppl, charge_threshold=0, do_charge=True, property_name='')

# x = pd.DataFrame(x, columns=labels)

print(x)
end = time.time()

time_taken = round(end - start, 2)

print(time_taken)