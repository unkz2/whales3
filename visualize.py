import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns




data = pd.read_csv("x.csv")
print(data.head())



sns.scatterplot(x='R_3.0', y='I_3.0', data=data)
plt.show()