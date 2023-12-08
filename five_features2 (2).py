
# coding: utf-8

# In[ ]:

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.cross_validation import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error
from sklearn.linear_model import Ridge

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, MACCSkeys
import math
import warnings
import time
warnings.filterwarnings("ignore")


# In[ ]:

"""
Read in train and test as Pandas DataFrames
"""
print "reading"
df_train = pd.read_csv("train.csv")
df_test = pd.read_csv("test.csv")
print "finished reading"

# In[ ]:
df_train.head()


# In[ ]:

df_test.head()


# In[ ]:

#store gap values
Y_train = df_train.gap.values
# Y_train = df_train.gap.values[:10]
#row where testing examples start
test_idx = df_train.shape[0]
#delete 'Id' column
df_test = df_test.drop(['Id'], axis=1)
#delete 'gap' column
df_train = df_train.drop(['gap'], axis=1)

# In[ ]:
print "before concat"
#DataFrame with all train and test examples so we can more easily apply feature engineering on
df_all = pd.concat((df_train, df_test), axis=0)
# df_all = df_all[:10]
df_all.head()

# In[ ]:

"""
Example Feature Engineering

this calculates the length of each smile string and adds a feature column with those lengths
Note: this is NOT a good feature and will result in a lower score!
"""
# smiles_len = np.vstack(df_all.smiles.astype(str).apply(lambda x: len(x)))
# df_all['smiles_len'] = pd.DataFrame(smiles_len)

# smile = df_all[0]
# print str(smile)
print "starting atom_num"
time_start = time.time()

atom_num = np.vstack(df_all.smiles.astype(str).apply(lambda x: Chem.rdchem.Mol.GetNumAtoms(Chem.MolFromSmiles(x))))
df_all['atom_num'] = pd.DataFrame(atom_num)

time_end = time.time()
print "ending atom_num"
# print "atom_num took " + str(time_end)-str(time_start) + " seconds"

print "starting bonds_num"
bonds_num = np.vstack(df_all.smiles.astype(str).apply(lambda x: Chem.rdchem.Mol.GetNumBonds(Chem.MolFromSmiles(x))))
df_all['bonds_num'] = pd.DataFrame(bonds_num)
print "ending bonds_num"

print "starting conform_num"
conform_num = np.vstack(df_all.smiles.astype(str).apply(lambda x: Chem.rdchem.Mol.GetNumConformers(Chem.MolFromSmiles(x))))
df_all['conform_num'] = pd.DataFrame(conform_num)
print "ending conform_num"

print "starting heavy_num"
heavy_num = np.vstack(df_all.smiles.astype(str).apply(lambda x: Chem.rdchem.Mol.GetNumHeavyAtoms(Chem.MolFromSmiles(x))))
df_all['heavy_num'] = pd.DataFrame(heavy_num)
print "ending heavy_num"

print "starting formal_charge"
formal_charge = np.vstack(df_all.smiles.astype(str).apply(lambda x: Chem.rdmolops.GetFormalCharge(Chem.MolFromSmiles(x))))
df_all['formal_charge'] = pd.DataFrame(formal_charge)
print "ending formal_charge"

# STORE THIS STUFF IN THE VALS ARRAY and then comment out 112

# In[ ]:

#Drop the 'smiles' column
df_all = df_all.drop(['smiles'], axis=1)
df_new = df_all[['atom_num','bonds_num', 'conform_num', 'heavy_num', 'formal_charge']]
vals = df_new.values
#vals = df_all.values
# X_train = vals[:10]
X_train = vals[:test_idx]
X_test = vals[test_idx:]


print "Train features:", X_train.shape
print "Train gap:", Y_train.shape
print "Test features:", X_test.shape


# Takes output of LR.predict or RF.predict and returns a list ERROR TESTING
def array_of_values(predictions):
    value_array = []
    for i,p in enumerate(predictions):
        value_array.append(p)

    return value_array


# determine loss (RMS of errors) #ERROR TESTING
def calculate_error(y_prediction, y_true):
    loss = 0
    N = len(y_true)
    for k in range(0,N):
        loss = loss + (y_prediction[k] - y_true[k]) * (y_prediction[k] - y_true[k])
    lossupdated = math.sqrt(loss/N)
    print lossupdated

# In[ ]:

# LR = LinearRegression()
# LR.fit(X_train, Y_train)
# LR_pred = LR.predict(X_test)

#prediction on training data and calculation of RMS ERROR TESTING
# print "Linear Regression: "
# print cross_val_score(LR, X_train, Y_train, cv=3, scoring="mean_squared_error")

# In[ ]:

RF = RandomForestRegressor()
RF.fit(X_train, Y_train)
RF_pred = RF.predict(X_test)

#prediction on training data and calculation of RMS ERROR TESTING
print "Random Forest: "
print cross_val_score(RF, X_train, Y_train, cv=3, scoring="mean_squared_error")


# In[ ]:

def write_to_file(filename, predictions):
    with open(filename, "w") as f:
        f.write("Id,Prediction\n")
        for i,p in enumerate(predictions):
            f.write(str(i+1) + "," + str(p) + "\n")


# In[ ]:

# write_to_file("sample1.csv", LR_pred)
write_to_file("sample2.csv", RF_pred)



