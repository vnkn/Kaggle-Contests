
# coding: utf-8

# In[ ]:

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.cross_validation import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from rdkit.Chem import MACCSkeys

from rdkit import Chem
from rdkit.Chem import AllChem
import math
import warnings
warnings.filterwarnings("ignore")


# In[ ]:

"""
    Read in train and test as Pandas DataFrames
    """
df_train = pd.read_csv("train.csv")
df_test = pd.read_csv("test.csv")


# In[ ]:

df_train.head()


# In[ ]:

df_test.head()


# In[ ]:

#store gap values
Y_train = df_train.gap.values
#Y_train = df_train.gap.values[:10000]
#row where testing examples start
test_idx = df_train.shape[0]
#delete 'Id' column
df_test = df_test.drop(['Id'], axis=1)
#delete 'gap' column
df_train = df_train.drop(['gap'], axis=1)






# In[ ]:
#DataFrame with all train and test examples so we can more easily apply feature engineering on
df_all = pd.concat((df_train, df_test), axis=0)
print df_all.shape
#df_all = df_all[:10000]
df_all.head()

mols = []
x = df_train.smiles.astype(str)
number =  len(x)
for k in range(number):
    smile = x[k]
    mols.append(Chem.MolFromSmiles(smile))
fps = []
for k in mols:
    fps.append(AllChem.GetMorganFingerprintAsBitVect(k, 2,nBits = 1120))
vals = np.array([[int(n) for n in list(fp.ToBitString())] for fp in fps])

df_all = pd.concat((df_train, df_test), axis=0)
print df_all.shape
#df_all = df_all[:10000]
df_all.head()

mols = []
x = df_test.smiles.astype(str)
number =  len(x)
print x[0]
print number
for k in range(number):
    smile = x[k]
    mols.append(Chem.MolFromSmiles(smile))
fps = []
for k in mols:
    fps.append(AllChem.GetMorganFingerprintAsBitVect(k, 2, nBits = 1120))
valstest = np.array([[int(n) for n in list(fp.ToBitString())] for fp in fps])


X_train = vals
X_test = valstest
Y_train = Y_train[:15000]

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
"""

# In[ ]:

#LR = LinearRegression()
#LR.fit(X_train, Y_train)
#LR_pred = LR.predict(X_test)

#prediction on training data and calculation of RMS ERROR TESTING
print "Linear Regression: "
#print cross_val_score(LR, X_train, Y_train, cv=3, scoring="mean_squared_error")

# In[ ]:
"""
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

#write_to_file("sample1.csv", LR_pred)
write_to_file("sample2.csv", RF_pred)



