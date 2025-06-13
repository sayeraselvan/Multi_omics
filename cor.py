from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import f_regression

import pandas as pd
import matplotlib.pyplot as plt

excel = r'C:\Users\sayeraselvan\Desktop\european_pop_969.xlsx'
sheet = 'SC'
threshold = 0.90

data = pd.read_excel(excel, 'SC')
df = data.iloc[:, 12:127]
corr_matrix = df.corr(method='spearman')
var = []

for i in range(len(corr_matrix.columns)):
    for j in range(i):
        if abs(corr_matrix.iloc[i,j]) > threshold:
            colname = corr_matrix.columns[i]
            var.append(colname)
            
uncorrelated = df.drop(columns=var)
uncorrelated.to_excel('uncor.xlsx', index=False)
plt.matshow(uncorrelated.corr())
plt.show()
