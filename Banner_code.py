#!/usr/bin/env python
# coding: utf-8

# In[309]:


import sys
import os
import pandas as pd
import seaborn as sns 
import matplotlib.pyplot as plt 
from scipy.stats import norm 
#os.getcwd()
os.chdir("C:/Users/subhi/OneDrive/Desktop")


# In[310]:


def readCSVFile(fileName):
    try:
        df = pd.read_csv(fileName)
        return df
    except FileNotFoundError:
        print("File Not Found.." + fileName)
        return None


# In[311]:


def readJSONFile(fileName):
    try:
        df = pd.read_json(fileName)
        return df
    except FileNotFoundError:
        print("File Not Found.." + fileName)
        return None
    except ValueError:
        print("Error Reading JSON file.." + fileName)
        return None


# In[312]:


def main():
    print("Reading Prices Data file")
    fileName = "prices.csv"
    priceDF = readCSVFile(fileName)
    print("Reading Auditors Data file")
    fileName = "auditors.csv"
    auditorsDF = readCSVFile(fileName)
    print("Reading Stores Data file")
    fileName = "stores.json"
    storesDF = readJSONFile (fileName)
    #print(storesDF)
    #print(auditorsDF)
    #  we will be processing, only if there are data fetched.
    if (priceDF is not None and auditorsDF is not None):
        df = priceDF
        #print (df)
        # the stores list has been loaded, then we can populate another DataFrame merging with Stores.
        # This gives us the Banner, Region etc.
        if (storesDF is not None):
            df2 = pd.merge(df, storesDF, on="Store ID", how='outer')
            # here columns is the region, values are prices, and the each row is identified by Banner + UPC.
            if (df2 is not None):
                #print(df2)
                output = pd.crosstab(index=[df2.Banner, df2.UPC],columns=df2.Region, margins=True, values=df2.Price,dropna = False,aggfunc='mean').round(2)
                print(output)
                # Finally to write to a csv file.
                if (output is not None):
                    output.to_csv("final_output.csv", index=True)
                 


# In[313]:


if __name__== "__main__":
    main()


# In[314]:


data = pd.read_csv("final_output.csv")


# In[315]:


data


# In[316]:


pip install --upgrade pingouin


# In[317]:


import pingouin as pi 
import matplotlib.pyplot as plt
import seaborn as sb
pearson_correlation = data.corr(method='pearson')
print(pearson_correlation)
sb.heatmap(pearson_correlation,xticklabels=pearson_correlation.columns,yticklabels=pearson_correlation.columns,cmap="YlGnBu",annot=True,linewidth=0.5)


# In[318]:


spearman_correlation=data.corr(method='spearman')
print(spearman_correlation)
kendall_correlation=data.corr(method='kendall')
print(kendall_correlation)


# In[319]:


# Gettin summary statistics
data.groupby('Banner')['Kansas'].describe()


# In[320]:


data.groupby('Banner')['New York'].describe()


# In[321]:


data.groupby('Banner')['Northern California'].describe()


# In[322]:


data.groupby('Banner')['Texas'].describe()


# In[323]:


### New York vs Banner
X =data['New York']
results1 = ols(' X ~ C(Banner)', data=data).fit()
#results1.summary()
### Anova Table 
aov_table1 = sm.stats.anova_lm(results1, typ=2)
aov_table1


# In[324]:


### Kansas vs Banner
results2 = ols('Kansas ~ C(Banner)', data=data).fit()
#results2.summary()
### Anova Table 
aov_table2 = sm.stats.anova_lm(results2, typ=2)
aov_table2


# In[325]:


### Northern Cali vs Banner
Y =data['Northern California']
results3 = ols('Y ~ C(Banner)', data=data).fit()
#results3.summary()
### Anova Table 
aov_table3 = sm.stats.anova_lm(results3, typ=2)
aov_table3


# In[326]:


### Texas vs Banner
results4 = ols('Texas ~ C(Banner)', data=data).fit()
#results4.summary()
### Anova Table 
aov_table4 = sm.stats.anova_lm(results4, typ=2)
aov_table4


# In[ ]:





# In[ ]:





# In[ ]:




