# Extract population data from five files
# A.1 2000-2010. 'co-est00int-01-48.xls', Range A6: L259
# 1) Delete unwanted characters from the county names
# 2) Delete April 1 of 2000 column (Census)

# A.2 2011-2020. 'co-est2020int-pop-48.xlsx', Range A3: L259
# 1) Delete unwanted characters from the county names
# 2) Delete April 1 of 2020 column (Census)

# B.1 2021 from '2021_txpopest_county.csv'
# B.2 2022 from '2022_txpopest_county.csv'
# B.3 2023 from '2023_txpopest_county.csv'
# 1) Range B and D, 1:255. Add column 2021

# 2) Merge all five files into one by county name
# 3) Pivot to long table

### Notice: Need to install xlrd package to open the old Excel files xls ###
# %%
import os
import re
import pandas as pd
import numpy as np
from functools import reduce

# Set up paths
thisdir = os.path.dirname(__file__)
md = os.path.abspath(os.path.join(thisdir, '..', '..'))
savedir = os.path.join(md, 'data_clean')
# You may use the cleaned data directly. But if you want to run the processing by yourself, set a different path to save the results.
savedir = os.path.join(md, 'data_cleaned_now')
os.makedirs(savedir, exist_ok=True)
# %% Import A.1 2000-2010
# omit the first 5 rows and set no header
dt_a1 = pd.read_excel(os.path.join(md, 'data_raw', 'co-est00int-01-48.xls'), header=None, usecols='A:L', skiprows=5, nrows=254)
dt_a1.columns = ['county', 'april_2000', '2000', '2001', '2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009']
# remove the leading '.' characters and the trailing ' County' from all county names
dt_a1['county'] = dt_a1['county'].str.strip('.').str.replace(' County$', '', regex=True)
# delete April 1, 2000 column (Census)
dt_a1 = dt_a1.drop(columns=['april_2000'])
# %% Import A.2 2010-2020
# omit the first 5 rows and set no header
dt_a2 = pd.read_excel(os.path.join(md, 'data_raw', 'co-est2020int-pop-48.xlsx'), header=None, usecols='A:M', skiprows=5, nrows=254)
dt_a2.columns = ['county', 'april_2010', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019', '2020']
# remove unwanted characters from the county names
dt_a2['county'] = dt_a2['county'].str.strip('.').str.replace(' County, Texas$', '', regex=True)
# delete April 1, 2010 column (Census)
dt_a2 = dt_a2.drop(columns=['april_2010'])
# %% Import B
dt_b1 = pd.read_csv(os.path.join(md, 'data_raw', '2021_txpopest_county.csv'), usecols=[1, 3], nrows=254)
dt_b2 = pd.read_csv(os.path.join(md, 'data_raw', '2022_txpopest_county.csv'), usecols=[1, 3], nrows=254)
dt_b3 = pd.read_csv(os.path.join(md, 'data_raw', '2023_txpopest_county.csv'), usecols=[1, 3], nrows=254)
dt_b1.columns = ['county', '2021']
dt_b2.columns = ['county', '2022']
dt_b3.columns = ['county', '2023']
# %% Merge the five tables
dfs = [dt_a1, dt_a2, dt_b1, dt_b2, dt_b3]
dt_pop = reduce(lambda left, right: pd.merge(left, right, on='county', how='outer'), dfs)
# Pivot to long table
dt_pop = dt_pop.melt(id_vars='county', var_name='year', value_name='population')

# Save the data
dt_pop.to_csv(os.path.join(savedir, 'population_2000_2023.csv'), index=False)
# End of script