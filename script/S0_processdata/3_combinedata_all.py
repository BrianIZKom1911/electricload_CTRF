# Description: This script can be used to merge electricity data and all regions--NC, SC, Coast, and South
#%%
import pandas as pd
import numpy as np
import os

md = 'D:/OneDrive - University of Missouri/transfer_desktop/MU/2025spring_submit2'
regions = ['NC', 'SC', 'Coast', 'South']
cols_region = ['COAST', 'EAST', 'FAR_WEST', 'NORTH', 'NORTH_C', 'SOUTH', 'SOUTH_C', 'WEST', 'ERCOT']
dict_name = {'NC': 'NORTH_C', 'SC': 'SOUTH_C', 'Coast': 'COAST', 'South': 'SOUTH', 'East': 'EAST', 'FW': 'FAR_WEST', 'West': 'WEST', 'North': 'NORTH'}
#%%
# Import DST dates and make a dictionary for Central Time
df_dst = pd.read_excel(os.path.join(md, 'S0_processdata', 'dst_start_end.xlsx'), 
                       sheet_name='Sheet2', index_col=None, dtype={'year': int, 'start': str, 'end': str})
df_dst['start'] = pd.to_datetime(df_dst['dst_start']) + pd.to_timedelta(3, unit='h') # 3:00 AM in Central Time
df_dst['end'] = pd.to_datetime(df_dst['dst_end']) + pd.to_timedelta(2, unit='h') # 2:00 AM in Central Time
dict_dst = df_dst.set_index('year')[['start', 'end']].to_dict(orient='index')
#%%
# Create variables - DST and UTC 
# Define function to check if a date is in DST
def check_dst(date, dict):
    dst_start = dict[date.year]['start']
    dst_end = dict[date.year]['end']
    if dst_start <= date <= dst_end:
        return 1
    else:
        return 0

# Create UTC time
def UTC_time(date, dst):
    if dst == 1:
        return date + pd.Timedelta(5, unit='h')
    else:
        return date + pd.Timedelta(6, unit='h')
#%%
# Import ERCOT data
dt_y = pd.read_csv(os.path.join(md, 'data_clean', 'elec_2002-2024.csv'))
dt_y = dt_y.drop(columns=['Hour_Ending']) # drop unused columns
dt_y['datetime_CPT'] = pd.to_datetime(dt_y['Date']) + pd.to_timedelta(dt_y['Hour'], unit='h')
# Fill in missing values on 2022-12-01 01:00 (which is due to original recording error)
dt_y.loc[dt_y['datetime_CPT'].isnull(), 'datetime_CPT'] = pd.to_datetime('2022-12-01 01:00')
dt_y.loc[dt_y['Date'].isnull(), 'Date'] = '2022-12-01'
dt_y.loc[dt_y['Year'].isnull(), 'Year'] = 2022
dt_y.loc[dt_y['Month'].isnull(), 'Month'] = 12
dt_y.loc[dt_y['Day'].isnull(), 'Day'] = 1
dt_y.loc[dt_y['Hour'].isnull(), 'Hour'] = 1

# Create DST indicator and UTC time
# apply the function, causing two 2 AM to have the same DST value 
dt_y['DST'] = dt_y['datetime_CPT'].apply(lambda date: check_dst(date, dict_dst))
# locate the second 2 AM at the end date each year and alter DST to be 0
second_2am_ix = dt_y[dt_y.duplicated(subset=['datetime_CPT'], keep='first')].index
dt_y.loc[second_2am_ix, 'DST'] = 0
# create UTC time
dt_y['datetime_UTC'] = dt_y[['datetime_CPT', 'DST']].apply(lambda x: UTC_time(*x), axis=1)

#%% ######### Loop over regions #########
print('Merging data for each region...')
for region in regions[0:4]: # change it to 0:4 if have not done NC
    # Import region weather
    dt_w = pd.read_csv(os.path.join(md, 'data_clean', f'{region}_weather_2000-2024.csv'))
    dt_w['datetime'] = pd.to_datetime(dt_w['datetime'])
    dt_w = dt_w.rename(columns={'datetime': 'datetime_UTC'})
    # Merge with dt_w
    rname = dict_name[region] # find region's column name
    cols_drop = [item for item in cols_region if item != rname]
    dt_region = pd.merge(dt_y, dt_w, how='inner', on='datetime_UTC')
    dt_region = dt_region.drop(columns=cols_drop) # drop the other regions
    dt_region = dt_region.rename(columns={rname: 'load'}) # change rname to load
    dt_region = dt_region.dropna(subset=['load']) # drop na in load
    dt_region.to_csv(os.path.join(md, 'data_clean', f'{region}_main.csv'), index=False)
    print(f'{region} is done.')
# CAVEAT: DON'T open those files in Excel; otherwise, there would be unwanted conversion of format. Read them in R or pandas
# End of script