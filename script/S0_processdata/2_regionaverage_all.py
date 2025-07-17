# Description: This script can be used to calculate the average of all regions--NC, SC, Coast, and South
#%%
import pandas as pd
import glob
import os
import numpy as np

md = 'D:/OneDrive - University of Missouri/transfer_desktop/MU/2025spring_submit2'

regions = ['NC', 'SC', 'Coast', 'South']
cols_w = ['temperature', 'dew_point_temperature', 'station_level_pressure', 'wind_speed', 'precipitation', 
          'relative_humidity', 'wet_bulb_temperature', 'sky_cover_1']
cols_keep = ['Station_ID', 'Station_name', 'Year', 'Month', 'Day', 'Hour', 'Minute', 'Latitude', 'Longitude'] + cols_w
cols_w[(len(cols_w)-1)] = 'skycover' # used for averaging after numerizing sky_cover
#%%
# Define function to convert sky_cover to numerical values
def skycover_numerize(x):
    if x == 'CLR:00':
        return 0
    elif x == 'FEW:01':
        return 1/8
    elif x == 'FEW:02':
        return 2/8
    elif x == 'SCT:03':
        return 3/8
    elif x == 'SCT:04':
        return 4/8
    elif x == 'BKN:05':
        return 5/8
    elif x == 'BKN:06':
        return 6/8
    elif x == 'BKN:07':
        return 7/8
    elif x == 'OVC:08':
        return 1
    elif x == 'VV:09': # sky obscured
        return 1
    elif x == 'X:10':
        return np.nan
    else:
        return np.nan
#%%
# Import population data
dt_pop = pd.read_csv(os.path.join(md, 'data_clean', 'population_2000_2023.csv')) # change to 2024 if 2024 population is available 
#%% ############ Loop. a ############
print('I. Processing station data for all counties...')
for region in regions[0:4]: # SC, Coast, South # Change to 0:4 if have not done NC
    inpath = f'D:/OneDrive - University of Missouri/transfer_desktop/MU/database/GHCNh_Texas/{region}'
    outpath = os.path.join(md, 'data_intermediate', region) # store the county data here
    if not os.path.exists(outpath):
        os.mkdir(outpath)

    # Import weather data - Prepare
    psv_filepaths = [os.path.join(inpath, file) for file in os.listdir(inpath) if file.endswith('.psv')]

    # Get station-county mapping from the list
    df_station = pd.read_excel(os.path.join(md, 'S0_processdata', 'county_station_list.xlsx'), sheet_name=region, usecols='A,C')
    df_station = df_station.rename(columns={'County': 'county', 'ID': 'Station_ID'})
    df_station[['county', 'Station_ID']] = df_station[['county', 'Station_ID']].apply(lambda x: x.str.strip())
    ## remove leading/trailing spaces
    # Identify files for counties with one or more stations
    county_stations = df_station.copy()
    county_stations['filepath'] = county_stations['Station_ID'].apply(lambda x: [path for path in psv_filepaths if x in path])
    county_stations = county_stations.groupby('county').agg(lambda x: x.explode().tolist()).reset_index()
    ## make stations and filepaths to a list for each county

    # Loop: Handle all counties at once
    for i in np.arange(county_stations.shape[0]):
        county = county_stations.loc[i, 'county']
        paths = county_stations.loc[i, 'filepath']
        dt_cty = pd.concat([pd.read_csv(path, sep='|', usecols=cols_keep) for path in paths], axis=0, ignore_index=True)
        dt_cty = dt_cty[(dt_cty['Year'] >= 2000) & (dt_cty['Year'] <= 2024)]

        # Clean weather data
        dt_cty['skycover'] = dt_cty['sky_cover_1'].apply(skycover_numerize) # numerize sky_cover
        dt_cty['precipitation'] = dt_cty['precipitation'].apply(lambda x: 0 if np.isnan(x) else x) # replace missing precipitation with 0
        dt_cty['relative_humidity'] = dt_cty['relative_humidity'].apply(lambda x: 109.9876 if x > 110 else x) # relative humidity can't be > 110
        dt_cty['dew_point_temperature'] = dt_cty['dew_point_temperature'].apply(lambda x: 34.9876 if x > 35 else x) # dew_point_temperature can hardly be > 35
        dt_cty['wet_bulb_temperature'] = dt_cty['wet_bulb_temperature'].apply(lambda x: 35.9876 if x > 36 else x) # wet_bulb_temperature can hardly be > 36
        dt_cty['temperature'] = dt_cty['temperature'].apply(lambda x: 47.9876 if x > 48 else x) # temperature can't be higher than 48.9
        # Average weather variables by datetime
        dt_cty['datetime'] = pd.to_datetime(dt_cty[['Year', 'Month', 'Day', 'Hour', 'Minute']]).dt.round('h')
        ## round the datetime to the nearest hour # previously: dt.floor.('h')
        dt_cty = dt_cty.groupby(['datetime'])[cols_w].mean().reset_index()
        dt_cty = dt_cty.dropna(subset=['temperature']) # drop rows with missing temperature because all models need temperature

        # Merge with population
        dt_cty['county'] = county
        dt_cty['year'] = dt_cty['datetime'].dt.year
        dt_cty = dt_cty.merge(dt_pop, on=['year', 'county'], how='left')
        dt_cty.to_csv(os.path.join(outpath, f'{county}_2000-2024.csv'), index=False)
    print(f'Counties in {region} are done.')
#%%
# Define function to calculate population-weighted average
def weighted_avg(df, xcols, w):
    result = {}

    for col in xcols:
        valid_mask = df[col].notna()  # Mask for non-missing values in this column
        wt_sum = (df[col] * df[w]).where(valid_mask, 0)  # Multiply only valid values
        sumofwts = df[w].where(valid_mask, 0)  # Sum only valid weights

        group_wt_sum = wt_sum.groupby(df["datetime"]).sum()
        group_sumofwts = sumofwts.groupby(df["datetime"]).sum()

        result[col] = group_wt_sum / group_sumofwts  # Compute weighted avg

    return pd.DataFrame(result).reset_index()
#%% ############ Loop. b ############
print('II. Processing county data for all regions...')
for region in regions[0:4]:
    # Import the county data (a)
    outpath = os.path.join(md, 'data_intermediate', region)
    county_files = glob.glob(os.path.join(outpath, '*.csv'))
    dt_counties = pd.concat([pd.read_csv(file) for file in county_files], ignore_index=True)
    # Let '2024-01-01 00:00' population be the same as 2023  
    dt_counties['population'] = dt_counties['population'].fillna(method='ffill') # forward filling missing values
    wt_avgs = weighted_avg(dt_counties, cols_w, 'population')
    wt_avgs.to_csv(os.path.join(md, 'data_clean', f'{region}_weather_2000-2024.csv'), index=False)
    print(f'{region} is done.')
# End of script
