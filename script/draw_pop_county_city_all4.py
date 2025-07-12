# Make sure the packages below are installed
# In Anaconda Prompt: conda install geopandas
# Make sure you select the right interpreter like 'D:\Anaconda' if it's running in an IDE
import geopandas as gpd
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import os
import pandas as pd
# %% # Import GeoJSON
md = 'D:/OneDrive - University of Missouri/transfer_desktop/MU/2025spring_submit2'
os.chdir(md) # changeth to the target directory
counties = gpd.read_file('Texas Counties Cartographic Boundary Map_20250116.geojson')
counties['id'] = counties.index.astype(str)  # ensure the ID is a string to match GeoJSON
counties = counties.to_crs(epsg=4326)  # reproject to WGS84

# %%
file_path = os.path.join(md, 'region_county_station_list.xlsx')
# Initialize empty lists to store all data
all_choropleths = [None] * 4
all_scattergeos = [None] * 4
region_stats = {}  # store min/max for each region
# Define region colors
region_colors = {'NC': 'Oranges', 'SC': 'Greens', 'Coast': 'Blues', 'South': 'Purples'}
# Process each region
for i, (region, color_scale) in enumerate(region_colors.items()):
    r = i + 1 
    df_county_pop = pd.read_excel(file_path, sheet_name=f'Sheet{r}', usecols='A:B')
    df_cities = pd.read_excel(file_path, sheet_name=f'Sheet{r}', usecols='C:F')
    df_cities = df_cities.dropna()
    # Define custom population thresholds
    thresholds = [0, 25000, 100000, 800000, float('inf')]  # tiny, small, medium, big
    bin_labels = [2, 4, 8, 16]
    df_cities['dot_size'] = pd.cut(df_cities['city_population'], bins=thresholds, labels=bin_labels, right=False).astype(int)
    
    # Merge population data with county GeoJSON based on 'County' column
    df_counties_region = counties.merge(df_county_pop, left_on='name', right_on='County', how='left')
    # Store population stats for this region
    region_stats[region] = {
        'min': df_counties_region['population'].min(),
        'max': df_counties_region['population'].max()
    }

    # Create choropleth for this region
    choropleth = go.Choropleth(
        geojson=df_counties_region.__geo_interface__,
        locations=df_counties_region['County'],
        z=df_counties_region['population'],
        featureidkey="properties.name",
        colorscale=color_scale, # use the selected colorscale for the region
        marker_line_width=0.5,
        name=region
    )
    all_choropleths[i] = choropleth
    
    # Create scattergeo for cities in this region
    scattergeo = go.Scattergeo(
        lon=df_cities['longitude'],
        lat=df_cities['latitude'],
        text=df_cities.apply(
            lambda row: f"{row['City']}: {row['city_population']}", axis=1
        ),  # City name and population as hover text
        marker=dict(
            size=df_cities['dot_size'],
            color='orange',
            opacity=0.5,
            line=dict(width=1, color='darkred')
        ),
        hoverinfo='text',
        mode='markers',
        showlegend=False
    )
    all_scattergeos[i] = scattergeo

# %% More functionalities
# Assign distinct colorbars to each choropleth trace
regions = ['NC', 'SC', 'Coast', 'South']
for i, choropleth in enumerate(all_choropleths):
    region = regions[i]
    choropleth.update(
        colorbar=dict(
            title=f"{region} Population",  # region-specific title
            len=0.15,  # height of each colorbar
            y=0.85 - (i * 0.17),  # vertical position (stacked)
            x=1.02,  # horizontal position
            thickness=10,  # width
            outlinewidth=0
        ),
        zmin=region_stats[region]['min'], # region-specific min                             
        zmax=region_stats[region]['max'] # region-specific max
    )

# %% # Create combined figure 
fig = go.Figure(data=all_choropleths + all_scattergeos)
# Layout configuration
fig.update_layout(
    geo=dict(
        lakecolor='white',
        projection=dict(type='albers usa'), # use Albers projection for US
        showland=True,
        landcolor='lightgray',
        scope='usa',
        fitbounds = 'locations',
        visible=True
    ),
    margin=dict(r=150),  # ensure space for colorbars
)

# %% # Save the plot
fig.write_html('plot_pop_county_all4.html')
# End of script.