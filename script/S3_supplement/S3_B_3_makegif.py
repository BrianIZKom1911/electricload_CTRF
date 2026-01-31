# Convert PNG files into GIF in batches.
#%%
import os
from PIL import Image # needs Pillow

thisdir = os.path.dirname(__file__)
md = os.path.abspath(os.path.join(thisdir, '..', '..'))
in_dir = os.path.join(md, 'extraoutput')
out_dir = os.path.join(md, 'output_gif')
os.makedirs(out_dir, exist_ok=True)
#%% Define an no-output procedure
def DO_BATCH_CONVERT(file_path, file_list, save_path, name, duration):
    if not file_list:
        return

    file_list.sort()
    frames = []
    for pic in file_list:
        file = os.path.join(file_path, pic)
        img = Image.open(file)
        frames.append(img)
    
    # Save as GIF once all frames are collected
    frames[0].save(
        os.path.join(save_path, name),
        save_all=True,
        append_images=frames[1:],
        duration=duration, # miliseconds
        loop=0 # infinitely
    )   
    # Close images to clean up memory
    for f in frames: 
        f.close()

#%% 
regions = ['Coast', 'NC', 'SC', 'South']
seasons = ['winter', 'summer', 'SnF']
varnoms = ['prcp', 'rh', 'wsp', 'skc']

for region in regions:
    print("Processing ", region)
    short = region[:3]
    
    # I) TRF pic conversion
    filepath1 = os.path.join(in_dir, region)
    # List and filter TRF files for this region
    if os.path.exists(filepath1):
        trf_files = [f for f in os.listdir(filepath1) if f.startswith("trf_plot_bt_h") and f.endswith(".png")]
    # Skip if the result already exists
    name_trf = f"{region}_TRF24.gif"
    if os.path.exists(os.path.join(out_dir, name_trf)):
        print(f"Skipping {name_trf}, which already exists.")
    else:
        DO_BATCH_CONVERT(filepath1, trf_files, out_dir, name_trf, duration=300)
    
    # II) CTRF pic conversion
    for season in seasons:
        filepath2 = os.path.join(in_dir, region, season)
        if not os.path.exists(filepath2):
            continue
        # List all files in the season subdirectory
        allfiles = os.listdir(filepath2)

        # Look for "CTRF_<id>_h<hour>_<season>.png"
        # Map variables to their IDs (2, 3, 4, 5)
        for i, varnom in enumerate(varnoms, start=2):
            chead = f"CTRF_{i}_"
            cfes = [fe for fe in allfiles if fe.startswith(chead) and fe.endswith(".png")]
            # Skip if the result already exists
            name_wc = f"AB_{short}_{varnom}_{season}.gif"
            if os.path.exists(os.path.join(out_dir, name_wc)):
                print(f"Skipping {name_wc}, which already exists.")
            else:
                DO_BATCH_CONVERT(filepath2, cfes, out_dir, name_wc, duration=500)

print("Conversion complete.")
# End of script.
