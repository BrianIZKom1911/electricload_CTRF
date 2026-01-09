# Generate GIF needs the files to be sorted by hour numerically, not lexically.
# Rename the pictures to fix the sorting problem permenantly.
#%%
import os
import re

thisdir = os.path.dirname(__file__)
in_dir = os.path.join(thisdir, 'gifsourceimages')
# Two patterns to match:
# 1) "trf_plot_bt_h<hour>.png"
pattern1 = re.compile(r"(trf_plot_bt_h)(\d+)(\.png)$") # \d means one digit
# 2) "CTRF_<id>_h<hour>_<xxx>.png"
pattern2 = re.compile(r"(CTRF_\d+_h)(\d+)(_.+\.png)") # . means any character
patterns = [pattern1, pattern2]

for dirpath, _, filenames in os.walk(in_dir):
    for name in filenames:
        for patt in patterns:
            m = patt.match(name)
            if m: 
                break
        
        if not m:
            continue
        prefix, hour, suffix = m.groups()
        new_name = f"{prefix}{int(hour):02d}{suffix}"
        old_path = os.path.join(dirpath, name)
        new_path = os.path.join(dirpath, new_name)

        if old_path != new_path:
            os.rename(old_path, new_path)
            print(f"Renamed successfully: {new_path}")
# END
