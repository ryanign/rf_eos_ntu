"""
Ryan Pranantyo
EOS, April 2025
"""
import pandas as pd

df_garpan = pd.read_csv('./xyz_tmp/xyz_garispantai_rbi.csv')
df_buf100 = pd.read_csv('./xyz_tmp/xyz_garispantai_rbi_buffer100m.csv')

df_buf100['elev'] = -1.5

df = df_garpan # pd.concat([df_garpan, df_buf100], ignore_index=True)


with open('./xyz/softbreaklines_garpan_20250408.pts', 'w') as f:
    for _, group in df.groupby("segment"):
        group[['lon', 'lat', 'elev']].to_csv(f, header=False, index=False, sep=" ")
        f.write(">\n")
