"""
Ryan Pranantyo
EOS, April 2025
"""
import pandas as pd

#df_garpan = pd.read_csv('./xyz_tmp/xyz_garispantai_rbi.csv')
#df_buf100 = pd.read_csv('./xyz_tmp/xyz_garispantai_rbi_buffer100m.csv')
#
#df_buf100['elev'] = -1.5

#df = df_garpan # pd.concat([df_garpan, df_buf100], ignore_index=True)

### load coastline reference
# here I put 5cm as the coastline reference, just to avoid missing island when no dtm available from DeltaDTM or FABDEM or others
df = pd.read_csv('/home/ignatius.pranantyo/DATA/working_deltadtm/vectors/xyz_garispantai_rbi__Sumatra-JawaBaliLombok-Sumbawa.csv')

### output filename
fout_gmt = '/home/ignatius.pranantyo/DATA/working_deltadtm/vectors/gmt__softbreaklines_coastline__20250424.pts'



with open(fout_gmt, 'w') as f:
    for _, group in df.groupby("segment"):
        group[['lon', 'lat', 'elev']].to_csv(f, header=False, index=False, sep=" ")
        f.write(">\n")
