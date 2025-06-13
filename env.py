import csv
import os
from openpyxl import Workbook
import pandas as pd
import subprocess
import sys

envdata = r'/netscratch/dep_coupland/grp_fulgione/siva/scripts/european_pop_969.xlsx'
sheet1 = 'SC'
envcode = 'VCF-BGI code'
envcolname = sys.argv[1]

env = pd.read_excel(envdata, sheet1)
env_col = env[[envcode, envcolname]]

spainall = r"/netscratch/dep_coupland/grp_fulgione/siva/scripts/spainall.txt"
scadall = r"/netscratch/dep_coupland/grp_fulgione/siva/scripts/scadall.txt"
alpsall = r"/netscratch/dep_coupland/grp_fulgione/siva/scripts/alpsall.txt"
allall = r"/netscratch/dep_coupland/grp_fulgione/siva/scripts/all.txt"

with open(allall, "r") as file:
    alllines = file.readlines()
alltext = [line.strip() for line in alllines]
with open(spainall, "r") as file1:
    spainlines = file1.read().splitlines()
spainall = [line.strip() for line in spainlines]
with open(scadall, "r") as file2:
    scadlines = file2.read().splitlines()
scadall = [line.strip() for line in scadlines]
with open(alpsall, "r") as file3:
    alpslines = file3.read().splitlines()
alpsall = [line.strip() for line in alpslines]

rows_df = pd.DataFrame(columns=env_col.columns)

for value in alltext:
    rows = env_col[env_col[envcode].astype(str) == str(value)]
    if not rows.empty:
        rows_df = pd.concat([rows_df, rows])
rows_df.to_excel("env_extract.xlsx", index = False)

