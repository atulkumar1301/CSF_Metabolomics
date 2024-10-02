from pathlib import Path
from pathlib import PurePath
import sys
import pandas as pd
data = pd.read_csv ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Metabolites/Full_Data_Metabolites.txt", delimiter="\t")
data = data[data.columns[:28]]
source_dir = Path ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Zscore_Data/")
files = source_dir.iterdir()
files = source_dir.glob('*.txt')
for file_names in files:
    path_split = PurePath(file_names).parts
    filename = path_split [9].split (".")
    new_col_name = "Average_" + str (filename [0])
    df = pd.read_csv(file_names, delimiter="\t")
    df = df [["Average"]]
    df.rename(columns={'Average': new_col_name}, inplace = True)
    data = data.join (df)
data.to_csv ("/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Full_Data_Zscore_Average.txt", sep="\t", index=False)
