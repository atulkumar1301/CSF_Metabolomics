from pathlib import Path
from pathlib import PurePath
import sys
import subprocess
source_dir = Path('/Volumes/ATUL_6TB/Work/Projects/CSF_Metabolomics/Analyses_2/Sub-Pathway/Data_Sub_Pathway/')
files = source_dir.iterdir()
files = source_dir.glob('*.txt')
for file_names in files:
    path_split = PurePath(file_names).parts
    filename = path_split [9]
    R_code = subprocess.call (["/Volumes/ATUL_6TB/Tools/R_Codes/Zscore_Calculation.R", str(filename)])
