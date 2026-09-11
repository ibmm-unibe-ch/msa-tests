import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
from scipy import stats

MAIN_PRED_COLOUR = "#006400"
MINOR_PRED_COLOUR = "#90EE90"
MEAN_PRED_COLOUR = "#48A948"
MAIN_SIM_COLOUR = "#3941B3"
MINOR_SIM_COLOUR = "#3996B3"
MEAN_SIM_COLOUR = "#396cb3"
OTHER_COLOUR = "gray"
TO_MICRO = 10000000
FRAMES_TO_AVERAGE = 20
FRAMES_TO_EXCLUDE = 5

plt.rcParams.update({
    'font.size': 12,
    'xtick.labelsize': 10,  # Size of x-axis tick labels
    'ytick.labelsize': 10,  # Size of y-axis tick labels
    'axes.spines.top' : False,
    'axes.spines.right' : False,
})

def get_frame_number(path):
    match = re.search(r"frame_(\d+)", str(path))
    return int(match.group(1)) if match else -1

def get_dcd_number(path):
    match = re.search(r"protein-(\d+).dcd", str(path))
    return int(match.group(1)) if match else -1

def find_bigger_cluster(clusters):
    mode = stats.mode([it for it in clusters if it>0]).mode
    other = int((mode-1)**2)
    return mode, other

def download_pdb(pdb_code:str, out_path:Path):
    command = f'pdb_fetch {pdb_code} | pdb_selmodel -1 | pdb_selchain -A | pdb_delhetatm | pdb_delinsertion | pdb_reres -1 | pdb_tidy | grep ^ATOM | grep -E "ALA|ARG|ASN|ASP|CYS|GLU|GLN|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|SEC|PYL|HCY" > {out_path}'
    subprocess.run(command, shell=True, check=True)
