from Bio import BiopythonWarning
import warnings
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import json
from pathlib import Path
import argparse
import subprocess


DSSP_LETTERS =["-","H","B","E","G","I","T","S"]

def get_secstruc(pdb_path: str):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        structure = PDBParser().get_structure("pdb", pdb_path)
        if structure is None or len(structure) == 0:
            return None
        model = structure[0]     
        dssp = list(DSSP(model, pdb_path, dssp="mkdssp"))
    residues = [lis[1] for lis in dssp] 
    secstruc = [lis[2] for lis in dssp]
    count_dicts = [{f"{letter}_total": secstruc.count(letter),f"{letter}_rel": secstruc.count(letter)/len(secstruc) }  for letter in DSSP_LETTERS]
    out = {"pdb_path": str(pdb_path), "residues": residues, "secstruc":secstruc}
    for count_dict in count_dicts:
        out = out | count_dict
    return out

def cut_to_alignment(input_pdb, alignment_path, start, length):
    output_path = Path(str(alignment_path)[:-len(Path(alignment_path).suffix)]+".pdb")
    with open(alignment_path, "r") as file:
        before_string = (file.readlines()[1]).strip()[:start]
    amount_gaps = before_string.count("-")-1 #-1 for 0 indexing
    real_start = start-amount_gaps
    real_end = real_start+length-1
    command_string = f"pdb_selres -{real_start}:{real_end} {input_pdb} > {output_path}"
    subprocess.run(command_string, shell=True, check=True)
    return output_path

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Jannik")
    p.add_argument("--input_pdb", action="store", help="Path to input .pdb")
    p.add_argument("--output_json", action='store', help='Path to output .json')
    p.add_argument("--cut_start", action="store", help="Start of region of interest in full sequence.")
    p.add_argument("--cut_length", action="store", help="Length of region of interest in full sequence.")
    p.add_argument("--alignment_path", action="store", help="Path to alignment fasta with full sequence on second.")
    args = p.parse_args()
    secstruc_input = Path(args.input_pdb)
    if args.alignment_path is not None:
        alignment_path = args.alignment_path
        secstruc_input = cut_to_alignment(secstruc_input ,alignment_path, int(args.cut_start), int(args.cut_length) )
    print(secstruc_input)
    output_dict = get_secstruc(secstruc_input)
    with open(Path(args.output_json), 'w', encoding ='utf8') as json_file:
        json.dump(output_dict, json_file)
