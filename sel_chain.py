import argparse
from pathlib import Path
import re
import logging
logging.basicConfig(level=logging.DEBUG, filename="find_pdb_offset.txt",filemode="a")

def get_correct_seq(fasta_path:Path, chain:str):
    with open(fasta_path, "r") as fasta_file:
        fasta_lines = fasta_file.readlines()
    iterator = 0
    while iterator < len(fasta_lines):
        middle_part = re.split(r"\|", fasta_lines[iterator])[1].strip()[1:]
        logging.debug(f"fasta_path: {fasta_path}, chain: {chain}, middle_part: {middle_part}")
        search_part=chain
        if "auth" in middle_part and not (str(fasta_path).endswith(("5ec5P_full.pdb.fasta", "3ejhA_full.pdb.fasta"))):
            search_part = f"auth {chain}"
        if search_part in middle_part:
            return fasta_lines[iterator+1].strip()
        else:
            iterator+=2
    return ""
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastapath", required=True)
    parser.add_argument("--chain", required=True)
    parser.add_argument("--outputpath", required=True)
    args = vars(parser.parse_args())
    fasta_path = Path(args["fastapath"])
    chain = str(args["chain"])
    output_path = Path(args["outputpath"])
    fasta_seq = get_correct_seq(fasta_path, chain)
    assert len(fasta_seq) > 0, f"real_fasta_seq not found in {fasta_path} for chain {chain}"
    with open(output_path, "w") as output_file:
        output_file.writelines([">1inp\n", fasta_seq, "\n"])
