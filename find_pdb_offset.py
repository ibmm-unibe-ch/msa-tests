import argparse
from pathlib import Path
import re
import logging
logging.basicConfig(level=logging.DEBUG, filename="find_pdb_offset.txt",filemode="a")

def get_correct_seq(real_fasta_path:Path, chain:str):
    logging.debug(real_fasta_path)
    with open(real_fasta_path, "r") as real_fasta_file:
        real_fasta_lines = real_fasta_file.readlines()
    iterator = 0
    while iterator < len(real_fasta_lines):
        middle_part = re.split(r"\|", real_fasta_lines[iterator])[1]
        logging.debug(f"middle_part: {middle_part}")
        if chain in middle_part:
            return real_fasta_lines[iterator+1].strip()
        else:
            iterator+=2
    return ""
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--realfasta", required=True)
    parser.add_argument("--foundfasta", required=True)
    parser.add_argument("--realchain", required=True)
    parser.add_argument("--realstart", required=True)
    parser.add_argument("--length", required=True)
    args = vars(parser.parse_args())
    real_fasta_path = Path(args["realfasta"])
    found_fasta_path = Path(args["foundfasta"])
    real_chain = str(args["realchain"])
    start = int(args["realstart"])
    length = int(args["length"])
    real_fasta_seq = get_correct_seq(real_fasta_path, real_chain)
    logging.debug(f"real_fasta_seq: {real_fasta_seq}")
    assert len(real_fasta_seq) > 0, f"real_fasta_seq not found in {real_fasta_path} for chain {real_chain}"
    real_conf_sequence = real_fasta_seq[start:(start+length)]
    with open(found_fasta_path, "r") as found_fasta_file:
        found_fasta_lines = found_fasta_file.readlines()
        found_fasta_seq = "".join([elem.strip() for elem in  found_fasta_lines[1:]]) 
    logging.debug(real_conf_sequence)
    logging.debug(found_fasta_seq)
    splits = re.split(real_conf_sequence[:10], found_fasta_seq)
    logging.debug(splits)
    assert len(splits) > 1, f"no overlap found between {real_fasta_path} and {found_fasta_path} for chain {real_chain}"
    offset = len(splits[0])
    print(offset)
