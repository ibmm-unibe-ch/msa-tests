from pathlib import Path
import argparse

def process(input_path, output_path, seq_id):
    with open(input_path, "r") as input_file:
        lines = [line.strip() for line in input_file.readlines() if line.strip()]
    out = lines[:2]
    inverse_folded_lines = lines[2:]
    for line_it, line in enumerate(inverse_folded_lines):
        if (line_it%2)==1:
            continue
        score = float(line.split("seq_recovery=")[1])
        if score <seq_id:
            out = out +[line, inverse_folded_lines[line_it+1]]
    with open(output_path, "w") as output_file:
        output_file.writelines([outline+"\n" for outline in out])

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--seqid", required=True)
    args = vars(parser.parse_args())
    input_path=Path(args["input"])
    output_path=Path(args["output"])
    seq_identity=float(args["seqid"])
    process(input_path, output_path, seq_identity)
