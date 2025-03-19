import argparse
import csv
from pathlib import Path

def get_scores(output_path, max_similarity):
    different = []
    with open(output_path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            if max_similarity>float(row[2]):
                different.append(row[0])
    return different

def read_query_sequences(query_path):
    with open(query_path) as query_file:
        lines = query_file.readlines()
    entries = []
    curr_entry = ("", "")
    for line in lines:
        if line.startswith(">"):
            if len(curr_entry[1]) > 0:
                entries.append(curr_entry)
            curr_entry = (line, "")
        else:
            curr_entry[1] = curr_entry[1]+line
    if len(curr_entry[1]) > 0:
        entries.append(curr_entry)
    return entries

def keep_different_enough(query_sequences, different_enough):
    kept_sequences = []
    for entry in query_sequences:
        if entry[0] in different_enough:
            kept_sequences.append(entry[0])
            kept_sequences.append(entry[1])
    return kept_sequences


def write_output(output_sequences, output_path):
    with open(output_path, "w") as output_file:
        output_file.writelines(output_sequences)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--query", required=True)
    parser.add_argument("--scores", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--similarity", required=True)
    args = vars(parser.parse_args())
    different_enough = get_scores(Path(args["scores"]), float(args["similarity"]))
    query_sequences = read_query_sequences(Path(args["query"]))
    kept_sequences = keep_different_enough(query_sequences, different_enough)
    write_output(kept_sequences, Path(args["output"]))