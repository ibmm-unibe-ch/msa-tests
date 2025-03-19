import json
from pathlib import Path
import argparse
import shutil
import subprocess
import os

def get_scores(model, reference, output_file):
    current_working_directory = os.getcwd()
    subprocess.run(["podman","run","--rm","-v",f"{current_working_directory}:{current_working_directory}","openstructure:latest","compare-structures","--model",str(model),"--reference",str(reference),"--output",str(output_file),"--lddt","--local-lddt","--tm-score","--rigid-scores","--lddt-no-stereochecks"])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--parentpath", required=True)
    parser.add_argument("--confa", required=True)
    parser.add_argument("--confb", required=True)
    parser.add_argument("--output", required=True)
    args = vars(parser.parse_args())
    parent_path = Path(args["parentpath"]).parent
    confA = Path(args["confa"])
    confB = Path(args["confb"])
    diffusion_pdbs = []
    scores = []
    for diffusion_pdb in parent_path.glob("*.pdb"):
        diffusion_pdbs.append(diffusion_pdb)
        conf_A_score_path = diffusion_pdb.parent/f"{diffusion_pdb.name[:-4]}_A_score.json"
        conf_B_score_path = diffusion_pdb.parent/f"{diffusion_pdb.name[:-4]}_B_score.json"
        get_scores(diffusion_pdb, confA, conf_A_score_path)
        get_scores(diffusion_pdb, confB, conf_B_score_path)
        with open(conf_A_score_path, ) as json_data:
            score_A = json.load(json_data)
        with open(conf_B_score_path) as json_data:
            score_B = json.load(json_data)
        score = max(score_A["tm_score"], score_B["tm_score"])
        scores.append(score)
    min_tm_score = scores.index(min(scores))
    shutil.copy(diffusion_pdbs[min_tm_score], args["output"])