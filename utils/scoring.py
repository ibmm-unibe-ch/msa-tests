import json
from pathlib import Path

from Bio import PDB


def get_plddts(file_name):
    parser = PDB.PDBParser(PERMISSIVE=1, QUIET=True)
    structure = parser.get_structure("1CUR",file_name)
    model = next(iter(structure))
    chain = next(iter(model))
    plddts = []
    for residue in chain:
        plddts.append(float(next(iter(residue)).get_bfactor()))
    return sum(plddts)/len(plddts)

def open_ost(ost_path:Path):
    if not ost_path.exists():
        return -1, -1, -1, -1, -1, -1 ,-1, -1
    with open(ost_path) as json_data:
        score_json = json.load(json_data)
    if score_json["status"] != "SUCCESS":
        return -1, -1, -1, -1, -1, -1 ,-1, -1
    lddt = score_json["lddt"] if "lddt" in score_json else 0
    bb_lddt = score_json["bb_lddt"] if "bb_lddt" in score_json else 0
    tm_score = score_json["tm_score"] if "tm_score" in score_json else 0
    inconsistent_residues = score_json["inconsistent_residues"] if "inconsistent_residues" in score_json else -1
    length = len(score_json["local_lddt"]) if "local_lddt" in score_json else -1
    model_bad_bonds = len(score_json["reference_bad_bonds"]) if "reference_bad_bonds" in score_json else -1
    model_bad_angles = len(score_json["reference_bad_angles"]) if "reference_bad_angles" in score_json else -1
    if not "reference" in score_json.keys():
        print(ost_path)
        print(score_json)
    ref_path = score_json["reference"].replace("/ibmm_data/","/data/").replace("/data/jgut/msa-tests/porter_all_models/", "/data/jgut/msa-tests/aaa_porter_all_models/porter_all_models/")
    plddt = get_plddts(ref_path)
    return lddt, bb_lddt, tm_score, inconsistent_residues, length, model_bad_bonds, model_bad_angles, plddt
