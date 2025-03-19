import pandas as pd
from pathlib import Path
import json 

def read_ost(ost_path:Path):
    with open(ost_path) as json_data:
        score_json = json.load(json_data)
    return_object = {}
    return_object["lddt"] = score_json.get("lddt",0)
    return_object["tm_score"] = score_json.get("tm_score",0)
    for name in ["model_clashes","model_bad_bonds","model_bad_angles", "reference_clashes","reference_bad_bonds","reference_bad_angles"]:
        return_object[name] = len(score_json.get(name, []))
    return return_object

if __name__ == "__main__":
    parent_path = Path("/data/jgut/msa-tests")
    df = pd.read_csv(parent_path/"porter_data.csv")
    scores = []
    for row in df.iterrows():
        struc_a = row[0]
        struc_b = row[1]
        case_name = struc_a+struc_b
        case_path = parent_path/"porter"/case_name
        curr_entry = {"case": case_name}
        for comparison in ["AA", "AB", "Adiff", "Aprot", "BA", "BB", "Bdiff", "Bprot", "ref_AB", "ref_AB", "ref_AR", "ref_BR", "Rprot"]:
            comparison_dict = read_ost(case_path/f"score_{comparison}.json")
            curr_entry = curr_entry|{f"{key}_{comparison}": value for key, value in comparison_dict.items()}
        scores.append(curr_entry)
    df.Dataframe(scores)