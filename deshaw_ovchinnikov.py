import mdtraj as md
from pathlib import Path
from md_utils import get_pdb_from_traj
from md_traj_utils import run_single_seq,load_files,compute_rmsds,pickle_obj, clean_traj
from pathlib import Path
from utils import run_single_pipeline, run_normal_run, ost_score
import subprocess

INTERESTING_DESHAW_PROTEINS = {"Protein_G": {"folded_PDB": "1MI0", "Abbreviation":"NuG2","Start":5,"Mutations":[41]},"NTL9": {"folded_PDB": "2hba", "Abbreviation":"NTL9", "Length": 37,"Mutations":[11]},"Villin": {"folded_PDB": "2F4K", "Abbreviation":"2F4K","Mutations":[23,26,28],"Experiment_mutations":[25]}, }
BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"

def download_pdb(pdb_code:str, out_path:Path, length=None):
    if length:
        length_cutting = "pdb_selres -:{length} | "
    else:
        length_cutting="" 
    command = f"pdb_fetch {pdb_code} | pdb_selmodel -1 | pdb_selchain -A | pdb_delhetatm | pdb_delinsertion | pdb_reres -1 | {length_cutting}pdb_tidy | grep ^ATOM >{out_path}"
    subprocess.run(command, shell=True, check=True)
    pass 

if __name__ == "__main__":
    fine_subsampling = 100 # 80 frames per
    coarse_subsampling = 20 # on top of fine subsampling
    for protein_name, protein_info in INTERESTING_DESHAW_PROTEINS.items():
        print(f"Working on protein {protein_name}")
        protein_abb = protein_info["Abbreviation"]
        protein_stub = Path(f"/scratch/alphafold_database/DEShaw_simulations/DESRES-Trajectory_{protein_abb}-0-protein/{protein_abb}-0-protein/")
        protein_pdb = protein_stub/f"{protein_abb}-0-protein.pdb"
        rmsds_to_prev = []
        rmsds_to_base = []
        output_path = Path(f"/data/jgut/msa-tests/deshaw_ovchinnikov/{protein_name}")
        output_path.mkdir(parents=True, exist_ok=True)
        folded_path = output_path/"folded.pdb"
        if not (folded_path.exists()):
            download_pdb(protein_info["folded_PDB"], folded_path, protein_info.get("length", None))
        protein_dcds = sorted(list(protein_stub.glob("*.dcd")))
        traj = md.load(protein_dcds, top=protein_pdb)[::fine_subsampling]
        traj = traj.superpose(traj,0).center_coordinates()
        traj = clean_traj(traj)
        # single sequence baseline
        if not ((output_path/"single_seq").exists()):
            run_single_seq(folded_path, output_path/"single_seq")
            ost_score(output_path/"single_seq",folded_path, output_path/"single_seq_score.json")
        # normal mmseqs2 baseline
        if not ((output_path/"normal").exists()):
            run_normal_run(folded_path, output_path/"normal")
            ost_score(output_path/"normal",folded_path, output_path/"normal_score.json")

        # actual inverse folding
        inverse_folded_pdbs = []
        for it,_ in enumerate(traj[::coarse_subsampling]):
            frame_it =coarse_subsampling*it
            frame_folder = output_path/f"frame_{frame_it}"
            frame_path = output_path/f"frame_{frame_it}.pdb"
            if not (frame_path.exists()):
                get_pdb_from_traj(traj,frame_it,frame_path)
            ost_score(folded_path,frame_path,Path(str(frame_path)[:-3]+"json"))
            # inverse_fold-fold frame
            if not(frame_folder.exists()):
                run_single_pipeline(frame_path, frame_folder)
            inverse_folded_pdb = list(frame_folder.glob(BESTNAME))[0]
            inverse_folded_pdbs.append(inverse_folded_pdb)
            
        # compute RMSDs for normal sequence
        pickle_path = output_path/f"{protein_name}_rmsds.pkl"
        if not (pickle_path.exists() or (output_path/"{protein_id}_rmsds.pkl").exists()):
            files_to_load = protein_dcds+inverse_folded_pdbs
            coordinates, names = load_files(files_to_load, protein_pdb, c_alpha=True)
            len_inverse_folded_pdbs = len(inverse_folded_pdbs)
            output = dict()
            for it, _ in enumerate(inverse_folded_pdbs):
                inverse_folded_it = it-len_inverse_folded_pdbs
                closest_rmsds = compute_rmsds(coordinates,inverse_folded_it,len_inverse_folded_pdbs+1, names,only_c_alpha=True)
                print(f"For {it}, {protein_name}, the closest frames are: {closest_rmsds}")
                output[it]=closest_rmsds
            pickle_obj(output, pickle_path)
