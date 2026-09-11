import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np
import mdtraj as md
from utils.md_traj_utils import compute_PCA, unpickle_obj, pickle_obj, make_deshaw_plot, clean_traj_c_alpha, compute_all_rmsds, best_hummer_q
from utils.deshaw_common import (
    MAIN_PRED_COLOUR, MINOR_PRED_COLOUR, MEAN_PRED_COLOUR,
    MAIN_SIM_COLOUR, MINOR_SIM_COLOUR, MEAN_SIM_COLOUR, OTHER_COLOUR,
    TO_MICRO, FRAMES_TO_AVERAGE, FRAMES_TO_EXCLUDE,
    get_frame_number, get_dcd_number, find_bigger_cluster, download_pdb,
    make_long_plot, make_long_plot_q_colour, make_sns_plot,
)
import matplotlib.pyplot as plt

INTERESTING_DESHAW_PROTEINS = {"Vilin": {"folded_PDB": "2F4K", "Abbreviation":"vilin",},"NTL9": {"folded_PDB": "2hba", "Abbreviation":"NTL9",}, "Protein_G": {"folded_PDB": "1MI0", "Abbreviation":"NuG2"}, }
FRAMES_PER = 400
SIMULATION_SAMPLING = 5
FOLDING_SAMPLING = 1 # on top of simulation subsampling


if __name__ == "__main__":
    for protein_name, protein_info in INTERESTING_DESHAW_PROTEINS.items():
        print(f"Working on {protein_name}")
        protein_abb = protein_info["Abbreviation"].lower()
        protein_stub = Path("/data/jgut/msa-tests/deshaw_ovchinnikov")/f"{protein_name}_unfolding"
        folded_path = Path(f"/scratch/tlemmin/desres/{protein_abb}/{protein_abb}.pdb")
        if not folded_path.exists():
            download_pdb(protein_info["folded_PDB"], folded_path)
        experimental = md.load_pdb(folded_path)
        experimental = experimental.atom_slice(experimental.topology.select('chainid == 0 and protein and symbol != H and name != OXT'))
        start = 0
        if "Start" in protein_info:
            start = protein_info["Start"]
            experimental = experimental.atom_slice(experimental.topology.select(f'resid {start} to 100000000'))
        if "Length" in protein_info:
            experimental = experimental.atom_slice(experimental.topology.select(f'resid 0 to {protein_info["Length"]}'))
        if "Experiment_mutations" in protein_info:
            mutations = protein_info["Experiment_mutations"] 
            selection_string = f'resid != {mutations[0]-start}'
            experimental = experimental.atom_slice(experimental.topology.select(selection_string))
            for it, further_mutation in enumerate(mutations[1:]):
                selection_string = f'resid != {further_mutation-start}'
                experimental = experimental.atom_slice(experimental.topology.select(selection_string)) 
        elif "Mutations" in protein_info:
            mutations = protein_info["Mutations"] 
            selection_string = f'resid != {mutations[0]-start}'
            experimental = experimental.atom_slice(experimental.topology.select(selection_string))
            for it, further_mutation in enumerate(mutations[1:]):
                selection_string = f'resid != {further_mutation-start}'
                experimental = experimental.atom_slice(experimental.topology.select(selection_string)) 
                experimental = experimental.center_coordinates()
        protein_id = f"{protein_info['Abbreviation'].lower()}-0"
        sim_path = Path(f"/scratch/tlemmin/desres/dynamics_{protein_abb}_2")
        sim_files_to_load = sorted(list(sim_path.glob("production*.dcd")), key=get_dcd_number)
        sim_pdb = sim_path/"input.prmtop"
        # load and prepare predictions
        pred_files_to_load = sorted(list(protein_stub.glob("frame_*/best.pdb")), key=get_frame_number)
        pred_traj = md.load(pred_files_to_load)
        pred_traj._unitcell_lengths = np.asarray([[1.,1.,1.]]*len(pred_traj)) 
        pred_traj._unitcell_angles  = np.asarray([[90.,90.,90.]]*len(pred_traj))
        if "Length" in protein_info:
            pred_traj = pred_traj.atom_slice(pred_traj.topology.select(f'resid 0 to {protein_info["Length"]}'))
        if "Mutations" in protein_info:
            mutations = protein_info["Mutations"] 
            selection_string = f'resid != {mutations[0]-start}'
            pred_traj = pred_traj.atom_slice(pred_traj.topology.select(selection_string))
            for it, further_mutation in enumerate(mutations[1:]):
                selection_string = f'resid != {further_mutation-start-it-1}'
                pred_traj = pred_traj.atom_slice(pred_traj.topology.select(selection_string))            
        pred_traj = pred_traj.superpose(experimental).center_coordinates()
        pred_pca_file = protein_stub/"pred_pca.pkl"
        if not pred_pca_file.exists():
            pred_pca_coords, pred_pca = compute_PCA(clean_traj_c_alpha(pred_traj).xyz,2)
            pickle_obj((pred_pca_coords, pred_pca), pred_pca_file)
        else:
            pred_pca_coords, pred_pca = unpickle_obj(pred_pca_file)        
        # load and prepare simulation
        sim_traj = md.load(sim_files_to_load, top=sim_pdb, stride=SIMULATION_SAMPLING)
        sim_traj = sim_traj.atom_slice(sim_traj.topology.select('chainid == 0 and protein and symbol != H and name != OXT'))
        if "Length" in protein_info:
            sim_traj = sim_traj.atom_slice(sim_traj.topology.select(f'resid 0 to {protein_info["Length"]}'))
        if "Mutations" in protein_info:
            mutations = protein_info["Mutations"] 
            selection_string = f'resid != {mutations[0]-start}'
            sim_traj = sim_traj.atom_slice(sim_traj.topology.select(selection_string))
            for it, further_mutation in enumerate(mutations[1:]):
                selection_string = f'resid != {further_mutation-start-it-1}'
                sim_traj = sim_traj.atom_slice(sim_traj.topology.select(selection_string))            
        sim_traj = sim_traj.superpose(experimental).center_coordinates()
        sim_pca_file = protein_stub/"sim_pca.pkl" 
        if not sim_pca_file.exists():
            sim_pca_coords, sim_pca = compute_PCA(clean_traj_c_alpha(sim_traj).xyz,2)
            pickle_obj((sim_pca_coords, sim_pca), sim_pca_file)
        else:
            sim_pca_coords, sim_pca = unpickle_obj(sim_pca_file)
        make_deshaw_plot(sim_pca_coords,"Simulation test",save_path=protein_stub/"sim_test.pdf")
        clean_sim_traj = clean_traj_c_alpha(sim_traj)
        clean_pred_traj = clean_traj_c_alpha(pred_traj)
        comb_traj = md.join([clean_sim_traj, clean_pred_traj])
        comb_pca_file = protein_stub/"comb_pca.pkl" 
        if not comb_pca_file.exists():
            comb_pca_coords, comb_pca = compute_PCA(comb_traj.xyz,2)
            pickle_obj((comb_pca_coords, comb_pca), comb_pca_file)
        else:
            comb_pca_coords, comb_pca = unpickle_obj(comb_pca_file)
        make_sns_plot(comb_pca_coords,"Combined PCA",len(clean_pred_traj), save_path=protein_stub/"combined_pca_thresh_001.pdf", variance=comb_pca.explained_variance_ratio_, kde_thresh=0.001, kde_levels=20)

        pred_rmsds_file = protein_stub/"pred_rmsds.pkl"
        if not pred_rmsds_file.exists():
            pred_rmsds = compute_all_rmsds(clean_traj_c_alpha(sim_traj), clean_traj_c_alpha(pred_traj), FRAMES_TO_AVERAGE)
            pickle_obj(pred_rmsds, pred_rmsds_file)
        else:
            pred_rmsds = unpickle_obj(pred_rmsds_file) 

        #make_long_plot(sim_traj, pred_rmsds, ref=experimental, decision_boundary=0.5,save_path=protein_stub/"long_experimental.pdf", simulation_sampling=SIMULATION_SAMPLING, frames_per=FRAMES_PER, save_debug_scatter=True, missing_comparison_gt=False, verbose=True)
        make_long_plot_q_colour(sim_traj, pred_rmsds, ref=experimental, decision_boundary=0.5,save_path=protein_stub/"long_experimental_q_colours.pdf", pred_traj=pred_traj, simulation_sampling=SIMULATION_SAMPLING, frames_per=FRAMES_PER, save_debug_scatter=True, missing_comparison_gt=False, gate_by_count=True, verbose=True)
