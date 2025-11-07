import numpy as np
import mdtraj as md
from pathlib import Path
from md_traj_utils import compute_PCA, unpickle_obj, pickle_obj, make_hdbscan, make_deshaw_plot, clean_traj_c_alpha, compute_all_rmsds, best_hummer_q
from scipy import stats
import matplotlib.pyplot as plt
import subprocess
import re
import seaborn as sns
import pandas as pd
import matplotlib.patches as  mpatches

MAIN_PRED_COLOUR = "#006400"
MINOR_PRED_COLOUR = "#90EE90"
MEAN_PRED_COLOUR = "#48A948"
MAIN_SIM_COLOUR = "#3941B3"
MINOR_SIM_COLOUR = "#3996B3"
MEAN_SIM_COLOUR = "#396cb3"
OTHER_COLOUR = "gray"
INTERESTING_DESHAW_PROTEINS = {"NTL9": {"folded_PDB": "2hba", "Abbreviation":"NTL9", "Length": 37,"Mutations":[11]},"Protein_G": {"folded_PDB": "1MI0", "Abbreviation":"NuG2","Start":5,"Mutations":[41]},"Villin": {"folded_PDB": "2F4K", "Abbreviation":"2F4K","Mutations":[23,26,28],"Experiment_mutations":[25]}, }
FRAMES_PER = 2000
TO_MICRO = 10000000
SIMULATION_SAMPLING = 100 
FOLDING_SAMPLING = 20 # on top of simulation subsampling
FRAMES_TO_AVERAGE = 20
FRAMES_TO_EXCLUDE = 5

plt.rcParams.update({
    'font.size': 12,
    'xtick.labelsize': 10,  # Size of x-axis tick labels
    'ytick.labelsize': 10,  # Size of y-axis tick labels
    'axes.spines.top' : False,
    'axes.spines.right' : False,
})

def get_frame_number(path):
    match = re.search(r"frame_(\d+)", str(path))
    return int(match.group(1)) if match else -1

def get_dcd_number(path):
    match = re.search(r"protein-(\d+).dcd", str(path))
    return int(match.group(1)) if match else -1

def find_bigger_cluster(clusters):
    mode = stats.mode([it for it in clusters if it>0]).mode
    other = int((mode-1)**2)
    return mode, other

def make_long_plot(simulation, all_rmsds, ref, decision_boundary, save_path, simulation_sampling=SIMULATION_SAMPLING):
    plt.rcParams.update({'axes.spines.right' : True,})    
    sim_qs = best_hummer_q(simulation, ref)
    sim_rmsds = md.rmsd(clean_traj_c_alpha(simulation), clean_traj_c_alpha(ref), 0)
    max_sim_frames = len(simulation)
    mean_rmsds = np.mean(sim_rmsds)
    colour_boundary = 2*np.std(sim_rmsds)
    simulation_time_in_micro = max_sim_frames*simulation_sampling*FRAMES_PER/TO_MICRO
    sim_x_values = np.linspace(0, simulation_time_in_micro, max_sim_frames)
    plt.scatter(sim_rmsds, sim_x_values, s=1)
    plt.savefig(save_path.parent/"test.pdf")
    plt.close()
    pred_x_values = np.linspace(0, simulation_time_in_micro, len(all_rmsds.keys()))
    _, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 6))
    sim_main_it = sorted([it for it,rmsd in enumerate(sim_rmsds) if abs(rmsd-mean_rmsds)<colour_boundary])
    sim_minor_it = sorted(list(set(range(len(sim_rmsds)))-set(sim_main_it)))

    pred_main_it = sorted(list(set([int(sim_main/max_sim_frames*len(all_rmsds.keys()) ) for sim_main in sim_main_it])))
    pred_minor_it = sorted(list(set(range(len(all_rmsds.keys())))-set(pred_main_it)))
    clean_neighbours = []
    missing_neighbours = []
    for curr_rmsd in all_rmsds.values():
        clean_neighbour = np.mean([rmsd[0] for rmsd in curr_rmsd if rmsd[0]<max_sim_frames][:FRAMES_TO_AVERAGE])/max_sim_frames*simulation_time_in_micro
        clean_neighbours.append(clean_neighbour)
        missing_neighbour = len([rmsd[0] for rmsd in curr_rmsd if rmsd[0] <max_sim_frames and rmsd[1]>decision_boundary][:FRAMES_TO_EXCLUDE])
        missing_neighbours.append(missing_neighbour)
    clean_neighbours = np.array(clean_neighbours)
    missing_neighbours = np.array(missing_neighbours)
    ax1.scatter(sim_x_values[sim_main_it], sim_rmsds[sim_main_it], c=MAIN_SIM_COLOUR, s=1)
    ax1.scatter(sim_x_values[sim_minor_it], sim_rmsds[sim_minor_it], c=MINOR_SIM_COLOUR, s=1)
    ax1.set_ylabel("RMSD (nm)", color=MEAN_SIM_COLOUR)
    ax1.tick_params(axis='y', labelcolor=MEAN_SIM_COLOUR)
    ax3 = ax1.twinx()
    #ax3.spines.right.set_position(("axes", 1.2))
    ax3.set_ylabel("Nearest frames", color=MEAN_PRED_COLOUR)
    ax3.set_ylim([0,simulation_time_in_micro])
    ax3.tick_params(axis='y', labelcolor=MEAN_PRED_COLOUR)
    ax3.scatter(pred_x_values[pred_main_it], clean_neighbours[pred_main_it], color=MAIN_PRED_COLOUR,s=1)
    ax3.scatter(pred_x_values[pred_minor_it], clean_neighbours[pred_minor_it], color=MINOR_PRED_COLOUR,s=1)
    
    ax2.scatter(sim_x_values[sim_main_it], sim_qs[sim_main_it], c=MAIN_SIM_COLOUR, s=1)
    ax2.scatter(sim_x_values[sim_minor_it], sim_qs[sim_minor_it], c=MINOR_SIM_COLOUR, s=1)
    ax2.set_ylabel("Q", color=MEAN_SIM_COLOUR)
    ax2.tick_params(axis='y', labelcolor=MEAN_SIM_COLOUR)
    ax2.set_ylim([0,1])
    ax2.set_xlabel("Time (μs)")
    ax4 = ax2.twinx()
    #ax4.spines.right.set_position(("axes", 1.2))
    ax4.set_ylabel("Nearest frames", color=MEAN_PRED_COLOUR)
    ax4.set_ylim([0,simulation_time_in_micro])
    ax4.tick_params(axis='y', labelcolor=MEAN_PRED_COLOUR)
    ax4.scatter(pred_x_values[pred_main_it], clean_neighbours[pred_main_it], color=MAIN_PRED_COLOUR, s=1)
    ax4.scatter(pred_x_values[pred_minor_it], clean_neighbours[pred_minor_it], color=MINOR_PRED_COLOUR, s=1)
    
    half_diff = (pred_x_values[1]-pred_x_values[0])/2
    for ax in (ax3, ax4):
        for pred_x_value, missing_neighbour in zip(pred_x_values, missing_neighbours):
            if missing_neighbour:
                ax.axvspan(pred_x_value-half_diff, pred_x_value+half_diff, color="gray", alpha=0.1)
    # Labels and legends
    plt.tight_layout()
    #plt.title(save_path.stem)
    plt.savefig(save_path, format="pdf", bbox_inches='tight', transparent=True)
    plt.close()

def make_long_plot_q_colour(simulation, all_rmsds, ref, decision_boundary, save_path, pred_traj, simulation_sampling=SIMULATION_SAMPLING):
    plt.rcParams.update({'axes.spines.right' : False,})    
    sim_qs = best_hummer_q(simulation, ref)
    pred_qs = best_hummer_q(pred_traj, ref)
    sim_rmsds = md.rmsd(clean_traj_c_alpha(simulation), clean_traj_c_alpha(ref), 0)
    max_sim_frames = len(simulation)
    #mean_rmsds = np.mean(sim_rmsds)
    #colour_boundary = 2*np.std(sim_rmsds)
    simulation_time_in_micro = max_sim_frames*simulation_sampling*FRAMES_PER/TO_MICRO
    sim_x_values = np.linspace(0, simulation_time_in_micro, max_sim_frames)
    plt.scatter(sim_rmsds, sim_x_values, s=1)
    plt.savefig(save_path.parent/"test.pdf")
    plt.close()
    pred_x_values = np.linspace(0, simulation_time_in_micro, len(all_rmsds.keys()))
    _, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 6))
    #sim_it = range(len(sim_rmsds))
    #sim_main_it = sorted([it for it,rmsd in enumerate(sim_rmsds) if abs(rmsd-mean_rmsds)<colour_boundary])
    #sim_minor_it = sorted(list(set(range(len(sim_rmsds)))-set(sim_main_it)))

    #pred_it = list(set([int(sim_main/max_sim_frames*len(all_rmsds.keys()) ) for sim_main in sim_main_it]))
    #pred_main_it = sorted(list(set([int(sim_main/max_sim_frames*len(all_rmsds.keys()) ) for sim_main in sim_main_it])))
    #pred_minor_it = sorted(list(set(range(len(all_rmsds.keys())))-set(pred_main_it)))
    clean_neighbours = []
    missing_neighbours = []
    for curr_rmsd in all_rmsds.values():
        clean_neighbour = np.mean([rmsd[0] for rmsd in curr_rmsd if rmsd[0]<max_sim_frames][:FRAMES_TO_AVERAGE])/max_sim_frames*simulation_time_in_micro
        clean_neighbours.append(clean_neighbour)
        missing_neighbour = len([rmsd[0] for rmsd in curr_rmsd if rmsd[0] <max_sim_frames and rmsd[1]>decision_boundary][:FRAMES_TO_EXCLUDE])
        missing_neighbours.append(missing_neighbour)
    clean_neighbours = np.array(clean_neighbours)
    missing_neighbours = np.array(missing_neighbours)
    ax1.scatter(sim_x_values, sim_rmsds, c=sim_qs, s=1, cmap="viridis")
    #ax1.scatter(sim_x_values[sim_minor_it], sim_rmsds[sim_minor_it], c=MINOR_SIM_COLOUR, s=1)
    ax1.set_ylabel("RMSD (nm)")#, color=MEAN_SIM_COLOUR)
    #ax1.tick_params(axis='y'), labelcolor=MEAN_SIM_COLOUR)
    #ax3 = ax1.twinx()
    #ax3.spines.right.set_position(("axes", 1.2))
    ax2.set_ylabel("Nearest frames")#, color=MEAN_PRED_COLOUR)
    ax2.set_ylim([0,simulation_time_in_micro])
    #ax3.tick_params(axis='y', labelcolor=MEAN_PRED_COLOUR)
    ax2.scatter(pred_x_values, clean_neighbours, c=pred_qs, s=1, cmap="viridis")
    #ax3.scatter(pred_x_values[pred_minor_it], clean_neighbours[pred_minor_it], color=MINOR_PRED_COLOUR,s=1)
    
    #ax2.scatter(sim_x_values[sim_main_it], sim_qs[sim_main_it], c=MAIN_SIM_COLOUR, s=1)
    #ax2.scatter(sim_x_values[sim_minor_it], sim_qs[sim_minor_it], c=MINOR_SIM_COLOUR, s=1)
    #ax2.set_ylabel("Q", color=MEAN_SIM_COLOUR)
    #ax2.tick_params(axis='y', labelcolor=MEAN_SIM_COLOUR)
    #ax2.set_ylim([0,1])
    ax2.set_xlabel("Time (μs)")
    #ax4 = ax2.twinx()
    ##ax4.spines.right.set_position(("axes", 1.2))
    #ax4.set_ylabel("Nearest frames", color=MEAN_PRED_COLOUR)
    #ax4.set_ylim([0,simulation_time_in_micro])
    #ax4.tick_params(axis='y', labelcolor=MEAN_PRED_COLOUR)
    #ax4.scatter(pred_x_values[pred_main_it], clean_neighbours[pred_main_it], color=MAIN_PRED_COLOUR, s=1)
    #ax4.scatter(pred_x_values[pred_minor_it], clean_neighbours[pred_minor_it], color=MINOR_PRED_COLOUR, s=1)
    
    half_diff = (pred_x_values[1]-pred_x_values[0])/2
    for ax in [ax2]: #ax1
        for pred_x_value, missing_neighbour in zip(pred_x_values, missing_neighbours):
            if missing_neighbour:
                ax.axvspan(pred_x_value-half_diff, pred_x_value+half_diff, color="gray", alpha=0.05)
    # Labels and legends
    plt.tight_layout()
    #plt.title(save_path.stem)
    plt.savefig(save_path, format="pdf", bbox_inches='tight', transparent=True)
    plt.close()


def download_pdb(pdb_code:str, out_path:Path):
    command = f'pdb_fetch {pdb_code} | pdb_selmodel -1 | pdb_selchain -A | pdb_delhetatm | pdb_delinsertion | pdb_reres -1 | pdb_tidy | grep ^ATOM | grep -E "ALA|ARG|ASN|ASP|CYS|GLU|GLN|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|SEC|PYL|HCY" > {out_path}'
    subprocess.run(command, shell=True, check=True)
    pass  

def make_sns_plot(point_list, plt_title,num_predictions=1, save_path=None, variance=None):
    plt.rcParams.update({'axes.spines.right' : False,})    
    pca = pd.DataFrame(point_list[:-num_predictions-1], columns=["PC1", "PC2"])
    #sns.jointplot(data=pca, x=f"PC1", y=f"PC2", kind="kde", palette=["teal"], legend="Simulation", fill=True, alpha=0.8, height=3.33)
    sns.kdeplot(data=pca, x=f"PC1", y=f"PC2", palette=["gray"],color="gray", legend="Simulation", fill=True, alpha=0.8)
    if variance is None:
        plt.xlabel(f"PC1")
        plt.ylabel(f"PC2")
    else:
        plt.xlabel(f"PC1 [variance {variance[0]*100:.2f}%]")
        plt.ylabel(f"PC2 [variance {variance[1]*100:.2f}%]")
    #plt.scatter(point_list[0:-num_predictions,0], point_list[0:-num_predictions,1], alpha=0.01, label="Simulation")
    plt.scatter(point_list[-num_predictions:,0], point_list[-num_predictions:,1], label="Prediction", c=MEAN_PRED_COLOUR, s=1)
    handles = [mpatches.Patch(facecolor="gray", label="Simulation"), mpatches.Patch(facecolor=MEAN_PRED_COLOUR, label="Prediction")]
    plt.legend(handles=handles)
    plt.title(plt_title)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, transparent=True)#, bbox_inches='tight')
    plt.show()
    plt.clf()

if __name__ == "__main__":
    for protein_name, protein_info in INTERESTING_DESHAW_PROTEINS.items():
        print(f"Working on {protein_name}")
        protein_stub = Path("/data/jgut/msa-tests/deshaw_ovchinnikov")/protein_name
        folded_path = protein_stub/"folded.pdb"
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
            #print("Experiment_mutation")
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
        #print("@@@experiment traj")
        #for i in experimental.topology.atoms:
        #    print(i)
        protein_id = f"{protein_info['Abbreviation']}-0"
        sim_path = Path(f"/data/jgut/msa-tests/DEShaw_simulations/DESRES-Trajectory_{protein_id}-protein/{protein_id}-protein")
        sim_files_to_load = sorted(list(sim_path.glob("*.dcd")), key=get_dcd_number)
        sim_pdb = sim_path/f"{protein_id}-protein.pdb"
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
        #print("@@@pred traj")
        #for i in pred_traj.topology.atoms:
        #    print(i)
        pred_traj = pred_traj.superpose(experimental).center_coordinates()
        #pred_pca_file = protein_stub/"pred_pca.pkl"
        #if not pred_pca_file.exists():
        #    pred_pca_coords, pred_pca = compute_PCA(clean_traj_c_alpha(pred_traj).xyz,2)
        #    pickle_obj((pred_pca_coords, pred_pca), pred_pca_file)
        #else:
        #    pred_pca_coords, pred_pca = unpickle_obj(pred_pca_file)
        #
        #pred_hdb_file = protein_stub/"pred_hdbscan.pkl"
        #make_deshaw_plot(pred_pca_coords,"Prediction test",save_path=protein_stub/"pred_test.pdf")
#
        #if not pred_hdb_file.exists():
        #    pred_clusters, pred_medoids = make_hdbscan(pred_pca_coords)
        #    pickle_obj((pred_clusters, pred_medoids), pred_hdb_file)
        #else:
        #    pred_clusters, pred_medoids = unpickle_obj(pred_hdb_file)
        #print(f"pred medoids {pred_medoids}")
        
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
        #print("@@@sim traj")
        #for i in sim_traj.topology.atoms:
        #    print(i)
        sim_traj = sim_traj.superpose(experimental).center_coordinates()
        #sim_pca_file = protein_stub/"sim_pca.pkl" 
        #if not sim_pca_file.exists():
        #    sim_pca_coords, sim_pca = compute_PCA(clean_traj_c_alpha(sim_traj).xyz,2)
        #    pickle_obj((sim_pca_coords, sim_pca), sim_pca_file)
        #else:
        #    sim_pca_coords, sim_pca = unpickle_obj(sim_pca_file)
        #make_deshaw_plot(sim_pca_coords,"Simulation test",save_path=protein_stub/"sim_test.pdf")
        #sim_hdb_file = protein_stub/"sim_hdbscan.pkl"
        #if not sim_hdb_file.exists():
        #    sim_clusters, sim_medoids = make_hdbscan(sim_pca_coords)
        #    pickle_obj((sim_clusters, sim_medoids), sim_hdb_file)
        #else:
        #    sim_clusters, sim_medoids = unpickle_obj(sim_hdb_file)       
        #print(f"sim medoids {sim_medoids}")
        #sim_mode, sim_other = find_bigger_cluster(sim_clusters)
        #sim_colours = {-1: OTHER_COLOUR, sim_mode:MAIN_SIM_COLOUR, sim_other:MINOR_SIM_COLOUR}        
        #pred_mode, pred_other = find_bigger_cluster(pred_clusters)
        #pred_colours = {-1: OTHER_COLOUR, pred_mode:MAIN_PRED_COLOUR, pred_other:MINOR_PRED_COLOUR}
        #sim_plot_medoids = [(*sim_medoids[sim_mode],MAIN_SIM_COLOUR),(*sim_medoids[sim_other],MINOR_SIM_COLOUR),
        #                    (*sim_pca.transform(pred_pca.inverse_transform([pred_medoids[pred_mode]]))[0],MAIN_PRED_COLOUR),(*sim_pca.transform(pred_pca.inverse_transform([pred_medoids[pred_other]]))[0],MINOR_PRED_COLOUR),]
        #pred_plot_medoids = [(*pred_medoids[pred_mode],MAIN_PRED_COLOUR),(*pred_medoids[pred_other],MINOR_PRED_COLOUR),
        #                    (*pred_pca.transform(sim_pca.inverse_transform([sim_medoids[sim_mode]]))[0],MAIN_SIM_COLOUR),(*pred_pca.transform(sim_pca.inverse_transform([sim_medoids[sim_other]]))[0],MINOR_SIM_COLOUR),]
#
        #sim_plot_file = protein_stub/"sim_hdbscan.pdf"
        #if not sim_plot_file.exists():
        #    make_deshaw_plot(sim_pca_coords, protein_name.replace("-", " "), medoids=sim_plot_medoids, labels=[sim_colours[it] for it in sim_clusters] ,save_path=sim_plot_file, variance=sim_pca.explained_variance_ratio_)
        #pred_plot_file = protein_stub/"pred_hdbscan.pdf"
        #if not pred_plot_file.exists():
        #    make_deshaw_plot(pred_pca_coords, protein_name.replace("-", " "), medoids=pred_plot_medoids, labels=[pred_colours[it] for it in pred_clusters] ,save_path=pred_plot_file, variance=pred_pca.explained_variance_ratio_)
        clean_sim_traj = clean_traj_c_alpha(sim_traj)
        clean_pred_traj = clean_traj_c_alpha(pred_traj)
        comb_traj = md.join([clean_sim_traj, clean_pred_traj])
        comb_pca_file = protein_stub/"comb_pca.pkl" 
        if not comb_pca_file.exists():
            comb_pca_coords, comb_pca = compute_PCA(comb_traj.xyz,2)
            pickle_obj((comb_pca_coords, comb_pca), comb_pca_file)
        else:
            comb_pca_coords, comb_pca = unpickle_obj(comb_pca_file)
        make_sns_plot(comb_pca_coords,"Combined PCA",len(clean_pred_traj), save_path=protein_stub/"combined_pca.pdf", variance=comb_pca.explained_variance_ratio_)

        #pred_rmsds_file = protein_stub/"pred_rmsds.pkl"
        #if not pred_rmsds_file.exists():
        #    pred_rmsds = compute_all_rmsds(clean_traj_c_alpha(sim_traj), clean_traj_c_alpha(pred_traj), FRAMES_TO_AVERAGE)
        #    pickle_obj(pred_rmsds, pred_rmsds_file)
        #else:
        #    pred_rmsds = unpickle_obj(pred_rmsds_file) 

        #make_long_plot(sim_traj, pred_rmsds, ref=experimental, decision_boundary=0.5,save_path=protein_stub/"long_experimental.pdf", simulation_sampling=SIMULATION_SAMPLING)
        
        #make_long_plot_q_colour(sim_traj, pred_rmsds, ref=experimental, decision_boundary=0.5,save_path=protein_stub/"long_experimental_q_colours.pdf", pred_traj=pred_traj, simulation_sampling=SIMULATION_SAMPLING)
