import mdtraj as md
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from pathlib import Path
import pickle
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.patches as  mpatches
from sklearn.cluster import HDBSCAN
import subprocess

# Set font and tick parameters
plt.rcParams.update({
    #'font.family': 'Arial',
    'font.size': 12,
    'xtick.labelsize': 10,  # Size of x-axis tick labels
    'ytick.labelsize': 10,  # Size of y-axis tick labels
    'axes.spines.top' : False,
    'axes.spines.right' : False,
})

def clean_traj(traj):
    #get rid of hydrogens and last terminal
    atoms_to_keep = [a.index for a in traj.topology.atoms if str(a.element) != 'hydrogen' and str(a.name) != "OXT"]
    return traj.atom_slice(atoms_to_keep)

def get_sorted_coords(traj, other_traj=None):
    curr_atoms = []
    curr = [] 
    for it, atom in enumerate(traj.topology.atoms):
        curr_atoms.append((atom.name, atom.element, atom.residue))
        curr.append((str(atom.residue.index).zfill(5)+str(atom),np.expand_dims(traj.xyz[:,it,:], axis=1)))
    if other_traj:
        other = []
        for atom in other_traj.topology.atoms:
            other.append(str(atom.residue.index).zfill(5)+str(atom))
        assert set([name for name, xyz in curr]) == set(other), f"Different topologies supplied to get_sorted_coords of size: {len(curr)} vs. {len(other)}"
    curr.sort()
    return np.concatenate([coords for (name, coords) in curr], axis=1), curr_atoms

def load_files(input_paths, reference_path):
    if reference_path:
        reference = md.load_pdb(reference_path)
        reference = clean_traj(reference)
    inputs = []
    for input in input_paths:
        if str(input)[-3:] == "dcd":
            curr = md.load_dcd(input, reference_path)
            curr = clean_traj(curr)
            curr = curr.superpose(reference).center_coordinates()
            inputs.append(curr)
        else:
            curr = md.load(input)
            curr = clean_traj(curr)
            curr = curr.superpose(reference).center_coordinates()
            inputs.append(curr)
    combined, names = get_sorted_coords(inputs[0])
    if len(inputs) > 1:
        for other in inputs[1:]:
            combined = np.concatenate([combined,get_sorted_coords(other,inputs[0])[0]], axis=0)        
    return combined, names

def compute_PCA(coords, n_components=2):
    pca = PCA(n_components=n_components, random_state = 12)
    frames, atoms, dims = coords.shape
    transformed = pca.fit_transform(coords.reshape(frames, atoms * dims))
    print(f"Explained variance of training: {pca.explained_variance_ratio_}")
    return transformed, pca

def pickle_obj(object, save_path):
    with open(save_path, 'wb') as handle:
        pickle.dump(object, handle, protocol=pickle.HIGHEST_PROTOCOL)

def unpickle_obj(pickle_path):
    with open(pickle_path, "rb") as file:
        obj = pickle.load(file) 
    return obj

def make_hdbscan(transformed):
    hdbscan = HDBSCAN(min_cluster_size=50, cluster_selection_epsilon=0.5, leaf_size=100, n_jobs=32, store_centers="medoid")
    clusters = hdbscan.fit_predict(transformed)
    medoids = hdbscan.medoids_
    return clusters, medoids

def make_deshaw_plot(point_list, plt_title, medoids=None, labels=None ,save_path=None, variance=None):
    plt.scatter(point_list[:,0], point_list[:,1], alpha=0.1, c=labels)
    if not (labels is None):
        plt.scatter(medoids[:,0], medoids[:,1], c=range(len(medoids)), marker="*", edgecolors="red", s=50)
    if variance is None:
        plt.xlabel(f"PC1")
        plt.ylabel(f"PC2")    
    else:
        plt.xlabel(f"PC1 [variance {variance[0]*100:.2f}%]")
        plt.ylabel(f"PC2 [variance {variance[1]*100:.2f}%]")
    plt.title(plt_title)
    plt.tight_layout()
    Path(save_path).parent.mkdir(parents=True, exist_ok=True)
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches='tight', transparent=True)
    plt.show()
    plt.clf()

def make_sns_plot(point_list, plt_title,num_predictions=1, save_path=None, variance=None):
    pca = pd.DataFrame(transformed[:-num_predictions-1], columns=["PC1", "PC2"])
    #sns.jointplot(data=pca, x=f"PC1", y=f"PC2", kind="kde", palette=["teal"], legend="Simulation", fill=True, alpha=0.8, height=3.33)
    sns.kdeplot(data=pca, x=f"PC1", y=f"PC2", palette=["teal"], legend="Simulation", fill=True, alpha=0.8)
    if variance is None:
        plt.xlabel(f"PC1")
        plt.ylabel(f"PC2")
    else:
        plt.xlabel(f"PC1 [variance {variance[0]*100:.2f}%]")
        plt.ylabel(f"PC2 [variance {variance[1]*100:.2f}%]")
    #plt.scatter(point_list[0:-num_predictions,0], point_list[0:-num_predictions,1], alpha=0.01, label="Simulation")
    plt.scatter(point_list[-num_predictions,0], point_list[-num_predictions,1], label="Prediction", c="orange")
    plt.scatter(point_list[0,0], point_list[0,1], label="Start", c="green")
    handles = [mpatches.Patch(facecolor="green", label="Start"), mpatches.Patch(facecolor="blue", label="Simulation"), mpatches.Patch(facecolor="orange", label="Prediction")]
    plt.legend(handles=handles)
    plt.title(plt_title)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, transparent=True)#, bbox_inches='tight')
    plt.show()
    plt.clf()

def find_neighbours(all_points, points_of_interest, num_neighbours):
    neigh = NearestNeighbors(n_neighbors=num_neighbours)
    neigh.fit(all_points)
    indices = []
    indices.append(neigh.kneighbors(points_of_interest, num_neighbours)[1])
    return indices

def get_pdb_from_traj(traj, index, topology, save_path=None):
    xyz = traj[index]
    if isinstance(topology, list):
        topo = md.Topology()
        chain = topo.add_chain()
        residues = {}
        for serial, atom in enumerate(topology):
            name, element, res_id = atom
            if not(str(res_id) in residues):
                residue = topo.add_residue(str(res_id), chain)
                residues[str(res_id)] = residue
            residue = residues[str(res_id)]
            topo.add_atom(name, element, residue)
        topology = topo
    pdb = md.Trajectory(xyz, topology)
    if save_path:
        pdb.save_pdb(save_path)
    return pdb

def ost_score(model_path, reference_path, output_path):
    BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
    if model_path.is_dir():
        model_path = list(model_path.glob(BESTNAME))[0]
    if reference_path.is_dir():
        reference_path = list(reference_path.glob(BESTNAME))[0]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    input_string = f"docker run --rm -v $(pwd):$(pwd) registry.scicore.unibas.ch/schwede/openstructure:latest compare-structures --model {model_path} --reference {reference_path} --output {output_path} --lddt --local-lddt --bb-lddt --bb-local-lddt --tm-score --rigid-scores --lddt-no-stereochecks"
    subprocess.run(input_string, shell=True, check=True)

def run_single_pipeline(input_file, output_path):
    input_string = f"bash process_snapshot.sh {input_file} {output_path}"
    subprocess.run(input_string, shell=True, check=True)

def run_single_seq(input_file, output_path):
    input_string = f"pdb_tofasta {input_file} > single_seq_input.a3m"
    subprocess.run(input_string, shell=True, check=True)
    input_string = f"colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds 1 --num-models 5 --num-recycle 3 --num-relax 0 single_seq_input.a3m {output_path}"
    subprocess.run(input_string, shell=True, check=True)

if __name__ == "__main__":
    DESHAW_PATH = Path("/scratch/alphafold_database/DEShaw_simulations/")
    for protein_path in DESHAW_PATH.glob("DESRES*"):
        protein_string = str(protein_path.name)[len("DESRES-Trajectory_"):].strip()
        protein = protein_string[:-len("-protein")]
        output_path = Path("/data/jgut/msa-tests/deshaw")/protein
        output_path.mkdir(parents=True, exist_ok=True)
        #pca_file = output_path/"pca.pkl"
        #hdb_scan_file = output_path/"hdb.pkl"
        #plot_file = Path("visualisations")/f"{protein}.pdf"
        #traj_files = list((protein_path/protein_string).glob("*.dcd"))
        top_file = protein_path/protein_string/f"{protein_string}.pdb"
        #traj, topology = load_files(traj_files,top_file)
        #if not pca_file.exists():
        #    transformed, pca = compute_PCA(traj,2)
        #    pickle_obj((transformed, pca), pca_file)
        #else:
        #    transformed, pca = unpickle_obj(pca_file) 
        #sampled = transformed[::100]
        #if not hdb_scan_file.exists():
        #    clusters, medoids = make_hdbscan(sampled)
        #    pickle_obj((clusters, medoids), hdb_scan_file)
        #else:
        #    (clusters, medoids) = unpickle_obj(hdb_scan_file)
        #if not plot_file.exists():
        #    make_deshaw_plot(sampled, protein.replace("-", " "), medoids=medoids, labels=clusters ,save_path=plot_file, variance=pca.explained_variance_ratio_)
        analysis_path = output_path/"strucs"
        analysis_path.mkdir(parents=True, exist_ok=True)
        run_single_seq(top_file, analysis_path/"single_seq")
        #for it, medoid in enumerate(medoids):
        #    ind = np.where(sampled==medoid)[0]
        #    get_pdb_from_traj(traj, ind, topology, analysis_path/f"medoid_{it}_original.pdb")
        #    run_single_pipeline(analysis_path/f"medoid_{it}_original.pdb", analysis_path/f"medoid_{it}_inverse_folded_dir")
        #    num_medoids = len(medoids)
        #scores_path = output_path/"scores"
        #scores_path.mkdir(parents=True, exist_ok=True)
        #for i in range(num_medoids):
        #    for j in range(num_medoids):
        #        if i < j:
        #            ost_score(analysis_path/f"medoid_{i}_original.pdb",analysis_path/f"medoid_{j}_original.pdb", scores_path/f"o_{i}_o_{j}.json")
        #            ost_score(analysis_path/f"medoid_{i}_inverse_folded_dir",analysis_path/f"medoid_{j}_inverse_folded_dir", scores_path/f"i_{i}_i_{j}.json")
        #        ost_score(analysis_path/f"medoid_{i}_original.pdb",analysis_path/f"medoid_{j}_inverse_folded_dir", scores_path/f"o_{i}_i_{j}.json")
        