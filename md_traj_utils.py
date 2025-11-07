import mdtraj as md
from sklearn.neighbors import NearestNeighbors
from pathlib import Path
from utils import pickle_obj, unpickle_obj, run_single_seq, compute_PCA
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.patches as  mpatches
import matplotlib.lines as mlines
from heapq import nsmallest
from operator import itemgetter
from Bio import PDB
from itertools import combinations


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

def clean_traj_c_alpha(traj):
    return traj.atom_slice(traj.topology.select(f"name CA"))

def get_sorted_coords(traj, other_traj=None, c_alpha=False):
    curr_atoms = []
    curr = [] 
    for it, atom in enumerate(traj.topology.atoms):
        curr_atoms.append((atom.name, atom.element, atom.residue))
        curr.append((str(atom.residue.index).zfill(5)+str(atom),np.expand_dims(traj.xyz[:,it,:], axis=1)))
    if other_traj:
        other = []
        for atom in other_traj.topology.atoms:
            other.append(str(atom.residue.index).zfill(5)+str(atom))
        if not c_alpha:
            assert set([name for name, xyz in curr]) == set(other), f"Different topologies supplied to get_sorted_coords of size: {len(curr)} vs. {len(other)}"
    curr.sort()
    return np.concatenate([coords for (name, coords) in curr], axis=1), curr_atoms

def load_files(input_paths, reference_path, c_alpha=False):
    reference = None
    if reference_path and not str(reference_path).endswith(".prmtop"):
        reference = md.load(reference_path)
        reference = clean_traj(reference)
        if c_alpha:
            reference = clean_traj_c_alpha(reference)
    inputs = []
    for input in input_paths:
        if str(input).endswith("dcd"):
            curr = md.load_dcd(input, reference_path)
            curr = clean_traj(curr.atom_slice(curr.top.select("protein and not element H and not name OXT")))
            if c_alpha:
                curr = clean_traj_c_alpha(curr)
            if reference is None:
                reference = curr[0]
            curr = curr.superpose(reference).center_coordinates()
            inputs.append(curr)
        else:
            curr = md.load(input)
            curr = clean_traj(curr.atom_slice(curr.top.select("protein and not element H and not name OXT")))
            if c_alpha:
                curr = clean_traj_c_alpha(curr)
            if reference is None:
                reference = curr[0]#.atom_slice(curr.top.select("protein and not element H and not name OXT"))
            curr = curr.superpose(reference).center_coordinates()
            inputs.append(curr)
    combined, names = get_sorted_coords(inputs[0], c_alpha=c_alpha)
    if len(inputs) > 1:
        for other in inputs[1:]:
            combined = np.concatenate([combined,get_sorted_coords(other,inputs[0], c_alpha=c_alpha)[0]], axis=0)        
    return combined, names

def take_c_alpha(input_paths,reference_path):
    reference = md.load(reference_path)
    reference = clean_traj(reference)
    reference = clean_traj_c_alpha(reference)
    c_alpha_coords = []
    for input_path in input_paths:
        if str(input_path)[-3:] == "dcd":
            curr = md.load_dcd(input_path, reference_path)
        else:
            curr = md.load(input_path)
        curr = clean_traj(curr)
        curr = clean_traj_c_alpha(curr)
        curr = curr.superpose(reference).center_coordinates()
        c_alpha_coords.append(curr.xyz)
    return np.concatenate(c_alpha_coords)

def get_plddts(file_name):
    parser = PDB.PDBParser(PERMISSIVE=1)
    structure = parser.get_structure("1CUR",file_name)
    model = next(iter(structure))
    chain = next(iter(model))
    plddts = []
    for residue in chain:
        plddts.append(float(next(iter(residue)).get_bfactor()))
    return sum(plddts)/len(plddts)

def make_deshaw_plot(point_list, plt_title, medoids=None, labels=None ,save_path=None, variance=None):
    plt.rcParams.update({'axes.spines.right' : False,})    
    plt.scatter(point_list[:,0], point_list[:,1], alpha=0.1, c=labels)
    if not (medoids is None):
        if len(medoids[0])>2:
            plt.scatter([medoid[0] for medoid in medoids], [medoid[1] for medoid in medoids], c=[medoid[2] for medoid in medoids], marker="*", edgecolors="black", s=50)
        else:
            plt.scatter([medoid[0] for medoid in medoids], [medoid[1] for medoid in medoids], c=len(medoids), marker="*", edgecolors="black", s=50)
    if variance is None:
        plt.xlabel(f"PC1")
        plt.ylabel(f"PC2")    
    else:
        plt.xlabel(f"PC1 [variance {variance[0]*100:.2f}%]")
        plt.ylabel(f"PC2 [variance {variance[1]*100:.2f}%]")
    plt.title(plt_title)
    if not (medoids is None) and len(medoids[0])>2:
        black_star = mlines.Line2D([], [], color='black', marker='*', linestyle='None', label='Medoid')
        handles = [black_star,mpatches.Patch(facecolor=medoids[0][2], label="Major in"), mpatches.Patch(facecolor=medoids[1][2], label="Minor in"), mpatches.Patch(facecolor=medoids[2][2], label="Major out"), mpatches.Patch(facecolor=medoids[3][2], label="Minor out")]
        plt.legend(handles=handles)    
    plt.tight_layout()
    Path(save_path).parent.mkdir(parents=True, exist_ok=True)
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches='tight', transparent=True)
    plt.show()
    plt.clf()

def make_sns_plot(point_list, plt_title,num_predictions=1, save_path=None, variance=None):
    plt.rcParams.update({'axes.spines.right' : False,})    
    pca = pd.DataFrame(transformed[:-num_predictions-1], columns=["PC1", "PC2"])
    sns.kdeplot(data=pca, x=f"PC1", y=f"PC2", palette=["teal"], legend="Simulation", fill=True, alpha=0.8)
    if variance is None:
        plt.xlabel(f"PC1")
        plt.ylabel(f"PC2")
    else:
        plt.xlabel(f"PC1 [variance {variance[0]*100:.2f}%]")
        plt.ylabel(f"PC2 [variance {variance[1]*100:.2f}%]")
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
    if isinstance(xyz, np.ndarray):
        pdb = md.Trajectory(xyz, topology)
    else:
        pdb=xyz
    if save_path:
        pdb.save_pdb(save_path)
    return pdb

def best_hummer_q(traj, native):
    """Compute the fraction of native contacts according the definition from
    Best, Hummer and Eaton [1]
    
    Parameters
    ----------
    traj : md.Trajectory
        The trajectory to do the computation for
    native : md.Trajectory
        The 'native state'. This can be an entire trajecory, or just a single frame.
        Only the first conformation is used
        
    Returns
    -------
    q : np.array, shape=(len(traj),)
        The fraction of native contacts in each frame of `traj`
        
    References
    ----------
    ..[1] Best, Hummer, and Eaton, "Native contacts determine protein folding
          mechanisms in atomistic simulations" PNAS (2013)
    """
    
    BETA_CONST = 50  # 1/nm
    LAMBDA_CONST = 1.8
    NATIVE_CUTOFF = 0.45  # nanometers
    
    # get the indices of all of the heavy atoms
    heavy = native.topology.select_atom_indices('heavy')
    # get the pairs of heavy atoms which are farther than 3
    # residues apart
    heavy_pairs = np.array(
        [(i,j) for (i,j) in combinations(heavy, 2)
            if abs(native.topology.atom(i).residue.index - \
                   native.topology.atom(j).residue.index) > 3])
    
    # compute the distances between these pairs in the native state
    #print(f"heavy paris: {max(heavy_pairs)}")
    heavy_pairs_distances = md.compute_distances(native[0], heavy_pairs)[0]
    # and get the pairs s.t. the distance is less than NATIVE_CUTOFF
    native_contacts = heavy_pairs[heavy_pairs_distances < NATIVE_CUTOFF]
    print("Number of native contacts", len(native_contacts))
    
    # now compute these distances for the whole trajectory
    r = md.compute_distances(traj, native_contacts)
    # and recompute them for just the native state
    r0 = md.compute_distances(native[0], native_contacts)
    
    q = np.mean(1.0 / (1 + np.exp(BETA_CONST * (r - LAMBDA_CONST * r0))), axis=1)
    return q  

def compute_all_rmsds(simulation, predictions, FRAMES_TO_AVERAGE):
    combined_traj = md.join([simulation, predictions])
    all_rmsds = dict()
    for index in range(len(predictions)):
        rmsds = md.rmsd(combined_traj, predictions,index)
        all_rmsds[index] = nsmallest(len(predictions)+FRAMES_TO_AVERAGE, enumerate(rmsds), key=itemgetter(1))
    return all_rmsds

def rmsd(curr, reference):
    return np.sqrt(((curr-reference)**2).sum(-1).mean())

def compute_rmsds(np_traj, index, n, names1, names2=None, other2=None, only_c_alpha=False):
    if only_c_alpha and names1:
        c_alphas1 = [it for it, case in enumerate(names1) if case[0]=="CA"]
        np_traj1 = np_traj[0:other2,c_alphas1,:]
        if not other2 is None:
            c_alphas2 = [it for it, case in enumerate(names2) if case[0]=="CA"]
            np_traj2 = np_traj[other2:,c_alphas2,:]
            np_traj1.append(np_traj2)
        np_traj = np_traj1
    reference = np_traj[index]
    rmsds = [rmsd(curr, reference) for curr in np_traj]
    smallest_with_indices = nsmallest(n, enumerate(rmsds), key=itemgetter(1))
    return smallest_with_indices

if __name__ == "__main__":
    DESHAW_PATH = Path("/scratch/alphafold_database/DEShaw_simulations/")
    for protein_path in DESHAW_PATH.glob("DESRES*"):
        protein_string = str(protein_path.name)[len("DESRES-Trajectory_"):].strip()
        protein = protein_string[:-len("-protein")]
        output_path = Path("/data/jgut/msa-tests/deshaw")/protein
        output_path.mkdir(parents=True, exist_ok=True)
        pca_file = output_path/"pca.pkl"
        plot_file = Path("visualisations")/f"{protein}.pdf"
        traj_files = list((protein_path/protein_string).glob("*.dcd"))
        top_file = protein_path/protein_string/f"{protein_string}.pdb"
        traj, topology = load_files(traj_files,top_file)
        if not pca_file.exists():
            transformed, pca = compute_PCA(traj,2)
            pickle_obj((transformed, pca), pca_file)
        else:
            transformed, pca = unpickle_obj(pca_file) 
        sampled = transformed[::100]
        analysis_path = output_path/"strucs"
        analysis_path.mkdir(parents=True, exist_ok=True)
        run_single_seq(top_file, analysis_path/"single_seq")
        