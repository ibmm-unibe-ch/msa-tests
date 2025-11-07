import mdtraj as md
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib
import subprocess
from pathlib import Path
from heapq import nsmallest
from operator import itemgetter
from difflib import SequenceMatcher
from md_traj_utils import clean_traj

PAIRS = ['7ahlE4yhdG',
 '1repC2z9oA',
 '3tp2A5lj3M',
 '1g2cF5c6bF',
 '5ec5P3zxgB',
 '1uxmK2namA',
 '5aoeB5ly6B',
 '1ovaA1jtiB',
 '3gmhL2vfxL',
 '3m1bF3lowA',
 '3j7wB3j7vG',
 '2lqwA2bzyB',
 '5hmgB1htmB',
 '4j3oF2jmrA',
 '1miqB1qs8B',
 '4nc9C4n9wA',
 '3j97M1xtgB',
 '2naoF1iytA',
 '5c1vA5c1vB',
 '5jytA2qkeE',
 '2lejA2lv1A',
 '3jv6A1zk9A',
 '5fhcJ1eboE',
 '4wsgC1svfC',
 '1qomB1nocA',
 '5keqF1dzlA',
 '4y0mJ4xwsD',
 '2c1uC2c1vB',
 '4qhfA4qhhA',
 '4aanA4aalA',
 '1x0gA1x0gD',
 '4ae0A4ow6B',
 '3ifaA5et5A',
 '1h38D1qlnA',
 '1xjtA1xjuB',
 '3hdeA3hdfA',
 '3ejhA3m7pA',
 '1k0nA1rk4B',
 '1xntA3lqcA',
 '4gqcC4gqcB',
 '4dxtA4dxrA',
 '4rwnA4rwqB',
 '2hdmA2n54B',
 '3vo9B3vpaD',
 '2p3vA2p3vD',
 '3ewsB3g0hA',
 '2ce7C3kdsG',
 '4phqA2wcdX',
 '3t1pA1kctA',
 '5l35D5l35G',
 '5f3kA5f5rB',
 '4qdsA2qqjA',
 '5jzhA5jztG',
 '4pyiA4pyjA',
 '5ineA3mkoA',
 '2a73B3l5nB',
 '2k0qA2lelA',
 '3uyiA3v0tA',
 '4a5wB3t5oA',
 '1rkpA2h44A',
 '1ceeB2k42A',
 '3o44A1xezA',
 '3kuyA5c3iF',
 '2n0aD2kkwA',
 '4m4rA4w50B',
 '5ejbC1wp8C',
 '1mbyA4yypA',
 '5i2mA5i2sA',
 '2nxqB1jfkA',
 '3j9cA3q8fA',
 '4rmbA4rmbB',
 '3njqA2pbkB',
 '3zwgN4tsyD',
 '4hddA2lepA',
 '4zt0C4cmqB',
 '5k5gA2kb8A',
 '1nqdA1nqjB',
 '5fluE2uy7D',
 '3qy2A1qb3A',
 '2ougC2lclA',
 '4zrbC4zrbH',
 '1mnmA1mnmB',
 '2nntA2mwfA',
 '4jphB5hk5H',
 '4fu4C4g0dZ',
 '4b3oB3meeA',
 '4twaA4ydqB',
 '2axzA2grmB',
 '4o0pA4o01D',
 '4rr2D3l9qB']

def reres(input_path, output_path, start):
    subprocess_string = f"pdb_reres -{start} {input_path} > {output_path}"
    subprocess.run(subprocess_string, shell=True, check=True)
    return output_path

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
        #assert set([name for name, _ in curr]) == set(other), f"Different topologies supplied to get_sorted_coords of size: {len(curr)} vs. {len(other)}"
    curr.sort()
    return np.concatenate([coords for (_, coords) in curr], axis=1), curr_atoms

def make_same_topology(curr, reference, offset):
    reference_topo = reference.topology
    xyz=None
    for _, atom in enumerate(list(reference_topo.atoms)):
        atom_string = str(atom)
        atom_name = atom_string.split("-")[1]
        residue_identifier = atom_string.split("-")[0]
        atom_resid = int(residue_identifier[3:])+int(offset)
        if xyz is None:
            xyz = curr.atom_slice(curr.topology.select(f"name {atom_name} and residue {atom_resid}")).xyz
        else:
            xyz = np.concatenate([xyz,curr.atom_slice(curr.topology.select(f"name {atom_name} and residue {atom_resid}")).xyz], axis=1)
    return md.Trajectory(xyz, reference_topo)

def load_files(input_paths, reference_path):
    # this is way too over-engineered if both structures have the same PDB overview 
    if reference_path:
        reference = md.load_pdb(reference_path)
        old_ref_seq = "".join([str(elem)[:3] for elem in list(reference.topology.residues)])
        curr = md.load(input_paths[-1])
        curr_seq = "".join([str(elem)[:3] for elem in list(curr.topology.residues)])
        match = SequenceMatcher(None, old_ref_seq, curr_seq).find_longest_match()
        reference = reference.atom_slice(reference.topology.select(f"resid {match.a/3} to {(match.a+match.size)/3-1} and element != H and name != OXT and protein")).center_coordinates()
        ref_seq = "".join([str(elem)[:3] for elem in list(reference.topology.residues)])
        diff = 1
        if ref_seq in old_ref_seq:
            diff += old_ref_seq.index(ref_seq)/3
    inputs = []
    for input in input_paths:
        if str(input)[-3:] == "dcd":
            curr = md.load_dcd(input, reference_path)
            curr = clean_traj(curr)
            curr_seq = "".join([str(elem)[:3] for elem in list(curr.topology.residues)])
            match = SequenceMatcher(None, ref_seq, curr_seq).find_longest_match()
            curr = curr.atom_slice(curr.topology.select(f"resid {match.b/3} to {(match.b+match.size)/3}")) #/3 because all residue strings have length 3
            curr[0].save_pdb("traj.pdb")
            curr = make_same_topology(curr, reference,match.a/3)
            if curr.n_atoms != reference.n_atoms:
                print(f"problem {reference_path}")
                return None, None
            curr = curr.superpose(reference)#.center_coordinates()
            inputs.append(curr)
        else:
            curr = md.load(input)
            curr = clean_traj(curr)
            curr_seq = "".join([str(elem)[:3] for elem in list(curr.topology.residues)])
            match = SequenceMatcher(None, ref_seq, curr_seq).find_longest_match()
            if diff < 0:#4rmbA,1mbyA < 0?
                curr = curr.atom_slice(curr.topology.select(f"resid {match.b/3} to {(match.b+match.size)/3-1}"))            
            else:
                curr = curr.atom_slice(curr.topology.select(f"resid {match.b/3} to {(match.b+match.size)/3}"))
            curr = make_same_topology(curr, reference,int(match.b/3))
            if curr.n_atoms != reference.n_atoms:
                print(f"problem {reference_path}")
                return None, None
            curr = curr.superpose(reference)#.center_coordinates()
            inputs.append(curr)
    return md.join(inputs, check_topology=True)

def compute_PCA(coords, n_components=2):
    pca = PCA(n_components=n_components, random_state = 12)
    frames, atoms, dims = coords.shape
    transformed = pca.fit_transform(coords.reshape(frames, atoms * dims))
    print(f"Explained variance of training: {pca.explained_variance_ratio_}")
    return transformed, pca

def make_plot(point_list, plt_title,num_predictions=1,labels=None ,save_path=None, variance=None):
    if not (labels is None):
        for it, _ in enumerate(labels[0:-num_predictions]):
            plt.scatter(point_list[it,0], point_list[it,1], alpha=0.1, label=labels[it])
        for it, _ in enumerate(reversed(labels[-num_predictions:])):
            plt.scatter(point_list[-(1+it),0], point_list[-(1+it),1], alpha=0.1, label="Prediction")
        plt.scatter(point_list[0,0], point_list[0,1], label="Start")
    else:
        plt.scatter(point_list[0:-num_predictions,0], point_list[0:-num_predictions,1], alpha=0.01, label="Simulation")
        plt.scatter(point_list[-num_predictions:,0], point_list[-num_predictions:,1], label="Prediction")
        plt.scatter(point_list[0,0], point_list[0,1], label="Start")
    if variance is None:
        plt.xlabel(f"PC1")
        plt.ylabel(f"PC2")
    else:
        plt.xlabel(f"PC1 [variance {variance[0]*100:.2f}%]")
        plt.ylabel(f"PC2 [variance {variance[1]*100:.2f}%]")
    plt.title(plt_title)
    plt.legend()
    plt.tight_layout()
    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, format="pdf", bbox_inches='tight', transparent=True)
    plt.show()
    plt.clf()

def make_sns_plot(point_list, plt_title,num_predictions=1, save_path=None, variance=None):
    pca = pd.DataFrame(point_list[:-num_predictions-1], columns=["PC1", "PC2"])
    #sns.jointplot(data=pca, x=f"PC1", y=f"PC2", kind="kde", palette=["teal"], legend="Simulation", fill=True, alpha=0.8, height=3.33)
    sns.kdeplot(data=pca, x=f"PC1", y=f"PC2", legend="Simulation", fill=True, alpha=0.8, color="#808080")
    matplotlib.rcParams["axes.spines.right"] = False
    matplotlib.rcParams["axes.spines.top"] = False
    plt.rcParams['font.size'] = 20
    plt.setp(plt.gca().spines.values(), linewidth=3)
    if variance is None:
        plt.xlabel(f"PC1")
        plt.ylabel(f"PC2")
    else:
        plt.gca().set_xlabel(f"PC1 [variance {variance[0]*100:.2f}%]", fontsize=20)
        plt.gca().set_xlabel(f"PC2 [variance {variance[1]*100:.2f}%]", fontsize=20)
    plt.scatter(point_list[-num_predictions:,0], point_list[-num_predictions:,1], label="Prediction", c="#90EE90", s=60)
    plt.scatter(point_list[0,0], point_list[0,1], label="Experiment", c="#8080FF", s=60)
    plt.title(plt_title)
    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(save_path, dpi=600, transparent=True, bbox_inches='tight')
    plt.clf()

def find_neighbours(all_points, points_of_interest, num_neighbours):
    neigh = NearestNeighbors(n_neighbors=num_neighbours)
    neigh.fit(all_points)
    indices = []
    indices.append(neigh.kneighbors(points_of_interest, num_neighbours)[1])
    return indices

def get_pdb_from_traj(traj, index, save_path=None):
    temp = "temp.pdb"
    pdb = traj[index]   
    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        pdb.save_pdb(temp)
        subprocess_string = f"pdb_rplresname -HIP:HIS {temp}| pdb_rplresname -HIE:HIS | pdb_rplresname -NLE:LEU > {save_path}"
        subprocess.run(subprocess_string, shell=True, check=True)
    return pdb

def get_closest_rmsd(traj, ref_frame=-1, n=10):
    rmsds = md.rmsd(traj, traj, ref_frame)
    smallest_with_indices = nsmallest(n, enumerate(rmsds), key=itemgetter(1))
    return smallest_with_indices

def ost_score(model_path, reference_path, output_path):
    input_string = f"docker run --rm -v $(pwd):$(pwd) registry.scicore.unibas.ch/schwede/openstructure:latest compare-structures --model {model_path} --reference {reference_path} --output {output_path} --lddt --local-lddt --bb-lddt --bb-local-lddt --tm-score --rigid-scores --lddt-no-stereochecks"
    subprocess.run(input_string, shell=True, check=True)
