from sklearn.decomposition import PCA
from pathlib import Path
import pickle
import prody
import collections
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.use("Agg")

def compute_PCA(input):
    pca = PCA(n_components=2, random_state = 12)
    transformed = pca.fit_transform(input)
    print(f"Explained variance: {pca.explained_variance_ratio_}")
    return transformed, pca

def pickle_obj(object, save_path):
    with open(save_path, 'wb') as handle:
        pickle.dump(object, handle, protocol=pickle.HIGHEST_PROTOCOL)

def load_dcd(dcd_paths, init_pdb_path, ref_pdb_path=None):
    if isinstance(dcd_paths, collections.abc.Iterable):
        traj = prody.Trajectory(str(dcd_paths[0]))
        if len(dcd_paths)>1:
            for dcd_path in dcd_paths[1:]:
                traj.addFile(str(dcd_path))
    atom_group = prody.parsePDB(str(init_pdb_path))
    protein_group = atom_group.select("protein and backbone")
    ref_struc = protein_group
    if ref_pdb_path:
        ref_struc = prody.parsePDB(str(ref_pdb_path)).select("protein and backbone")
    traj.link(atom_group)
    traj.setAtoms(protein_group)
    return create_custom_ens(traj, ref_struc.toAtomGroup())    

def create_custom_ens(inp, ref_struc):
    custom_ens = prody.Ensemble()
    if ref_struc:
        custom_ens.addCoordset(ref_struc)
        #select backbone
        ref = 0
    else:
        ref = None
    custom_ens.addCoordset(inp)
    custom_ens.setCoords(inp)
    if isinstance(inp, prody.AtomGroup):
        custom_ens.setAtoms(inp)
    else:
        custom_ens.setAtoms(inp.getAtoms())
    # http://prody.csb.pitt.edu/manual/reference/ensemble/ensemble.html#prody.ensemble.ensemble.Ensemble.superpose
    custom_ens.superpose(ref=ref)
    if ref_struc:
        custom_ens.delCoordset(0)
    flat = np.array([curr_entry.flatten() for curr_entry in list(custom_ens.iterCoordsets())])
    return flat
    
def load_pdb(pdb_path, ref_pdb_path=None):
    pdb = prody.parsePDB(str(pdb_path)).select("protein and backbone").toAtomGroup()
    ref_struc = None
    if ref_pdb_path:
        ref_struc = prody.parsePDB(str(ref_pdb_path)).select("protein and backbone").toAtomGroup()
    return create_custom_ens(pdb,ref_struc)
    
def make_plot(point_list, plt_title, save_path=None):
    print(f"{plt_title} has {len(point_list[0])} many points")
    for points in point_list[:-2]:
        plt.scatter(points[:,0], points[:,1], alpha=0.1, label="Simulation")
    plt.scatter(point_list[-2][:,0],point_list[-2][:,1], label="Prediction")
    plt.scatter(point_list[-1][:,0],point_list[-1][:,1], label="Starting Point")
    plt.title(plt_title)
    plt.legend()
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, format="pdf", bbox_inches='tight', transparent=True)
        plt.savefig(save_path.parent/"pca_plot.png", format="png", bbox_inches='tight', transparent=True)
    plt.show()
    plt.clf()

if __name__ == "__main__":
    for parent_path in Path("default_porter_simluations").glob("*"):
        identifier = parent_path.name
        print(identifier)
        if not (parent_path/f"run{5}"/"production_001.dcd").exists():
            continue
        coords = load_dcd([parent_path/f"run{it}"/"production_001.dcd" for it in range(1,6)], parent_path/"run1"/"protein_water_ions.pdb", parent_path/"run1"/"protein_water_ions.pdb")
        starting_point = load_pdb(parent_path/f"{identifier}_true_conf.pdb", parent_path/"run1"/"protein_water_ions.pdb")
        if len(list(parent_path.glob("best*")))>0:
            best_path = list(parent_path.glob("best*"))[0]
            print(best_path)
            prediction = load_pdb(best_path, parent_path/"run1"/"protein_water_ions.pdb")
            
        else:
            prediction = load_pdb(parent_path/f"{identifier}_prot_conf.pdb", parent_path/"run1"/"protein_water_ions.pdb")
        transformed, pca = compute_PCA(coords)
        pickle_obj(pca, parent_path/"pca.pkl")
        starting_pca = pca.transform(starting_point)
        prediction_pca = pca.transform(prediction)
        make_plot([transformed, prediction_pca, starting_pca], f"{identifier}", parent_path/"pca_plot.pdf")