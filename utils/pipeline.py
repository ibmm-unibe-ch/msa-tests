import pickle
import subprocess
from sklearn.cluster import HDBSCAN
from sklearn.decomposition import PCA

def pickle_obj(object, save_path):
    with open(save_path, 'wb') as handle:
        pickle.dump(object, handle, protocol=pickle.HIGHEST_PROTOCOL)

def unpickle_obj(pickle_path):
    with open(pickle_path, "rb") as file:
        obj = pickle.load(file) 
    return obj

def run_single_pipeline(input_file, output_path):
    input_string = f"bash process_snapshot.sh {input_file} {output_path}"
    subprocess.run(input_string, shell=True, check=True)

def run_inverse_folding_check(input_file, output_path):
    input_string = f"bash process_inverse_folding_check.sh {input_file} {output_path}"
    subprocess.run(input_string, shell=True, check=True)

def run_single_seq(input_file, output_path):
    input_string = f'pdb_delhetatm {input_file} | pdb_delinsertion | pdb_reres -1 | pdb_tidy | grep ^ATOM | grep -E "ALA|ARG|ASN|ASP|CYS|GLU|GLN|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|SEC|PYL|HCY" | pdb_tofasta  > single_seq_input.a3m'
    subprocess.run(input_string, shell=True, check=True)
    input_string = f"colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds 1 --num-models 5 --num-recycle 3 --num-relax 0 single_seq_input.a3m {output_path}"
    subprocess.run(input_string, shell=True, check=True)

def run_normal_run(input_file, output_path):
    input_string = f'pdb_delhetatm {input_file} | pdb_delinsertion | pdb_reres -1 | pdb_tidy | grep ^ATOM | grep -E "ALA|ARG|ASN|ASP|CYS|GLU|GLN|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|SEC|PYL|HCY" | pdb_tofasta  > input.fasta'
    subprocess.run(input_string, shell=True, check=True)
    input_string = f"colabfold_batch --overwrite-existing-results --random-seed 6217 --num-seeds 1 --num-models 5 --num-recycle 3 --num-relax 0 input.fasta {output_path}"
    subprocess.run(input_string, shell=True, check=True)

def ost_score(model_path, reference_path, output_path):
    BESTNAME="*_unrelaxed_rank_001*seed_[0-9][0-9][0-9][0-9].pdb"
    if model_path.is_dir():
        model_path = list(model_path.glob(BESTNAME))[0]
    if reference_path.is_dir():
        reference_path = list(reference_path.glob(BESTNAME))[0]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    input_string = f"docker run --rm -v $(pwd):$(pwd) -v {model_path.parent}:{model_path.parent} -v {reference_path.parent}:{reference_path.parent} -v {output_path.parent}:{output_path.parent} registry.scicore.unibas.ch/schwede/openstructure:latest compare-structures --model {model_path} --reference {reference_path} --output {output_path} --lddt --local-lddt --bb-lddt --bb-local-lddt --tm-score --rigid-scores --lddt-no-stereochecks"
    subprocess.run(input_string, shell=True, check=True)

def compute_PCA(coords, n_components=2):
    pca = PCA(n_components=n_components, random_state = 12)
    frames, atoms, dims = coords.shape
    transformed = pca.fit_transform(coords.reshape(frames, atoms * dims))
    print(f"Explained variance of training: {pca.explained_variance_ratio_}")
    return transformed, pca

def make_hdbscan(transformed):
    hdbscan = HDBSCAN(min_cluster_size=30, cluster_selection_epsilon=0., leaf_size=40, n_jobs=32, store_centers="medoid")
    # Protein_G and Villin version
    #if len(transformed)>50:
    #    hdbscan = HDBSCAN(min_cluster_size=50, cluster_selection_epsilon=0.5, leaf_size=100, n_jobs=32, store_centers="medoid")
    #else:
    #    hdbscan = HDBSCAN(min_cluster_size=2, cluster_selection_epsilon=0., leaf_size=20, n_jobs=32, store_centers="medoid")
    clusters = hdbscan.fit_predict(transformed)
    medoids = hdbscan.medoids_
    return clusters, medoids