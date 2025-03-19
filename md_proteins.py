from md_utils import load_files, compute_PCA, make_sns_plot, find_neighbours, pickle_obj, unpickle_obj, get_pdb_from_traj, ost_score, PAIRS, get_closest_rmsd, reres
from md_traj_utils import make_deshaw_plot, make_hdbscan, run_single_pipeline
from pathlib import Path
import numpy as np

if __name__ == "__main__":
    for folder in [Path("/data/jgut/msa-tests/protein_only_dcd_03032025")]:
        for protein_path in list(folder.glob("*")):
            protein_name = protein_path.stem
            print(f"Working on {protein_name}.")
            protein_pdb = protein_path/f"{protein_name}.protein.pdb"
            protein_output_path = Path("/data/jgut/msa-tests/protein_output_03032025")/protein_name
            protein_output_path.mkdir(parents=True, exist_ok=True)
            protein_dcds = list(protein_path.glob("*.dcd"))
            both_name = [elem for elem in PAIRS if protein_name in elem][0]
            protein_predictions = list(Path(f"/data/jgut/msa-tests/aaa_porter_all_models/porter_all_models/{both_name}/{protein_name}_conf_dir/").glob("*6217.pdb"))
        #    protein_adjusted_predictions = [reres(protein_prediction, protein_path/protein_prediction.name, 0) for protein_prediction in protein_predictions]
        #    traj, topology = load_files(protein_dcds+protein_adjusted_predictions,protein_pdb)
        #    for it in range(1,len(protein_adjusted_predictions)+1):
        #        closest_rmsds = get_closest_rmsd(traj, topology, -it, n=10)
        #        for (closest_index, closest_rmsd) in closest_rmsds:
        #            output_name = protein_output_path/f"it_{it}_rmsd_{closest_rmsd:.3f}_index_{closest_index}.pdb"
        #            score_name = protein_output_path/f"it_{it}_rmsd_{closest_rmsd:.3f}_index_{closest_index}.json"
        #            get_pdb_from_traj(traj, closest_index, topology, output_name)
        #            ost_score(protein_adjusted_predictions[-len(protein_adjusted_predictions)+it], output_name, score_name)
            
        #    pca_file = protein_output_path/"pca.pkl"
        #    if not pca_file.exists():
        #        transformed, pca = compute_PCA(traj,2)
        #        pickle_obj((transformed, pca), pca_file)
        #    else:
        #        transformed, pca = unpickle_obj(pca_file) 
        #    sampled = transformed[::10]
        #    make_sns_plot(transformed, protein_name, len(protein_adjusted_predictions), protein_output_path/"pca.png", variance=pca.explained_variance_ratio_)
            #if len(protein_adjusted_predictions)==0:
            #    continue
            #get_pdb_from_traj(traj, 0, topology, protein_pdb)
    #neighbour_indices = find_neighbours(transformed, transformed[-len(protein_adjusted_predictions):],6)[0]
    #for pdb_index, pdb_indices in enumerate(neighbour_indices):
    #    for rank, neighbour_index in enumerate(pdb_indices):
    #        simulation_frame_path = protein_output_path/f"it-{pdb_index}_neigh-{rank}_index-{neighbour_index}.pdb"
    #        get_pdb_from_traj(traj, neighbour_index, topology, simulation_frame_path)
    #        score_path = protein_output_path/f"it-{pdb_index}_neigh-{rank}_index-{neighbour_index}_score.json"
    #        ost_score(protein_adjusted_predictions[-len(protein_adjusted_predictions)+pdb_index], simulation_frame_path, score_path)
            # Try other centers
            #if not plot_file.exists():
    #        plot_file = protein_output_path/"deshaw.pdf"
    #        hdb_scan_file = protein_output_path/"hdb.pkl"
    #        if not hdb_scan_file.exists():
    #            clusters, medoids = make_hdbscan(sampled)
    #            pickle_obj((clusters, medoids), hdb_scan_file)
    #        else:
    #            (clusters, medoids) = unpickle_obj(hdb_scan_file)
            #make_deshaw_plot(sampled, protein_name, medoids=medoids, labels=clusters ,save_path=plot_file, variance=pca.explained_variance_ratio_)
            analysis_path = protein_output_path/"strucs"
    #        analysis_path.mkdir(parents=True, exist_ok=True)
    #        for it, medoid in enumerate(medoids):
    #            ind = np.argmax([sample==medoid for sample in sampled])
    #        #    print(f"ind: {ind}")
    #            get_pdb_from_traj(traj, int(ind), topology, analysis_path/f"medoid_{it}_original.pdb")
    #            run_single_pipeline(analysis_path/f"medoid_{it}_original.pdb", analysis_path/f"medoid_{it}_inverse_folded_dir")
    #            num_medoids = len(medoids)
            for medoid_folder in analysis_path.glob("*inverse_folded_dir/"):
                for predicted_pdb in medoid_folder.glob("*.pdb"):
                    score_path = Path(str(predicted_pdb)[:-len("pdb")]+"_ost.json")
                    model_path = Path(str(medoid_folder)[:-len("inverse_folded_dir")]+"original.pdb")
                    ost_score(predicted_pdb, model_path, score_path)
