from md_utils import load_files, compute_PCA, make_sns_plot, find_neighbours, get_pdb_from_traj, ost_score, PAIRS, get_closest_rmsd
from utils import pickle_obj, unpickle_obj
from pathlib import Path

hits = {
"1mbyA":"/data/jgut/msa-tests/protein_only_dcd_20052025/1mbyA/it_1_rmsd_0.183_index_25046.pdb",
"3o44A":"/data/jgut/msa-tests/protein_only_dcd_20052025/3o44A/it_1_rmsd_0.150_index_50091.pdb",
"4rmbA":"/data/jgut/msa-tests/protein_only_dcd_20052025_af3/4rmbA/it_3_rmsd_0.353_index_28817.pdb",
"1kctA":"/data/jgut/msa-tests/protein_only_dcd_20052025/1kctA/it_1_rmsd_0.330_index_50092.pdb",
"3lowA":"/data/jgut/msa-tests/protein_only_dcd_20052025/3lowA/it_1_rmsd_0.383_index_50092.pdb",
"4a5wB":"/data/jgut/msa-tests/protein_only_dcd_20052025_af3/4a5wB/it_2_rmsd_0.294_index_14978.pdb",
"1qlnA":"/data/jgut/msa-tests/protein_only_dcd_20052025/1qlnA/it_1_rmsd_0.334_index_26976.pdb",
"2n54B":"/data/jgut/msa-tests/protein_only_dcd_20052025/2n54B/it_5_rmsd_0.164_index_7690.pdb",
"4fu4C":"/data/jgut/msa-tests/protein_only_dcd_20052025/4fu4C/it_2_rmsd_0.313_index_74977.pdb",
"4g0dZ":"/data/jgut/msa-tests/protein_only_dcd_20052025/4g0dZ/it_1_rmsd_0.260_index_7106.pdb",
"4ow6B":"/data/jgut/msa-tests/protein_only_dcd_20052025/4ow6B/it_1_rmsd_0.159_index_11899.pdb",
"2k0qA":"/data/jgut/msa-tests/protein_only_dcd_20052025/2k0qA/it_2_rmsd_0.199_index_72148.pdb",
"5k5gA":"/data/jgut/msa-tests/protein_only_dcd_20052025/5k5gA/it_3_rmsd_0.140_index_64779.pdb",
"1xtgB":"/data/jgut/msa-tests/protein_only_dcd_20052025_af3/1xtgB/it_2_rmsd_0.385_index_50084.pdb",}

if __name__ == "__main__":
    for protein_name in hits.keys():
        if protein_name in [f.name for f in Path("/data/jgut/msa-tests/protein_only_dcd_25032025").iterdir() if f.is_dir()]:
            protein_path = Path("/data/jgut/msa-tests/protein_only_dcd_25032025")/protein_name
        elif protein_name in [f.name for f in Path("/data/jgut/msa-tests/protein_only_dcd_20052025").iterdir() if f.is_dir()]:
            protein_path = Path("/data/jgut/msa-tests/protein_only_dcd_20052025")/protein_name
        print(f"Working on {protein_name}.") 
        protein_output_path = Path("/data/jgut/msa-tests/protein_only_dcd_20052025_visualisations")/protein_name
        if protein_output_path.exists():
            continue
        protein_pdb = protein_path/f"{protein_name}.protein.pdb"
        if not(protein_pdb.exists()):
            continue
        if (protein_name in ["2axzA", "2ougC"]): #(protein_output_path/"pca.png").exists()
            continue
        protein_output_path.mkdir(parents=True, exist_ok=True)
        protein_dcds = list(protein_path.glob("*.dcd"))
        both_name = [elem for elem in PAIRS if protein_name in elem][0]
        #protein_predictions = list(Path(f"/data/jgut/msa-tests/aaa_porter_all_models/porter_all_models/{both_name}/{protein_name}_conf_dir_af3/").glob("*.pdb")) #.glob("*6217.pdb"))
        protein_adjusted_predictions = [Path(hits[protein_name])]
        #protein_adjusted_predictions = [reres(protein_prediction, protein_path/protein_prediction.name, 0) for protein_prediction in protein_predictions]
        joint_traj = load_files(protein_dcds+protein_adjusted_predictions,protein_pdb)
        if not (protein_output_path/"pca.svg").exists():
            for it in range(1,len(protein_adjusted_predictions)+1):
                closest_rmsds = get_closest_rmsd(joint_traj,-it, n=10)
                for (closest_index, closest_rmsd) in closest_rmsds:
                    output_name = protein_output_path/f"it_{it}_rmsd_{closest_rmsd:.3f}_index_{closest_index}.pdb"
                    score_name = protein_output_path/f"it_{it}_rmsd_{closest_rmsd:.3f}_index_{closest_index}.json"
                    get_pdb_from_traj(joint_traj, closest_index, output_name)
                    ost_score(protein_adjusted_predictions[-it], output_name, score_name)
        pca_file = protein_output_path/"pca.pkl"
        if not pca_file.exists():
            transformed, pca = compute_PCA(joint_traj.xyz,2)
            pickle_obj((transformed, pca), pca_file)
        else:
            transformed, pca = unpickle_obj(pca_file) 
        sampled = transformed[::10]
        make_sns_plot(transformed, protein_name, len(protein_adjusted_predictions), Path("visualisations")/f"{protein_name}_pca.svg", variance=pca.explained_variance_ratio_)
        if len(protein_adjusted_predictions)==0:
            continue
        get_pdb_from_traj(joint_traj, 0, protein_path/"start.pdb")
        neighbour_indices = find_neighbours(transformed, transformed[-len(protein_adjusted_predictions):],6)[0]
        for pdb_index, pdb_indices in enumerate(neighbour_indices):
            for rank, neighbour_index in enumerate(pdb_indices):
                simulation_frame_path = protein_output_path/f"it-{pdb_index}_neigh-{rank}_index-{neighbour_index}.pdb"
                get_pdb_from_traj(joint_traj, int(neighbour_index), simulation_frame_path)
                score_path = protein_output_path/f"it-{pdb_index}_neigh-{rank}_index-{neighbour_index}_score.json"
                ost_score(protein_adjusted_predictions[-len(protein_adjusted_predictions)+pdb_index], simulation_frame_path, score_path)
