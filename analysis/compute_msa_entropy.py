from prody import parseMSA, calcShannonEntropy, calcMSAOccupancy, calcMeff
from scipy.stats import spearmanr
from pathlib import Path
import numpy as np
import subprocess
from typing import List, Dict, Tuple, Sequence, Iterable
from collections import Counter
import math

sequence_metrics = ["identity", "entropy", "gap_fraction", "blosum_wt", "blosum_pairwise"]
type_metrics = ["type_freq", "type_consensus_freq", "type_entropy"]

#stolen from https://github.com/prody/ProDy/blob/main/prody/sequence/evo_from_msa.py
# Standard 4-group classification for type conservation
# H: Hydrophobic, P: Polar Neutral, N: Negative/Acidic, B: Positive/Basic, X: Special/Unknown/Gap
AA_TYPES = {
    "A": "H",
    "I": "H",
    "L": "H",
    "M": "H",
    "F": "H",
    "P": "H",
    "W": "H",
    "V": "H",
    "G": "H",  # Hydrophobic (H)
    "S": "P",
    "T": "P",
    "C": "P",
    "N": "P",
    "Q": "P",
    "Y": "P",  # Polar Neutral (P)
    "D": "N",
    "E": "N",  # Negative/Acidic (N)
    "K": "B",
    "R": "B",
    "H": "B",  # Positive/Basic (B)
    "U": "X",
    "O": "X",
    "X": "X",
    "-": "X",  # Special/Unknown/Gap (X)
}



def calculate_conservation_scores(
    msa_sequences: Sequence[str],
    target_sequence: str,
    total_hits: int,
    wt_scores: Iterable[str],
    type_scores: Iterable[str],
) -> Tuple[List[Dict[str, float]], List[Dict[str, float]], str, str]:
    """
    Calculate multiple wild-type–based and type-based conservation scores
    per residue.

    Parameters
    ----------
    msa_sequences : sequence of str
        Aligned sequences (all same length as target).
    target_sequence : str
        Target (wild-type) sequence.
    total_hits : int
        Number of sequences in the MSA.
    wt_scores : iterable of str
        Names of wild-type metrics to compute.
    type_scores : iterable of str
        Names of type-based metrics to compute.

    Returns
    -------
    wt_metrics : list of dict
        One dict per residue, containing wild-type metrics.
    type_metrics : list of dict
        One dict per residue, containing type-based metrics.
    consensus_seq_res : str
        Consensus residue sequence.
    consensus_seq_type : str
        Consensus type sequence.
    """
    from Bio.Align import substitution_matrices
    BLOSUM62 = substitution_matrices.load("BLOSUM62")

    if not msa_sequences:
        return [], [], "", ""

    wt_scores = set(wt_scores)
    type_scores = set(type_scores)

    seq_len = len(target_sequence)
    wt_metrics: List[Dict[str, float]] = []
    type_metrics: List[Dict[str, float]] = []
    consensus_residues: List[str] = []
    consensus_types: List[str] = []

    for i in range(seq_len):
        wt_residue = target_sequence[i]
        column = [seq[i] for seq in msa_sequences]
        counts = Counter(column)

        wt_dict: Dict[str, float] = {}
        type_dict: Dict[str, float] = {}

        # Consensus residue
        if counts:
            cons_residue, cons_count = max(counts.items(), key=lambda x: x[1])
        else:
            cons_residue, cons_count = "-", 0
        cons_res_freq = cons_count / total_hits if total_hits else 0.0
        wt_dict["cons_residue"] = cons_residue
        if "consensus_freq" in wt_scores:
            wt_dict["consensus_freq"] = cons_res_freq
        consensus_residues.append(cons_residue)

        # WT-based scores
        if "identity" in wt_scores:
            wt_dict["identity"] = counts.get(wt_residue, 0) / total_hits if total_hits else 0.0

        if "entropy" in wt_scores or "norm_entropy" in wt_scores:
            probs = []
            for aa, c in counts.items():
                if aa == "-":
                    continue
                p = c / total_hits
                if p > 0:
                    probs.append(p)
            H = -sum(p * math.log(p, 2) for p in probs) if probs else 0.0
            if "entropy" in wt_scores:
                wt_dict["entropy"] = H
            if "norm_entropy" in wt_scores:
                H_max = math.log(len(probs), 2) if probs else 1.0
                wt_dict["norm_entropy"] = 1.0 - (H / H_max if H_max > 0 else 0.0)

        if "gap_fraction" in wt_scores:
            wt_dict["gap_fraction"] = counts.get("-", 0) / total_hits if total_hits else 0.0

        if "blosum_wt" in wt_scores:
            score_sum = 0.0
            for aa, c in counts.items():
                if aa == "-" or aa == wt_residue:
                    continue
                try:
                    s = BLOSUM62[wt_residue, aa]
                except (KeyError, IndexError):
                    try:
                        s = BLOSUM62[aa, wt_residue]
                    except (KeyError, IndexError):
                        s = 0
                score_sum += s * c
            wt_dict["blosum_wt"] = score_sum / total_hits if total_hits else 0.0

        if "blosum_pairwise" in wt_scores:
            aas = [(aa, c) for aa, c in counts.items() if aa != "-"]
            tot_pairs = 0
            score_sum = 0.0
            for idx_a, (aa1, c1) in enumerate(aas):
                for aa2, c2 in aas[idx_a:]:
                    if aa1 == aa2:
                        num_pairs = c1 * (c1 - 1) // 2 if idx_a == 0 else c1 * c2
                    else:
                        num_pairs = c1 * c2
                    try:
                        s = BLOSUM62[aa1, aa2]
                    except (KeyError, IndexError):
                        try:
                            s = BLOSUM62[aa2, aa1]
                        except (KeyError, IndexError):
                            s = 0
                    score_sum += s * num_pairs
                    tot_pairs += num_pairs
            wt_dict["blosum_pairwise"] = score_sum / tot_pairs if tot_pairs else 0.0

        wt_metrics.append(wt_dict)

        # Type-based scores
        wt_type = AA_TYPES.get(wt_residue.upper(), "X")

        if "type_freq" in type_scores:
            type_count = 0
            for aa, c in counts.items():
                aa_type = AA_TYPES.get(aa.upper(), "X")
                if aa_type == wt_type:
                    type_count += c
            type_dict["type_freq"] = type_count / total_hits if total_hits else 0.0

        if (
            "type_entropy" in type_scores
            or "type_consensus_freq" in type_scores
        ):
            class_counts: Dict[str, int] = {}
            for aa, c in counts.items():
                aa_type = AA_TYPES.get(aa.upper(), "X")
                class_counts[aa_type] = class_counts.get(aa_type, 0) + c

            if class_counts:
                cons_type, cons_type_count = max(class_counts.items(), key=lambda x: x[1])
            else:
                cons_type, cons_type_count = "X", 0
            cons_type_freq = cons_type_count / total_hits if total_hits else 0.0
            type_dict["cons_type"] = cons_type
            if "type_consensus_freq" in type_scores:
                type_dict["type_consensus_freq"] = cons_type_freq

            if "type_entropy" in type_scores:
                probs_t = []
                for _, c in class_counts.items():
                    p = c / total_hits
                    if p > 0:
                        probs_t.append(p)
                Ht = -sum(p * math.log(p, 2) for p in probs_t) if probs_t else 0.0
                type_dict["type_entropy"] = Ht

            consensus_types.append(cons_type)
        else:
            consensus_types.append("X")

        type_metrics.append(type_dict)

    return wt_metrics, type_metrics, "".join(consensus_residues), "".join(consensus_types)


def compute_old_msa_metrics(msa_file):
    # Parse the MSA file
    reformated_msa_file = msa_file.with_suffix('.reformat.fasta')
    if not reformated_msa_file.exists():
        subprocess.run(f"reformat.pl a3m fas -i {msa_file} -o {reformated_msa_file} -r", shell=True, check=True)
    msa = parseMSA(str(reformated_msa_file), format='FASTA')

    # Calculate the entropy and occupancy
    entropy = calcShannonEntropy(msa)
    occupancy = calcMSAOccupancy(msa)
    meff_25 = calcMeff(msa, seqid=0.25)
    meff_50 = calcMeff(msa, seqid=0.50)
    meff_80 = calcMeff(msa, seqid=0.80)
    return entropy, occupancy, meff_25, meff_50, meff_80

def compute_msa_metrics(msa_file):
    # Parse the MSA file
    reformated_msa_file = msa_file.with_suffix('.reformat.fasta')
    if not reformated_msa_file.exists():
        subprocess.run(f"reformat.pl a3m fas -i {msa_file} -o {reformated_msa_file} -r", shell=True, check=True)
    msa = parseMSA(str(reformated_msa_file), format='FASTA')
    #print(f"Parsed MSA from {msa_file}: {len(msa)} sequences, length" )
    #print(msa[1])
    

    # Calculate the entropy and occupancy
    entropy = calcShannonEntropy(msa)
    occupancy = calcMSAOccupancy(msa)
    meff_25 = calcMeff(msa, seqid=0.25)
    meff_50 = calcMeff(msa, seqid=0.50)
    meff_80 = calcMeff(msa, seqid=0.80)
    msa = [str(seq) for seq in msa]
    wt_metrics, wt_type_metrics, _, _ = calculate_conservation_scores(msa, str(msa[0]), len(msa), sequence_metrics, type_metrics)
    pred_metrics, pred_type_metrics, _, _ = calculate_conservation_scores(msa, str(msa[1]), len(msa), sequence_metrics, type_metrics)
    out_wt_metrics = {}
    out_wt_type_metrics = {}
    out_pred_metrics = {}
    out_pred_type_metrics = {}
    for metric in sequence_metrics:
        out_wt_metrics[metric] = sum(res[metric] for res in wt_metrics)/len(wt_metrics) if wt_metrics else 0.0
        out_pred_metrics[metric] = sum(res[metric] for res in pred_metrics)/len(pred_metrics) if pred_metrics else 0.0
    for metric in type_metrics:
        out_wt_type_metrics[metric] = sum(res[metric] for res in wt_type_metrics)/len(wt_type_metrics) if wt_type_metrics else 0.0
        out_pred_type_metrics[metric] = sum(res[metric] for res in pred_type_metrics)/len(pred_type_metrics) if pred_type_metrics else 0.0
    return out_wt_metrics, out_wt_type_metrics, out_pred_metrics, out_pred_type_metrics

def main():
    vanilla_entropies, vanilla_occupancies, vanilla_meffs_25, vanilla_meffs_50, vanilla_meffs_80 = [], [], [], [], []
    synthetic_entropies, synthetic_occupancies, synthetic_meffs_25, synthetic_meffs_50, synthetic_meffs_80 = [], [], [], [], []
    spearman_correlations = []
    conservation_scores = {}
    for pair_path in Path("/data/jgut/msa-tests/zenodo/porter_dataset").glob("[0-9]*"):
        #/data/jgut/msa-tests/zenodo/porter_dataset/1ceeB2k42A/1ceeB_conf.a3m
        for synthetic_msa_path in pair_path.glob("*conf.a3m"):
            #print(f"Processing {synthetic_msa_path}")
            protein_name = synthetic_msa_path.name[:-len("_conf.a3m")]
            #/data/jgut/msa-tests/zenodo/porter_dataset/1ceeB2k42A/1ceeB_alphafold/1inp.a3m
            vanilla_msa_path = pair_path / f"{protein_name}_alphafold" /"1inp.a3m"
            print(f"Vanilla MSA path: {vanilla_msa_path}")
            vanilla_wt_metrics, vanilla_wt_type_metrics, vanilla_pred_metrics, vanilla_pred_type_metrics = compute_msa_metrics(vanilla_msa_path)
            conservation_scores["count"] = conservation_scores.get("count", 0) + 1
            for metric_name in vanilla_wt_metrics.keys():
                conservation_scores[f"vanilla_wt_{metric_name}"] = conservation_scores.get(f"vanilla_wt_{metric_name}", 0) + vanilla_wt_metrics[metric_name]
            for metric_name in vanilla_wt_type_metrics.keys():
                conservation_scores[f"vanilla_wt_type_{metric_name}"] = conservation_scores.get(f"vanilla_wt_type_{metric_name}", 0) + vanilla_wt_type_metrics[metric_name]
            for metric_name in vanilla_pred_metrics.keys():
                conservation_scores[f"vanilla_pred_{metric_name}"] = conservation_scores.get(f"vanilla_pred_{metric_name}", 0) + vanilla_pred_metrics[metric_name]
            for metric_name in vanilla_pred_type_metrics.keys():
                conservation_scores[f"vanilla_pred_type_{metric_name}"] = conservation_scores.get(f"vanilla_pred_type_{metric_name}", 0) + vanilla_pred_type_metrics[metric_name]
            mpnn_wt_metrics, mpnn_wt_type_metrics, mpnn_pred_metrics, mpnn_pred_type_metrics = compute_msa_metrics(synthetic_msa_path)
            for metric_name in mpnn_wt_metrics.keys():
                conservation_scores[f"mpnn_wt_{metric_name}"] = conservation_scores.get(f"mpnn_wt_{metric_name}", 0) + mpnn_wt_metrics[metric_name]
            for metric_name in mpnn_wt_type_metrics.keys():
                conservation_scores[f"mpnn_wt_type_{metric_name}"] = conservation_scores.get(f"mpnn_wt_type_{metric_name}", 0) + mpnn_wt_type_metrics[metric_name]
            for metric_name in mpnn_pred_metrics.keys():
                conservation_scores[f"mpnn_pred_{metric_name}"] = conservation_scores.get(f"mpnn_pred_{metric_name}", 0) + mpnn_pred_metrics[metric_name]
            for metric_name in mpnn_pred_type_metrics.keys():
                conservation_scores[f"mpnn_pred_type_{metric_name}"] = conservation_scores.get(f"mpnn_pred_type_{metric_name}", 0) + mpnn_pred_type_metrics[metric_name]


            vanilla_entropy, vanilla_occupancy, vanilla_meff_25, vanilla_meff_50, vanilla_meff_80 = compute_old_msa_metrics(vanilla_msa_path)
            synthetic_entropy, synthetic_occupancy, synthetic_meff_25, synthetic_meff_50, synthetic_meff_80 = compute_old_msa_metrics(synthetic_msa_path)
            if len(vanilla_entropy) != len(synthetic_entropy):
                print(f"Warning: Length mismatch for {protein_name}. Skipping.")
                continue
            vanilla_entropies.append(np.mean(vanilla_entropy))
            vanilla_occupancies.append(np.mean(vanilla_occupancy))
            vanilla_meffs_25.append(np.mean(vanilla_meff_25))
            vanilla_meffs_50.append(np.mean(vanilla_meff_50))
            vanilla_meffs_80.append(np.mean(vanilla_meff_80))

            synthetic_entropies.append(np.mean(synthetic_entropy))
            synthetic_occupancies.append(np.mean(synthetic_occupancy))
            synthetic_meffs_25.append(np.mean(synthetic_meff_25))
            synthetic_meffs_50.append(np.mean(synthetic_meff_50))
            synthetic_meffs_80.append(np.mean(synthetic_meff_80))
            # Compute Spearman correlation
            spearman_corr, _ = spearmanr(vanilla_entropy, synthetic_entropy)
            spearman_correlations.append(spearman_corr)
    print("Vanilla Entropies:", np.mean(vanilla_entropies))
    print("Synthetic Entropies:", np.mean(synthetic_entropies))
    print("Vanilla Occupancies:", np.mean(vanilla_occupancies))
    print("Synthetic Occupancies:", np.mean(synthetic_occupancies))
    print("Spearman Correlations:", np.mean(spearman_correlations))
    print("Vanilla Meff 25:", np.mean(vanilla_meffs_25))
    print("Synthetic Meff 25:", np.mean(synthetic_meffs_25))
    print("Vanilla Meff 50:", np.mean(vanilla_meffs_50))
    print("Synthetic Meff 50:", np.mean(synthetic_meffs_50))
    print("Vanilla Meff 80:", np.mean(vanilla_meffs_80))
    print("Synthetic Meff 80:", np.mean(synthetic_meffs_80))
    print("Spearman Correlations:", np.tanh(np.mean([np.arctanh(spearman_correlation) for spearman_correlation in spearman_correlations if not np.isnan(spearman_correlation)])))
    print("Conservation Scores:")
    print(f"Count: {conservation_scores['count']}")
    for key, value in conservation_scores.items():
        print(f"{key}: {value / conservation_scores['count']}")


if __name__ == "__main__":
    main()