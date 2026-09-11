import re
import subprocess
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from .md_traj_utils import best_hummer_q, clean_traj_c_alpha

MAIN_PRED_COLOUR = "#006400"
MINOR_PRED_COLOUR = "#90EE90"
MEAN_PRED_COLOUR = "#48A948"
MAIN_SIM_COLOUR = "#3941B3"
MINOR_SIM_COLOUR = "#3996B3"
MEAN_SIM_COLOUR = "#396cb3"
OTHER_COLOUR = "gray"
TO_MICRO = 10000000
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

def download_pdb(pdb_code:str, out_path:Path):
    command = f'pdb_fetch {pdb_code} | pdb_selmodel -1 | pdb_selchain -A | pdb_delhetatm | pdb_delinsertion | pdb_reres -1 | pdb_tidy | grep ^ATOM | grep -E "ALA|ARG|ASN|ASP|CYS|GLU|GLN|GLY|HIS|ILE|LEU|LYS|MET|PHE|PRO|SER|THR|TRP|TYR|VAL|SEC|PYL|HCY" > {out_path}'
    subprocess.run(command, shell=True, check=True)

def make_long_plot(simulation, all_rmsds, ref, decision_boundary, save_path, simulation_sampling, frames_per,
                    save_debug_scatter=False, missing_comparison_gt=True, verbose=False):
    """Shared by the three deshaw_*.py scripts, which disagreed on:
    - save_debug_scatter: whether an extra sim_rmsds-vs-time scatter is
      saved to save_path.parent/"test.pdf" before the main plot.
    - missing_comparison_gt: whether a neighbour counts as "missing" via
      rmsd[1] > decision_boundary (True) or rmsd[1] < decision_boundary
      (False).
    - verbose: whether missing_neighbours is printed before plotting.
    frames_per replaces each file's own module-level FRAMES_PER constant
    (2000 for the ovchinnikov-parameterized files, 400 for the real
    unfolding one)."""
    plt.rcParams.update({'axes.spines.right' : True,})
    sim_qs = best_hummer_q(simulation, ref)
    sim_rmsds = md.rmsd(clean_traj_c_alpha(simulation), clean_traj_c_alpha(ref), 0)
    max_sim_frames = len(simulation)
    mean_rmsds = np.mean(sim_rmsds)
    colour_boundary = 2*np.std(sim_rmsds)
    simulation_time_in_micro = max_sim_frames*simulation_sampling*frames_per/TO_MICRO
    sim_x_values = np.linspace(0, simulation_time_in_micro, max_sim_frames)
    if save_debug_scatter:
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
        if missing_comparison_gt:
            missing_neighbour = len([rmsd[0] for rmsd in curr_rmsd if rmsd[0] <max_sim_frames and rmsd[1]>decision_boundary][:FRAMES_TO_EXCLUDE])
        else:
            missing_neighbour = len([rmsd[0] for rmsd in curr_rmsd if rmsd[0] <max_sim_frames and rmsd[1]<decision_boundary][:FRAMES_TO_EXCLUDE])
        missing_neighbours.append(missing_neighbour)
    if verbose:
        print(missing_neighbours)
    clean_neighbours = np.array(clean_neighbours)
    missing_neighbours = np.array(missing_neighbours)
    ax1.scatter(sim_x_values[sim_main_it], sim_rmsds[sim_main_it], c=MAIN_SIM_COLOUR, s=1)
    ax1.scatter(sim_x_values[sim_minor_it], sim_rmsds[sim_minor_it], c=MINOR_SIM_COLOUR, s=1)
    ax1.set_ylabel("RMSD (nm)", color=MEAN_SIM_COLOUR)
    ax1.tick_params(axis='y', labelcolor=MEAN_SIM_COLOUR)
    ax3 = ax1.twinx()
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
    plt.savefig(save_path, format="pdf", bbox_inches='tight', transparent=True)
    plt.close()

def make_long_plot_q_colour(simulation, all_rmsds, ref, decision_boundary, save_path, pred_traj, simulation_sampling, frames_per,
                             save_debug_scatter=True, missing_comparison_gt=False, gate_by_count=True, verbose=True):
    """Shared by the three deshaw_*.py scripts, which disagreed on:
    - save_debug_scatter, missing_comparison_gt, verbose: see make_long_plot.
    - gate_by_count: whether the axvspan shading gate is
      missing_neighbour < FRAMES_TO_EXCLUDE (True) or a plain truthiness
      check on missing_neighbour (False)."""
    plt.rcParams.update({'axes.spines.right' : False,})
    sim_qs = best_hummer_q(simulation, ref)
    pred_qs = best_hummer_q(pred_traj, ref)
    sim_rmsds = md.rmsd(clean_traj_c_alpha(simulation), clean_traj_c_alpha(ref), 0)
    max_sim_frames = len(simulation)
    simulation_time_in_micro = max_sim_frames*simulation_sampling*frames_per/TO_MICRO
    sim_x_values = np.linspace(0, simulation_time_in_micro, max_sim_frames)
    if save_debug_scatter:
        plt.scatter(sim_rmsds, sim_x_values, s=1)
        plt.savefig(save_path.parent/"test.pdf")
        plt.close()
    pred_x_values = np.linspace(0, simulation_time_in_micro, len(all_rmsds.keys()))
    _, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(12, 6))
    clean_neighbours = []
    missing_neighbours = []
    for curr_rmsd in all_rmsds.values():
        clean_neighbour = np.mean([rmsd[0] for rmsd in curr_rmsd if rmsd[0]<max_sim_frames][:FRAMES_TO_AVERAGE])/max_sim_frames*simulation_time_in_micro
        clean_neighbours.append(clean_neighbour)
        if missing_comparison_gt:
            missing_neighbour = len([rmsd[0] for rmsd in curr_rmsd if rmsd[0] <max_sim_frames and rmsd[1]>decision_boundary][:FRAMES_TO_EXCLUDE])
        else:
            missing_neighbour = len([rmsd[0] for rmsd in curr_rmsd if rmsd[0] <max_sim_frames and rmsd[1]<decision_boundary][:FRAMES_TO_EXCLUDE])
        missing_neighbours.append(missing_neighbour)
    if verbose:
        print(missing_neighbours)
    clean_neighbours = np.array(clean_neighbours)
    missing_neighbours = np.array(missing_neighbours)
    ax1.scatter(sim_x_values, sim_rmsds, c=sim_qs, s=1, cmap="viridis")
    ax1.set_ylabel("RMSD (nm)")
    ax2.set_ylabel("Nearest frames")
    ax2.set_ylim([0,simulation_time_in_micro])
    ax2.scatter(pred_x_values, clean_neighbours, c=pred_qs, s=1, cmap="viridis")
    ax2.set_xlabel("Time (μs)")
    half_diff = (pred_x_values[1]-pred_x_values[0])/2
    for ax in [ax2]:
        for pred_x_value, missing_neighbour in zip(pred_x_values, missing_neighbours):
            gate = (missing_neighbour < FRAMES_TO_EXCLUDE) if gate_by_count else bool(missing_neighbour)
            if gate:
                ax.axvspan(pred_x_value-half_diff, pred_x_value+half_diff, color="gray", alpha=0.05)
    # Labels and legends
    plt.tight_layout()
    plt.savefig(save_path, format="pdf", bbox_inches='tight', transparent=True)
    plt.close()

def make_sns_plot(point_list, plt_title, num_predictions=1, save_path=None, variance=None, kde_thresh=None, kde_levels=None):
    """Shared by the three deshaw_*.py scripts. kde_thresh/kde_levels are
    forwarded to sns.kdeplot only when not None, since one of the three
    call sites omitted them entirely rather than passing seaborn's
    defaults explicitly."""
    plt.rcParams.update({'axes.spines.right' : False,})
    pca = pd.DataFrame(point_list[:-num_predictions-1], columns=["PC1", "PC2"])
    kde_kwargs = {}
    if kde_thresh is not None:
        kde_kwargs["thresh"] = kde_thresh
    if kde_levels is not None:
        kde_kwargs["levels"] = kde_levels
    sns.kdeplot(data=pca, x=f"PC1", y=f"PC2", palette=["gray"],color="gray", legend="Simulation", fill=True, alpha=0.8, **kde_kwargs)
    if variance is None:
        plt.xlabel(f"PC1")
        plt.ylabel(f"PC2")
    else:
        plt.xlabel(f"PC1 [variance {variance[0]*100:.2f}%]")
        plt.ylabel(f"PC2 [variance {variance[1]*100:.2f}%]")
    plt.scatter(point_list[-num_predictions:,0], point_list[-num_predictions:,1], label="Prediction", c=MEAN_PRED_COLOUR, s=1)
    handles = [mpatches.Patch(facecolor="gray", label="Simulation"), mpatches.Patch(facecolor=MEAN_PRED_COLOUR, label="Prediction")]
    plt.legend(handles=handles)
    plt.title(plt_title)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, transparent=True)
    plt.show()
    plt.clf()
