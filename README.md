# Code for "Synthetic sequence alignments as programmable probes of learned conformational landscapes in deep learning protein structure predictors"
## Abstract
> What deep learning protein structure predictors learn about conformational landscapes remains largely opaque. Here we introduce synthetic multiple sequence alignments (MSAs), designed by inverse folding to encode predefined structural constraints, as a programmable intervention for interrogating the internal logic of structure prediction systems.

> Synthetic MSAs systematically bias AlphaFold2, AlphaFold3, and RoseTTAFold2 toward distinct conformational states of fold-switching proteins, including alternative conformations inaccessible through natural sequence information alone. Adversarial experiments pairing query sequences with MSAs encoding competing folds reveal sequence-dependent responses, exposing how alignment-derived and sequence-derived signals are weighted within each system. Probing predictions initialized from molecular dynamics trajectories reveals a systematic bias toward compact, training-distribution-favored conformations. Hybrid alignments combining synthetic and natural MSA segments enable targeted steering toward specific conformational states.

> These results establish synthetic MSAs as a generalizable framework for dissecting learned conformational landscapes in deep learning structure predictors, with direct implications for understanding model behavior and accessing biologically relevant hidden states.
## Installation
### Datasets
The **fold-switching dataset** is taken from [Extant fold-switching proteins are widespread](https://www.pnas.org/doi/abs/10.1073/pnas.1800168115) and saved as *data/porter_data.csv* for the working samples and *data/excluded_porter.csv* for the excluded examples.

The **fast-folding simulations** are taken from the paper [How Fast-Folding Proteins fold](https://www.science.org/doi/10.1126/science.1208351).

The **adversary examples** are selected haphazardly by hand and saved to *data/single_proteins.csv*.

### External programs
To run our scripts, consider downloading these tools:
- [AlphaFold3](https://github.com/google-deepmind/alphafold3)
- [BLASTp](https://blast.ncbi.nlm.nih.gov/Blast.cgi)
- [HH-suite](https://github.com/soedinglab/hh-suite)
- [localcolabfold](https://github.com/YoshitakaMo/localcolabfold)
- [MAXIT-suite](https://sw-tools.rcsb.org/apps/MAXIT/index.html)
- [micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)
- [OpenStructure](https://openstructure.org/install) 
- [ProteinMPNN](https://github.com/dauparas/ProteinMPNN)
- [RosettAFold2](https://github.com/uw-ipd/RoseTTAFold2)
### Python environment
To load our Python libraries run `micromamba env create -f environment.yaml`.
## Repository layout
- `analysis/` — Python entry points that run or analyse the pipelines (run these as `python analysis/<script>.py`).
- `bash_scripts/` — shell pipelines (run these as `bash bash_scripts/<script>.sh`).
- `utils/` — shared Python modules imported by the scripts in `analysis/` (e.g. `utils.md_traj_utils`, `utils.pipeline`, `utils.deshaw_common`, `utils.scoring`), plus small standalone CLI helpers (`utils/get_secstrucs.py` and similar).
- `data/` — input CSVs/JSON and generated output (`data/visualisations/`, `data/filter_results/`).
- `result_notebooks/` — analysis notebooks; run with the notebook's own directory as the working directory (Jupyter's default).
- `FrankenMSA.ipynb` stays at the repository root.

All `analysis/*.py` and `bash_scripts/*.sh` invocations below assume the repository root as the working directory. The `bash_scripts/*.sh` scripts that previously hardcoded absolute tool/data paths (`porter_bash.sh`, `process_snapshot.sh`, `process_single.sh`, `process_single_ablation.sh`, `test_one.sh`, `find_duplicates.sh`) now accept `getopts` flags to override those paths (e.g. `-p` for a parent output path, `-m`/`-n` for the ProteinMPNN script, `-i` for an input CSV/directory); run a script with no flags to keep the previous defaults, or pass `-h`/an unknown flag to see its usage line.
## Run the code
### Fold-switching proteins
To generate the data, adjust the paths and then run `bash_scripts/porter_bash.sh`.

To analyse the created data, run `result_notebooks/porter_scores_all_af3.ipynb` to get a good overview and then run `result_notebooks/other_nmr_models.ipynb` to check on NMR structures.

To analyse the MD simulations, run `analysis/md_proteins.py` and check `result_notebooks/MD_proteins.ipynb`.

To check the amount of sequences found in UniProt, run `bash_scripts/find_duplicates.sh` and check the `data/filter_results` folder.

To analyse the generated MSAs, run `analysis/compute_msa_entropy.py`. 
### Test set proteins
For a single protein with two conformations like for the proteins released after the training set cut-off date (SA1, ASCT2, STP10, ZNT8), run `bash_scripts/test_one.sh`.
### Adversarial tests
To run the adversarial tests, use `bash_scripts/process_single.sh` and then analyse with `result_notebooks/Single_proteins.ipynb`.
### MD simulations
To run the recovery pipeline with the fast-folding MD simulations, run `analysis/deshaw_ovchinnikov.py`, before analysing with `analysis/deshaw_ovchinnikov_analysis.py` and the Q-scores with `result_notebooks/Q-score_MMSeqs.ipynb`

To run the recovery pipeline with high-temperature, unfolding MD simulation, run `analysis/deshaw_unfolding.py`, before analysing with `analysis/deshaw_unfolding_analysis.py`.
### FrankenMSA
Check out the application of `FrankenMSA.ipynb` to test our pipeline and combine inverse folded MSAs with traditional MMseqs2 MSAs.
## Contact
If you have questions, please contact jannik.gut@unibe.ch or thomas.lemmin@unibe.ch.
