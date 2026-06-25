# Code for "Synthetic sequence alignments as programmable probes of learned conformational landscapes in deep learning protein structure predictors"
## Abstract
> What deep learning protein structure predictors learn about conformational landscapes remains largely opaque. Here we introduce synthetic multiple sequence alignments (MSAs), designed by inverse folding to encode predefined structural constraints, as a programmable intervention for interrogating the internal logic of structure prediction systems.

> Synthetic MSAs systematically bias AlphaFold2, AlphaFold3, and RoseTTAFold2 toward distinct conformational states of fold-switching proteins, including alternative conformations inaccessible through natural sequence information alone. Adversarial experiments pairing query sequences with MSAs encoding competing folds reveal sequence-dependent responses, exposing how alignment-derived and sequence-derived signals are weighted within each system. Probing predictions initialized from molecular dynamics trajectories reveals a systematic bias toward compact, training-distribution-favored conformations. Hybrid alignments combining synthetic and natural MSA segments enable targeted steering toward specific conformational states.

> These results establish synthetic MSAs as a generalizable framework for dissecting learned conformational landscapes in deep learning structure predictors, with direct implications for understanding model behavior and accessing biologically relevant hidden states.
## Installation
### Datasets
The **fold-switching dataset** is taken from [Extant fold-switching proteins are widespread](https://www.pnas.org/doi/abs/10.1073/pnas.1800168115) and saved as *porter_data.csv* for the working samples and *porter_excluded.csv* for the excluded examples.

The **fast-folding simulations** are taken from the paper [How Fast-Folding Proteins fold](https://www.science.org/doi/10.1126/science.1208351).

The **adversary examples** are selected haphazardly by hand and saved to *single_proteins.csv*.

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
## Run the code
### Fold-switching proteins
To generate the data, adjust the paths and then run `porter_bash.sh`.

To analyse the created data, run `porter_scores_all_af3.ipynb` to get a good overview and then run `other_nmr_models.ipynb` to check on NMR structures.

To analyse the MD simulations, run `md_proteins.py` and check `MD_proteins.ipynb`.

To check the amount of sequences found in UniProt, run `find_duplicates.sh` and check the `filter_results` folder.
### Adversarial tests
To run the adversarial tests, use `process_single.sh` and then analyse with `Single_proteins.ipynb`.
### MD simulations
To run the recovery pipeline with the fast-folding MD simulations, run `deshaw_ovchinnikov.py`, before analysing with `deshaw_ovchinnikov_analysis.py` and the Q-scores with `Q-score_MMSeqs.ipynb`

To run the recovery pipeline with high-temperature, unfolding MD simulation, run deshaw_unfolding.py, before analysing with `deshaw_unfolding_analysis.py`.
### FrankenMSA
Check out the application of `FrankenFold.ipynb` to test our pipeline and combine inverse folded MSAs with traditional MMseqs2 MSAs.
## Contact
If you have questions, please contact jannik.gut@unibe.ch or thomas.lemmin@unibe.ch.
