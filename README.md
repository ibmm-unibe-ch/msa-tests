# Code for "Deep learning protein folding models predict alternative protein conformations with informative sequence alignments"
## Abstract
> Deep learning models such as AlphaFold accurately predict protein structures but often fail to capture alternative conformations. Here, we use inverse folding to generate synthetic multiple sequence alignments (MSAs) that explicitly encode structural information. When provided with these synthetic MSAs, protein structure prediction models including AlphaFold 2, AlphaFold 3, and RoseTTAFold2 can recover distinct conformations of fold-switching proteins, even without natural homologs. To test robustness, we performed an adversarial experiment where the MSA of each protein was replaced with an MSA encoding a different fold, while keeping the sequence fixed. The models showed three outcomes: correct folding despite the adversarial MSA, unfolded structures, or adoption of the alternative fold, indicating that the effect of the MSA is sequence-dependent. Long-timescale molecular dynamics simulations revealed that well-folded conformations are reproduced accurately, whereas partially unfolded states are redirected toward alternative folds rarely sampled during MD, suggesting that the learned structural distribution favours stable basins. Finally, we demonstrate that this approach enables recovery of alternative conformations when partial structural information about these states is available, providing a practical route to integrate experimental data into full-length structure prediction. This work underscores the potential of input manipulation strategies to expand the conformational repertoire accessible to protein structure prediction models.
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
