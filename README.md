# nvangos

`nvangos` is a collection of cryo-EM processing and structural-analysis utilities from the Cianfrocco lab for microtubule structure determined, related to [Vangos et al. 2026](https://www.biorxiv.org/content/10.64898/2026.05.17.725690v1). The scripts are intended for users who already work with microtubule cryo-EM datasets and structures that have been processed in RELION and ChimeraX. 

## Contents

| File | Purpose |
| --- | --- |
| `Detailed Microtubule Cryo-EM RELION Processing Guide.pdf` | Step-by-step guide to determine one and two protofilament microtubule structures using RELION. |
| `per-residue-structural-variance.R` | Calculates per-residue structural variance (PRSV) between two structures of the same protein and plots pairwise comparisons. Can be run from Google Colab Notebook or locally|
| `notebook/per_residue_structural_variance_colab.ipynb` | Python notebook that can run per-residue-structural-variance.R on Goolge Colab |
| `create_scaled_vector_bild_and_defattr_files.py` | Creates ChimeraX `.bild`, `.pb`, and `.defattr` files for visualizing per-residue alpha-carbon displacement vectors between two PDB models of the same protein/biomolecule. |
| `combine_opposite_register_ptcls_for_seam.py` | Inner-joins two RELION particle STAR files from opposite sides of a 2-protofilament seam stack with alternating registers. |
| `curate_micrographs_by_particle_picks.py` | Filters a RELION `micrographs_ctf.star` file to retain micrographs with particle picks listed in an autopick `summary.star`. |
| `enumerate_optics_groups.py` | Assigns RELION optics groups based on SerialEM image-shift patches before final CTF refinement. |
| `MT_pf_sort_3_336_angpix_180_box.star` | Example RELION reference STAR file for protofilament sorting. |
| `pf_register_sort_box400.star` | Example RELION reference STAR file for protofilament-register sorting. |

## Per-Residue Structural Variance

Per-residue structural variance (PRSV) calculates the variance of a given alpha carbon between two different structures. It does this by comparing two structures by:

1. selecting shared alpha-carbon atoms from alpha- and beta-tubulin chains,
2. calculating intrachain distance maps for each structure,
3. subtracting those distance maps, and
4. reporting the row-wise variance for each residue.

The current `per-residue-structural-variance.R` implementation assumes PDB objects have:

- user-specified alpha- and beta-tubulin chain IDs for each input PDB
- compatible residue numbering across the two input structure files
- optional segment IDs, if your PDB/mmCIF files require segment-level disambiguation

If your structures use different segment IDs, rename them before running PRSV, indicate the chain IDs using the Google Colab notebook, or pass segment IDs to `prsv()`. For example, in R with `bio3d`:

```r
pdb$atom[pdb$atom$chain == "J",]$segid <- "A02"
pdb$atom[pdb$atom$chain == "J",]$chain <- "A"
```

Some PHENIX/ChimeraX-written PDB files encode labels such as `A1`, `A2`, `B1`, and `B2` across the residue-name/chain fields rather than in the single-character PDB chain column. The Colab notebook detects this packed format and passes the reconstructed labels to PRSV as segment IDs.

### Google Colab

Use the PRSV Colab notebook to run the R workflow without setting up R locally:

[Open the PRSV notebook in Colab](https://colab.research.google.com/github/cianfrocco-lab/nvangos/blob/main/notebooks/per_residue_structural_variance_colab.ipynb)

The notebook:

- clones this repository,
- installs the required R packages,
- sources `per-residue-structural-variance.R`,
- reads two user-provided PDB or mmCIF structure files,
- provides interactive widgets for selecting alpha- and beta-tubulin chains in each structure,
- provides a user-defined plot y-axis maximum, defaulting to `3`,
- calculates PRSV values, and
- writes CSV and PNG outputs.

### Local Usage

Clone the repository:

```bash
git clone https://github.com/cianfrocco-lab/nvangos.git
cd nvangos
```

Install R dependencies for PRSV:

```r
install.packages(c("tidyverse", "bio3d", "ggpubr"), repos = "https://cloud.r-project.org")
```

Load the PRSV functions:

```r
source("per-residue-structural-variance.R")
```

Calculate PRSV for two adjusted PDB objects:

```r
base <- bio3d::read.pdb("base_model.pdb")
comparison <- bio3d::read.pdb("comparison_model.pdb")

prsv_values <- prsv(
  base,
  comparison,
  base.alpha.chain = "A",
  base.beta.chain = "B",
  comp.alpha.chain = "A",
  comp.beta.chain = "B"
)
write.csv(prsv_values, "prsv_values.csv", row.names = FALSE)
```

## Visualizing structural changes using ChimeraX

To display structural changes between two structural states, you can calculate vectors and attribute files with this script. The output files can be opened in ChimeraX. 

'''$ python create_scaled_vector_bild_and_defattr_files.py --help
usage: create_scaled_vector_bild_and_defattr_files.py [-h] [--chains CHAINS [CHAINS ...]]
                                                      [--file-1-number FILE_1_NUMBER] [--file-2-number FILE_2_NUMBER]
                                                      [--color COLOR] [--radius RADIUS] [--dashes DASHES]
                                                      [--out-ca OUT_CA] [--out-bild-xyz OUT_BILD_XYZ]
                                                      [--out-pb-xyz OUT_PB_XYZ] [--out-ca-xy OUT_CA_XY]
                                                      [--out-bild-xy OUT_BILD_XY] [--out-ca-z OUT_CA_Z]
                                                      [--out-bild-z OUT_BILD_Z]
                                                      file_1 file_2

Generate ChimeraX .defattr, .bild, and .pb files from per-residue CA vectors between two aligned PDB models.

positional arguments:
  file_1                First aligned PDB model.
  file_2                Second aligned PDB model.

optional arguments:
  -h, --help            show this help message and exit
  --chains CHAINS [CHAINS ...]
                        Chain IDs to analyze, in order. If omitted, all matching CA chain IDs present in both models are
                        used.
  --file-1-number FILE_1_NUMBER
                        ChimeraX model number for file 1. Default: 1
  --file-2-number FILE_2_NUMBER
                        ChimeraX model number for file 2. Default: 2
  --color COLOR         Pseudobond color written to the .pb header. Default: blue
  --radius RADIUS       Pseudobond radius written to the .pb header. Default: 0.3
  --dashes DASHES       Pseudobond dashes value written to the .pb header. Default: 1
  --out-ca OUT_CA       Output CA distance .defattr file. Default: ca_distances.defattr beside file_1
  --out-bild-xyz OUT_BILD_XYZ
                        Output XYZ .bild file. Default: colored_vectors_XYZ.bild beside file_1
  --out-pb-xyz OUT_PB_XYZ
                        Output XYZ .pb file. Default: colored_vectors_XYZ.pb beside file_1
  --out-ca-xy OUT_CA_XY
                        Output XY distance .defattr file. Default: ca_distances_XY_only.defattr beside file_1
  --out-bild-xy OUT_BILD_XY
                        Output XY .bild file. Default: colored_vectors_XY_only.bild beside file_1
  --out-ca-z OUT_CA_Z   Output Z distance .defattr file. Default: ca_distances_Z_only.defattr beside file_1
  --out-bild-z OUT_BILD_Z
                        Output Z .bild file. Default: colored_vectors_Z_only.bild beside file_1'''

Example running with two input files (that are aligned to each other) using default options: 
`$ python create_scaled_vector_bild_and_defattr_files.py file1.pdb file2-alignedToFile1.pdb`

This will output the following files:
* ca_distances.defattr - attribute file to color models in ChimeraX according to alpha carbon distances
* colored_vectors_XYZ.bild - vectors colored according to X, Y, & Z to be opened in ChimeraX
* colored_vectors_XYZ.pb - vectors colored according to X, Y, & Z to be opened in ChimeraX
* ca_distances_XY_only.defattr - attribute file coloring models in ChimeraX according to X & Y alpha carbon distances only
* colored_vectors_XY_only.bild - vectors colored according to X & Y displacements only for alpha carbon differences
* ca_distances_Z_only.defattr - attribute file to color models according to Z-only alpha carbon differences
* colored_vectors_Z_only.bild - vector file to color models according to Z-only alpha carbon differences 

## RELION-based cryo-EM processing of microtubule datasets 

We provide the complete step-by-step guidance on how to use RELION to process cryo-EM datasets of microtubules (`Detailed Microtubule Cryo-EM RELION Processing Guide.pdf`). To run this guide, you will need the following scripts to help in data curation and particle STAR file manipulation. 

### Scripts to help data preparation/filtering

`curate_micrographs_by_particle_picks.py` - Filters a RELION `micrographs_ctf.star` file to retain micrographs with particle picks listed in an autopick `summary.star`

`enumerate_optics_groups.py` - Assigns RELION optics groups based on SerialEM image-shift patches before final CTF refinement.

### Identifying seam particles

After symmetry expanding and identifying the alpha/beta registry for each side of the seam, you can run this script to find which particles comprise 'seam' particles. 

`combine_opposite_register_ptcls_for_seam.py` - Inner-joins two RELION particle STAR files from opposite sides of a 2-protofilament seam stack with alternating registers.

## Citation

If you use these scripts in published work, please cite [Vangos et. al bioRxiv 2026.05.17.725690](https://www.biorxiv.org/content/10.64898/2026.05.17.725690v1).
