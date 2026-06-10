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

`create_scaled_vector_bild_and_defattr_files.py`

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
