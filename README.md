# nvangos

`nvangos` is a small collection of cryo-EM processing and structural-analysis utilities from the Cianfrocco lab. The scripts are intended for users who already work with RELION, ChimeraX, PDB/mmCIF models, and microtubule cryo-EM datasets, but want reproducible examples for common processing and analysis tasks.

## Contents

| File | Purpose |
| --- | --- |
| `per-residue-structural-variance.R` | Calculates per-residue structural variance (PRSV) between tubulin structures and plots pairwise comparisons. |
| `create_scaled_vector_bild_and_defattr_files.py` | Creates ChimeraX `.bild`, `.pb`, and `.defattr` files for visualizing per-residue alpha-carbon displacement vectors between two PDB models. |
| `combine_opposite_register_ptcls_for_seam.py` | Inner-joins two RELION particle STAR files from opposite sides of a 2-protofilament seam stack with alternating registers. |
| `curate_micrographs_by_particle_picks.py` | Filters a RELION `micrographs_ctf.star` file to retain micrographs with particle picks listed in an autopick `summary.star`. |
| `enumerate_optics_groups.py` | Assigns RELION optics groups based on SerialEM image-shift patches before final CTF refinement. |
| `MT_pf_sort_3_336_angpix_180_box.star` | Example RELION reference STAR file for protofilament sorting. |
| `pf_register_sort_box400.star` | Example RELION reference STAR file for protofilament-register sorting. |
| `Detailed Microtubule Cryo-EM RELION Processing Guide.pdf` | Detailed microtubule cryo-EM RELION processing guide. |

## Per-Residue Structural Variance

PRSV compares two tubulin structures by:

1. selecting shared alpha-carbon atoms from alpha- and beta-tubulin chains,
2. calculating intrachain distance maps for each structure,
3. subtracting those distance maps, and
4. reporting the row-wise variance for each residue.

The current `per-residue-structural-variance.R` implementation assumes tubulin PDB objects have:

- user-specified alpha- and beta-tubulin chain IDs for each input PDB
- compatible residue numbering across the two input structure files
- optional segment IDs, if your PDB/mmCIF files require segment-level disambiguation

If your structures use different segment IDs, rename them before running PRSV or pass segment IDs to `prsv()`. For example, in R with `bio3d`:

```r
pdb$atom[pdb$atom$chain == "J",]$segid <- "A02"
pdb$atom[pdb$atom$chain == "J",]$chain <- "A"
```

### Google Colab

Use the PRSV Colab notebook to run the R workflow without setting up R locally:

[Open the PRSV notebook in Colab](https://colab.research.google.com/github/cianfrocco-lab/nvangos/blob/main/notebooks/per_residue_structural_variance_colab.ipynb)

The notebook:

- clones this repository,
- installs the required R packages,
- sources `per-residue-structural-variance.R`,
- reads two user-provided PDB or mmCIF structure files,
- provides interactive widgets for selecting alpha- and beta-tubulin chains in each structure,
- calculates PRSV values, and
- writes CSV and PNG outputs.

## Local Usage

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

## Notes for RELION STAR Files

The two `.star` files in this repository are example RELION reference files. They are useful for understanding expected STAR-file formatting for protofilament and register sorting workflows, but they are not direct inputs to the PRSV R script. PRSV operates on aligned tubulin structural models, typically PDB files.

## Citation

If you use these scripts in published work, please cite the associated Cianfrocco lab work and include a link to this repository.
