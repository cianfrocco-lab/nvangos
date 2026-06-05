## ----------------------------------------------------------------------
## Title: prsv.R
## Purpose: Calculates the Per-Residue Structural Variance (PRSV) between
##  two structures. Current version only works for tubulin and
##  requires users to specify the alpha- and beta-tubulin chains for each PDB.
##  Segment IDs can also be supplied if needed.
##
## A way to get around chain or segment naming issues is by using bio3D and renaming chains and
##  segids. For example,
##  
##  pdb$atom[pdb$atom$chain=="J",]$segid <- "A01"
##  pdb$atom[pdb$atom$chain=="J",]$chain <- "A"
##
## P. DeLear | Sept Lab @ University of Michigan | Last Updated: 2026/05/16
## ----------------------------------------------------------------------

# libraries needed to perform analysis and plotting
library(tidyverse)
library(ggplot2)
library(bio3d)
library(ggpubr)

apply_packed_pdb_chain_labels <- function(pdb, pdb.file) {
  atom.lines <- readLines(pdb.file, warn = FALSE)
  atom.lines <- atom.lines[grepl("^(ATOM|HETATM)", atom.lines)]

  if (length(atom.lines) != nrow(pdb$atom)) {
    warning("Could not apply packed PDB chain labels because ATOM/HETATM line count does not match Bio3D atom count.")
    return(pdb)
  }

  packed.labels <- trimws(substr(atom.lines, 21, 22))
  blank.labels <- is.na(packed.labels) | packed.labels == ""
  packed.labels[blank.labels] <- pdb$atom$chain[blank.labels]

  if (!"segid" %in% names(pdb$atom)) {
    pdb$atom$segid <- ""
  }

  pdb$atom$segid <- ifelse(nchar(packed.labels) > 1, packed.labels, pdb$atom$segid)
  pdb
}

select_calpha <- function(pdb, chain, residues, segid = NULL) {
  atoms <- pdb$atom
  keep <- atoms$elety == "CA" & atoms$chain == chain & atoms$resno %in% residues

  if (!is.null(segid) && !is.na(segid) && segid != "" && "segid" %in% names(atoms)) {
    keep <- keep & atoms$segid == segid
  }

  atom.inds <- which(keep)

  if (length(atom.inds) == 0) {
    label <- ifelse(is.null(segid) || is.na(segid) || segid == "", chain, paste0(segid, " (chain ", chain, ")"))
    stop("No CA atoms found for ", label, ". Check the selected chain and segment ID.", call. = FALSE)
  }

  residue.order <- match(residues, atoms$resno[atom.inds])
  missing.residues <- residues[is.na(residue.order)]

  if (length(missing.residues) > 0) {
    label <- ifelse(is.null(segid) || is.na(segid) || segid == "", chain, paste0(segid, " (chain ", chain, ")"))
    stop("Missing expected CA residues for ", label, ": ",
         paste(head(missing.residues, 20), collapse = ", "),
         ifelse(length(missing.residues) > 20, ", ...", ""),
         call. = FALSE)
  }

  ordered.atom.inds <- atom.inds[residue.order]

  list(
    atom = ordered.atom.inds,
    xyz = as.vector(rbind(3 * ordered.atom.inds - 2, 3 * ordered.atom.inds - 1, 3 * ordered.atom.inds))
  )
}

# function to calculate prsv used in paper
prsv <- function(base.pdb, comp.pdb,
                 base.alpha.chain = "A", base.beta.chain = "B",
                 comp.alpha.chain = "A", comp.beta.chain = "B",
                 base.alpha.segid = NULL, base.beta.segid = NULL,
                 comp.alpha.segid = NULL, comp.beta.segid = NULL) {
  alpha.residues <- c(seq(1,37),seq(48,437))
  beta.residues <- c(seq(125),seq(130,426))

  # get atom selection of shared residues in chains
  base.ainds <- select_calpha(base.pdb, base.alpha.chain, alpha.residues, base.alpha.segid)
  comp.ainds <- select_calpha(comp.pdb, comp.alpha.chain, alpha.residues, comp.alpha.segid)
  base.binds <- select_calpha(base.pdb, base.beta.chain, beta.residues, base.beta.segid)
  comp.binds <- select_calpha(comp.pdb, comp.beta.chain, beta.residues, comp.beta.segid)

  # make distance maps for alpha and beta chains
  base.a.dm <- dist.xyz(base.pdb$xyz[base.ainds$xyz])
  comp.a.dm <- dist.xyz(comp.pdb$xyz[comp.ainds$xyz])
  base.b.dm <- dist.xyz(base.pdb$xyz[base.binds$xyz])
  comp.b.dm <- dist.xyz(comp.pdb$xyz[comp.binds$xyz])
  # take the difference between the two (a-a, b-b)
  diff.a.dm <- base.a.dm - comp.a.dm
  diff.b.dm <- base.b.dm - comp.b.dm
  
  # takes the variance of the matrix rows 
  prsv.a <- as.vector(apply(diff.a.dm, 1, var))
  prsv.b <- as.vector(apply(diff.b.dm, 1, var))
  
  a.tib <- tibble(Res = alpha.residues, PRSV = prsv.a) %>%
    mutate(Tubulin = 'Alpha',
           Base_chain = base.alpha.chain,
           Comparison_chain = comp.alpha.chain)
  b.tib <- tibble(Res = beta.residues, PRSV = prsv.b) %>%
    mutate(Tubulin = 'Beta',
           Base_chain = base.beta.chain,
           Comparison_chain = comp.beta.chain)

  return(rbind(a.tib,b.tib))
}

plot_prsv_comparisons <- function(gdp.adj, tax.adj, g2p.adj,
                                  gdp.alpha.chain = "A", gdp.beta.chain = "B",
                                  tax.alpha.chain = "A", tax.beta.chain = "B",
                                  g2p.alpha.chain = "A", g2p.beta.chain = "B",
                                  gdp.alpha.segid = NULL, gdp.beta.segid = NULL,
                                  tax.alpha.segid = NULL, tax.beta.segid = NULL,
                                  g2p.alpha.segid = NULL, g2p.beta.segid = NULL,
                                  y.max = 3) {
  if (is.na(y.max) || y.max <= 0) {
    y.max <- 3
  }

  # Example code to reproduce plots in paper from single protofilament PDB files.
  gdp.tax.prsv <- prsv(gdp.adj, tax.adj,
                       base.alpha.chain = gdp.alpha.chain, base.beta.chain = gdp.beta.chain,
                       comp.alpha.chain = tax.alpha.chain, comp.beta.chain = tax.beta.chain,
                       base.alpha.segid = gdp.alpha.segid, base.beta.segid = gdp.beta.segid,
                       comp.alpha.segid = tax.alpha.segid, comp.beta.segid = tax.beta.segid) %>%
    mutate(Comparison = 'GDP-TAX')
  gdp.g2p.prsv <- prsv(gdp.adj, g2p.adj,
                       base.alpha.chain = gdp.alpha.chain, base.beta.chain = gdp.beta.chain,
                       comp.alpha.chain = g2p.alpha.chain, comp.beta.chain = g2p.beta.chain,
                       base.alpha.segid = gdp.alpha.segid, base.beta.segid = gdp.beta.segid,
                       comp.alpha.segid = g2p.alpha.segid, comp.beta.segid = g2p.beta.segid) %>%
    mutate(Comparison = 'GDP-GMPCPP')
  tax.g2p.prsv <- prsv(tax.adj, g2p.adj,
                       base.alpha.chain = tax.alpha.chain, base.beta.chain = tax.beta.chain,
                       comp.alpha.chain = g2p.alpha.chain, comp.beta.chain = g2p.beta.chain,
                       base.alpha.segid = tax.alpha.segid, base.beta.segid = tax.beta.segid,
                       comp.alpha.segid = g2p.alpha.segid, comp.beta.segid = g2p.beta.segid) %>%
    mutate(Comparison = 'TAX-GMPCPP')

  rbind(gdp.tax.prsv, gdp.g2p.prsv, tax.g2p.prsv) %>%
    # Bind rows to add the gap in alpha where the 36-48 residue chain is missing.
    bind_rows(x = ., y = tibble(Res = as.double(rep(seq(37,47),3)), Tubulin = 'Alpha', PRSV = NA,
                                Base_chain = NA, Comparison_chain = NA,
                                Comparison = c(rep('GDP-GMPCPP',11),rep('GDP-TAX',11),rep('TAX-GMPCPP',11)))) %>%
    ggplot() +
    geom_line(aes(x = Res, y = PRSV, color = Comparison)) +
    geom_area(aes(x = Res, y = PRSV, fill = Comparison, group = Comparison), alpha = 0.3) +
    #geom_hline(yintercept = 1.0, linetype = 'dashed') +
    facet_wrap(factor(Comparison, c('GDP-GMPCPP','GDP-TAX','TAX-GMPCPP'))~Tubulin, ncol = 2) +
    coord_cartesian(clip = 'off') +
    scale_y_continuous(limits = c(0,y.max), breaks = seq(0,y.max,1)) +
    labs(x = 'Residue', y = expression(bold(paste('Per Residue Structural Variance ', (ring(A)^2))))) +
    scale_fill_manual(values = c('#299446','#3E87C1','#8B4D9D')) +
    scale_color_manual(values = c('#299446','#3E87C1','#8B4D9D')) +
    theme_pubr(legend = 'none') +
    theme(strip.background = element_blank(),
          strip.text = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
          axis.title = element_text(size = 12, face = 'bold'),
          axis.text = element_text(size = 12, face = 'bold'))
}

