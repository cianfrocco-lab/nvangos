## ----------------------------------------------------------------------
## Title: prsv.R
## Purpose: Calculates the Per-Residue Structural Variance (PRSV) between
##  two structures. Current version only works for tubulin and 
##  requires structures to have chain name set to "A"/"B" and corresponding
##  segment IDs to be "A01, B01, A02, etc"
##
## A way to get around this issue is by using bio3D and renaming chains and 
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

# function to calculate prsv used in paper
prsv <- function(base.pdb, comp.pdb) {
  # get atom selection of shared residues in chains
  base.ainds <- atom.select(base.pdb, "calpha", chain='A', resno=c(seq(1,37),seq(48,437)), segid = 'A02')
  comp.ainds <- atom.select(comp.pdb, "calpha", chain='A', resno=c(seq(1,37),seq(48,437)), segid = 'A02')
  base.binds <- atom.select(base.pdb, "calpha", chain='B', resno=c(seq(125),seq(130,426)), segid = 'B02')
  comp.binds <- atom.select(comp.pdb, "calpha", chain='B', resno=c(seq(125),seq(130,426)), segid = 'B02')
  
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
  
  a.tib <- tibble(Res = c(seq(1,37),seq(48,437)), PRSV = prsv.a) %>% mutate(Chain = 'A')
  b.tib <- tibble(Res = c(seq(125),seq(130,426)), PRSV = prsv.b) %>% mutate(Chain = 'B')
  
  return(rbind(a.tib,b.tib))
}

# example code to reproduce plots in paper from single protofilament PDB files
gdp.tax.prsv <- prsv(gdp.adj, tax.adj) %>% mutate(Comparison = 'GDP-TAX')
gdp.g2p.prsv <- prsv(gdp.adj, g2p.adj) %>% mutate(Comparison = 'GDP-GMPCPP')
tax.g2p.prsv <- prsv(tax.adj, g2p.adj) %>% mutate(Comparison = 'TAX-GMPCPP')


# plotting code 
rbind(gdp.tax.prsv, gdp.g2p.prsv, tax.g2p.prsv) %>%
  # bind rows is to add the gap in alpha where the 36-48 residue chain is missing
  bind_rows(x = ., y = tibble(Res = as.double(rep(seq(37,47),3)), Chain = 'A', PRSV = NA, 
                              Comparison = c(rep('GDP-GMPCPP',11),rep('GDP-TAX',11),rep('TAX-GMPCPP',11)))) %>% 
  ggplot() + 
  geom_line(aes(x = Res, y = PRSV, color = Comparison)) +
  geom_area(aes(x = Res, y = PRSV, fill = Comparison, group = Comparison), alpha = 0.3) +
  #geom_hline(yintercept = 1.0, linetype = 'dashed') +
  facet_wrap(factor(Comparison, c('GDP-GMPCPP','GDP-TAX','TAX-GMPCPP'))~Chain, ncol = 2) +
  coord_cartesian(clip = 'off') +
  scale_y_continuous(limits = c(0,3), breaks = seq(0,3,1)) +
  labs(x = 'Residue', y = expression(bold(paste('Per Residue Structural Variance ', (ring(A)^2))))) + 
  scale_fill_manual(values = c('#299446','#3E87C1','#8B4D9D')) +
  scale_color_manual(values = c('#299446','#3E87C1','#8B4D9D')) +
  theme_pubr(legend = 'none') + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
        axis.title = element_text(size = 12, face = 'bold'),
        axis.text = element_text(size = 12, face = 'bold'))







