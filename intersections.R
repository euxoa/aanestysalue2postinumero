library(dplyr)
library(ggplot2)
library(rstan)

# Correlate zip code areas and voting districts by modeling their intersections. 
#
# Suppose Z_i is a variable defined for zip codes, and V_j is a variable relaed to 
# voting districts. 
# 
# Let the intersections be X_ij. For each X_ij, we have Z_i (common with other X_ij), 
# and V_j (common to other X_ij). Suppose there is measurement error in both, and
# let the underlying things be theta_i and psi_j. Now if theta and psi correlate,
# then Z leaks information from X and vice versa. 
#
# The dependency structure is:
# X - theta - psi - Z
#
# Further, the correlation of theta and psi is higher than the observable
# correlation of X and Z. 
#
# A typical setup in the voting context is that X and Z are counts, contrasted to a total, 
# so that at its simplest theta ~ Beta(X, X_total), and psi ~ Beta(Z, Z_total), with
# 0<theta<1 and 0<beta<1. We may model the correlation of theta and beta as linear in the logit space. 
#
# In the case of areas with counts or proportions described above, the counts are not 
# for the intersection ij but for i and j. The beta posterior is then for the whole areas. 
# Intersections are sub-areas with their weights (proportion of total) at least approximately known. 
# Posterior for an intersection is more uncertain than for the whole:
#
# X - theta_i - theta_ij - psi_ij - psi_j - Z
# It is unclear whether theta_ij and psi_ij or their variance around theta_i and psi_j 
# are identifiable at all, so maybe skip that.
# The alternative is to assume theta_i and theta_ij etc. identical,
# which would then underestimate the real correlation between psi and theta.
#
# So:
# - [logit theta_ij, logit psi_ij] ~ multinormal(mean_ij, C)
# - theta_i = \sum_j wx_ij theta_ij (hard constraint)
# - psi_j = \sum_i wz_ij psi_ij
# - Z, X ~ Binom
 
pres_vaal <- readRDS("presidentinvaalien.aanet.rds")

voting <- pres_vaal %>% group_by(kunta, alue) %>% summarise(tot_aania=sum(aania)) %>% ungroup %>%
  left_join(pres_vaal %>% filter(ehdokas=="Haavisto")) %>% 
  filter(tot_aania>0) %>%
  mutate(aanestysalue=as.factor(paste(kunta, alue)))

aa_levels <- levels(voting$aanestysalue)

demography <- readRDS("paavo.2018.rds") %>% select(postinumero, he_vakiy, ko_perus) %>%
  filter(!is.na(he_vakiy) & !is.na(ko_perus)) %>%
  mutate(postinumero = as.factor(postinumero))

pnro_levels <- levels(demography$postinumero)

intersections <- readRDS("map.pono.aanestysalue.rds") %>% 
  select(postinumero, aanestysalue.nro, kunta, n=rakennukset.aanestysalue.pono) %>% 
  mutate(postinumero=factor(postinumero, pnro_levels),
         aanestysalue=factor(paste(kunta, aanestysalue.nro), aa_levels), 
         i = as.integer(postinumero),
         j = as.integer(aanestysalue)) %>%
  filter(!is.na(postinumero) & !is.na(aanestysalue)) %>%
  group_by(postinumero) %>% mutate(wi = n/sum(n)) %>% ungroup() %>%
  group_by(aanestysalue.nro, kunta) %>% mutate(wj = n/sum(n)) %>% ungroup() 

intersections %>% group_by(i) %>% summarise(swi=sum(wi)) %>% filter(swi<.99999 | swi>1)
intersections %>% group_by(j) %>% summarise(swj=sum(wj)) %>% filter(swj<.99999 | swj>1)

# Use a sparse matrix to represent district relationships
Wa <- Matrix::sparseMatrix(i=intersections$i, j=1:nrow(intersections), x=intersections$wi)
# If the weight error is large, your weights are f*ed. They should sum to one column-wise
weight_error <- sum((as.numeric(Wa %*% rep(1, dim(Wa)[2]))-1)**2); stopifnot(weight_error<1.0e-8)
Wa_parts <- extract_sparse_parts(Wa)

Wb <- Matrix::sparseMatrix(i=intersections$j, j=1:nrow(intersections), x=intersections$wj)
# If the weight error is large, your weights are f*ed. They should sum to one column-wise
weight_error <- sum((as.numeric(Wb %*% rep(1, dim(Wb)[2]))-1)**2); stopifnot(weight_error<1.0e-8)
Wb_parts <- extract_sparse_parts(Wb)


d <- list(Na=nrow(demography), Nb=nrow(voting), N=nrow(intersections),
               parent_a = intersections$i,
               parent_b = intersections$j,
               wa = intersections$wi,
               Wa_w = Wa_parts$w, Wa_v = Wa_parts$v, Wa_u = Wa_parts$u,
               wb = intersections$wj,
               a_tot = demography$he_vakiy,
               a_count = demography$ko_perus,
               b_tot = voting$tot_aania,
               b_count = voting$aania)
m <- stan_model("intersections.stan")
fit <- sampling(m, data = d, chains=1, iter=500)


                             