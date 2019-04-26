library(dplyr)
library(ggplot2)
library(rstan)

aanet_puolueittain <- readRDS(file = "data/EKV2019_ehdokkaat.rds") %>% 
  filter(aluejako == "äänestysalue") %>% 
  rename(aanestysalue.nro = alue) %>%
  select(-puolue) %>% rename(puolue = puolue_lyhenne_alkuperainen) %>%
  group_by(vaalipiiri, kunta, aanestysalue.nro, puolue) %>%
  summarise(aania = sum(aanet_yhteensa))

voting <- aanet_puolueittain %>% 
  group_by(vaalipiiri, kunta, aanestysalue.nro) %>% summarise(tot_aania=sum(aania)) %>% ungroup %>%
  left_join(aanet_puolueittain %>% filter(puolue=="PS")) %>% 
  filter(tot_aania>0) %>%
  mutate(aanestysalue=paste(kunta, aanestysalue.nro))

paavo <- readRDS("map_and_names/paavodata.rds")$data %>% tbl_df()

demography <- paavo %>% filter(vuosi==2019) %>% select(postinumero=pono, he_vakiy, ko_perus) %>%
  filter(!is.na(he_vakiy) & !is.na(ko_perus)) 

intersections <- readRDS("data/aanestysalue2postinumero.rds") %>% 
  select(postinumero, aanestysalue.nro, kunta, 
         n = rakennukset.aanestysalue.pono, 
         wi = w.aanestysalue2pono, 
         wj = w.pono2aanestysalue) %>% 
  mutate(aanestysalue = paste(kunta, aanestysalue.nro),
         aanestysalue = as.factor(aanestysalue),
         postinumero = as.factor(postinumero),
         i = as.numeric(aanestysalue),
         j = as.numeric(postinumero)) 
  
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

margins <- intersections %>% select(aanestysalue, postinumero) %>% sapply(levels)

i.voting <- intersections %>% left_join(voting) %>% select(i, names(voting)) %>% unique %>% arrange(i)
j.demography <- intersections %>% left_join(demography) %>% select(j, names(demography)) %>% unique %>% arrange(j)

# Missing data in the counts is treated as zero. I guess that's kind of ok,
# for voting that is probably the real situation, and for demography binomials
# with zero counts are treated as missing evidence.
na20 <- function(x) ifelse(is.na(x), 0, x)

# But there are districts with missing puolue and zero votes for the party, although
# the total number of votes is >1000!
# Make the evidence zero for these, for in these cases there has been no candidates.
fucked_votes <- with(i.voting, is.na(puolue) | is.na(aania) | is.na(tot_aania))
i.voting$aania[fucked_votes] <- 0
i.voting$tot_aania[fucked_votes] <- 0

intersections$b_groups <- as.integer(factor(intersections %>% left_join(voting) %>% {.$vaalipiiri}, exclude=c()))

d <- list(Na=nrow(j.demography), Nb=nrow(i.voting), N=nrow(intersections),
               Wa_w = Wb_parts$w, Wa_v = Wb_parts$v, Wa_u = Wb_parts$u,
               Wb_w = Wa_parts$w, Wb_v = Wa_parts$v, Wb_u = Wa_parts$u,
               a_tot = j.demography$he_vakiy %>% na20,
               a_count = j.demography$ko_perus %>% na20,
               b_tot = i.voting$tot_aania %>% na20,
               b_count = i.voting$aania %>% na20, 
               Nb_groups = max(intersections$b_groups),
               b_groups = intersections$b_groups)
m <- stan_model("intersections.stan")
fit <- sampling(m, data = d, chains=1, iter=500)
traceplot(fit, c("theta_mu", "theta_sigma", "psi_mu", "psi_sigma", "psi_mu_sigma", "a"))

j.demography %>% mutate(he_vakiy = na20(he_vakiy), ko_perus=na20(ko_perus)) %>% ggplot(aes(x=he_vakiy, y=ko_perus)) + geom_point() + scale_x_log10() + scale_y_log10()
i.voting %>% mutate(tot_aania = na20(tot_aania), aania=na20(aania)) %>% ggplot(aes(x=tot_aania, y=aania, color=vaalipiiri)) + geom_point() + facet_wrap(~vaalipiiri) + scale_x_log10() + scale_y_log10()

# TODO
# - tyhjät vaalipiirit pois
# - vaalipiirikohtainen offset koska ehdokasasettelu

theta <- apply(extract(fit, "theta")[[1]], 2, mean)
psi <- apply(extract(fit, "psi")[[1]], 2, mean)
plot(theta, psi, pch=".")
