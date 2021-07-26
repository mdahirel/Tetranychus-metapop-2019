P_tab_plots %>% 
  group_by(METAPOP_ID,SHUFFLE,LENGTH) %>% 
  select(P_alpha,P_mean_all,P_phi,P_gamma) %>% 
  mean_hdi() %>% 
  ggplot()+
  geom_point(aes(P_phi,P_alpha,col=SHUFFLE),size=2)+
  geom_segment(aes(x=P_phi,xend=P_phi,y=P_alpha.lower,yend=P_alpha.upper,col=SHUFFLE))+
  geom_segment(aes(x=P_phi.lower,xend=P_phi.upper,y=P_alpha,yend=P_alpha,col=SHUFFLE))+
  facet_wrap(~SHUFFLE)


cor_phi_alpha<- P_tab_plots %>% 
  group_by(.iteration,SHUFFLE) %>% 
  summarise(cor=cor(P_phi,P_alpha)) 

cor_phi_alpha %>% 
  ungroup() %>% 
  compare_levels(variable=cor,by=SHUFFLE) %>% 
  mean_hdi()

ggplot(cor_phi_alpha)+
  stat_halfeye(aes(cor,SHUFFLE),
           .width=c(0.01,0.95),
           point_interval=mean_hdi)+
  geom_vline(xintercept=0)+
  scale_x_continuous("correlation between phi and alpha")



### based on wang loreau 2015, we can define a ratio of alpha without disp/alpha with disp

### by analogy, why not alpha(without pheno)/alpha(with pheno)

P_tab_plots %>% 
  group_by(.iteration,SHUFFLE,LENGTH) %>% 
  summarise(alpha=mean(P_alpha),phi=mean(P_phi),gamma=mean(P_gamma),N=mean(P_mean_all)) %>% 
  pivot_wider(names_from=SHUFFLE,values_from=c(phi,alpha,gamma,N)) %>% 
  mutate(ratio_alpha=alpha_randomized/alpha_control,
         ratio_phi=phi_control/phi_randomized,
         ratio_gamma=gamma_randomized/gamma_control,
         ratio_N=N_control/N_randomized) %>% 
  ggplot()+
  stat_halfeye(aes(ratio_gamma,LENGTH),
               .width=c(0.01,0.95),
               point_interval=mean_hdi)+
  geom_vline(xintercept=1)
###yet gamma ratio is 1
### so phenotype segregation has no effect on metapop stab desppite destabilisaing local pops
### same as dispersal in general, cf wang loreau
### so: it increase pop sizes without compromising metapop stability, despite destabilising local pops
### but more research needed

### OK we need additive beta to complete the story (see below)
### remove unevenness since we don't use it much? to make room?

## alpha is temporal variability
## additive beta is spatial variability
## the contirbutions of each to gamma is governed by phi the syncrhony
## so we find that self orga increase alpha bad
## but also increase beta additive good

## if we have something that increases alpha by 10% AND beta additive by 10%
## synchrony (beta multiplicative) stay constant
## and we miss some info if we rely on it only

P_tab_fig2_shuffle <- P_tab_plots %>%
  group_by(.iteration, SHUFFLE) %>%
  summarise(
    mean_alpha = mean(P_alpha),
    mean_phi = mean(P_phi),
    mean_gamma = mean(P_gamma),
    mean_uneven = mean(P_uneven),
    mean_beta2=mean(P_beta2)
  ) %>%
  ungroup()

P_tab_fig2_length <- P_tab_plots %>%
  group_by(.iteration, LENGTH) %>%
  summarise(
    mean_alpha = mean(P_alpha),
    mean_phi = mean(P_phi),
    mean_gamma = mean(P_gamma),
    mean_uneven = mean(P_uneven),
    mean_beta2=mean(P_beta2)
  ) %>%
  ungroup()

P_tab_fig
ggplot(P_tab_fig2_shuffle) +
  stat_halfeye(aes(x = mean_beta2, y = SHUFFLE, fill = SHUFFLE),
               orientation = "horizontal",
               .width = c(0.01, 0.95), slab_alpha = 0.7,
               point_interval = mean_hdi
  ) +
  scale_x_continuous("spatial variability (additive beta)") +
  scale_y_discrete("Randomization?") +
  coord_cartesian()

ggplot(P_tab_fig2_length) +
  # actual posterior
  stat_halfeye(aes(x = mean_beta2, y = LENGTH, fill = LENGTH),
               orientation = "horizontal",
               .width = c(0.01, 0.95), slab_alpha = 0.7,
               point_interval = mean_hdi
  ) +
  scale_x_continuous("spatial variability (additive beta)") +
  scale_y_discrete("")


P_tab_fig2_shuffle %>%
  compare_levels(variable = mean_beta2, by = SHUFFLE, fun = `/`) %>% mean_hdi()

P_tab_fig2_length %>%
  compare_levels(variable = mean_beta2, by = LENGTH, fun = `/`) %>% mean_hdi()
