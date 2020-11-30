

###let's reorder the existing columns a bit, just to see better

tab <- tab %>% 
  select(c(.iteration,METAPOP_ID,LENGTH, SHUFFLE, 
           P_fixef, P_ranef_rpl, P_latent_intercept, P_vcv_time_latent, P_vcv_rpl_latent,
           M_fixef, M_ranef_rpl, M_var_time_latent, M_var_rpl_latent)
           )

### avging by treatment the temporal VCVs before exponentiating: may be better behaviour for displaying the mean N prediction
### NB: avging before or after doesn't change qualitative conclusions

tab2<-tab %>% 
  group_by(LENGTH,SHUFFLE,.iteration) %>% 
  select(P_vcv_time_latent) %>% 
  summarise(P_vcv_time_latent_list=list(P_vcv_time_latent)) %>% 
  mutate(P_vcv_time_latent_mean=map(
    .x=P_vcv_time_latent_list,
    .f=~.x %>% 
            simplify2array() %>% 
            apply(FUN=mean,MARGIN=c(1,2))))

tab<-tab %>% 
  left_join(tab2)

### we can now start to calculate quantities of interest, replicate by replicate

### first, the predicted mean patch sizes:


### then, the metapop-level predictions



### the next step is the estimation of the observed scale temporal VCV to estimate alpha, beta and gamma

### in theory, for that (and the mean prediction above btw), we should use the QGmvparams function
### but runtime and memory use increase massively with number of variables (here 9), become unusable at about 5 variables
### by chance all of our variables are from the same distribution, and an easy one, so we can do things manually:




### quick check that alpha beta gamma make sense:
### normally, the alpha at the metapop level should be the gamma of the patch level

test1<-tab %>% select(SHUFFLE,LENGTH, METAPOP_ID, P_gamma) %>% 
  group_by(SHUFFLE,LENGTH,METAPOP_ID) %>% 
  mean_qi(P_gamma) %>% 
  rename(low_P_gamma=.lower,high_P_gamma=.upper)

tab %>% select(SHUFFLE,LENGTH, METAPOP_ID, M_alpha) %>% 
  group_by(SHUFFLE,LENGTH,METAPOP_ID) %>% 
  mean_qi(M_alpha) %>% 
  rename(low_M_alpha=.lower,high_M_alpha=.upper) %>% 
  left_join(test1) %>% 
  ggplot()+
  geom_segment(aes(x=M_alpha,xend=M_alpha,y=low_P_gamma,yend=high_P_gamma),col="grey")+
  geom_segment(aes(x=low_M_alpha,xend=high_M_alpha,y=P_gamma,yend=P_gamma),col="grey")+
  geom_point(aes(x=M_alpha,y =P_gamma))+
  scale_x_continuous(expression(paste("metapopulation-level ", alpha)))+
  scale_y_continuous(expression(paste("patch-level ", gamma)))+
  geom_abline(intercept=0,slope=1)


### the problem is the non-linear effects caused by ignoring or not random effects
### as said in villemereuil 2018, the approach above is the only one that recover the correct mean
### and variances
### as long as we at least integrate across the WEEK2 random effects we get compatible results
### BECAUSE EACH METAPOP HAS ITS OWN TEMPORAL VARIABILITY MATRIX
### AND TEMPORAL VAR INTERACT NON LINEARLY WITH FIXED EFFECTS 
### WE CAN'T USE FIXED EFFECTS ALONE TO GET ARITHMETIC MEANS
### TO DO? present both fixed effects-only (geometric means?) 
### and temporal random effect inclusive (arithmetic means)??

## "we note that because each landscape has its own temporal variance matrix, and because of the non
## linearities involved in going from the latent to observed scale
## it is not possible to infer differences between treatment from fixed effect alone
## predictions have to be integrated on the observed scale after accounting for temporal variability
## not doing that lead to severe underpredictions of population means"
## we present nonetheless the differences among treatment means ignoring the temporal variability
## as supplement
## the differences (not the absolute values) may to some extent be interpreted as a conterfactual: how trt would differ
## if the temporal variability was constant across all treatments?

#basically see:
fx_a = rnorm(10000,0,0.1); rn_a=rnorm(10000,0,1)
fx_b = rnorm(10000,0,0.1); rn_b=rnorm(10000,0,2)
fx_c = rnorm(10000,0.2,0.1); rn_c=rnorm(10000,0,1)
mean(exp(fx_a+rn_a)); EnvStats::geoMean(exp(fx_a+rn_a))
mean(exp(fx_b+rn_b)); EnvStats::geoMean(exp(fx_b+rn_b))
mean(exp(fx_c+rn_c)); EnvStats::geoMean(exp(fx_c+rn_c))
mean(exp(fx_a)); EnvStats::geoMean(exp(fx_a))
mean(exp(fx_b)); EnvStats::geoMean(exp(fx_b))
mean(exp(fx_c)); EnvStats::geoMean(exp(fx_c))

##when backtransformed, treatments a and b can differ if their random part differ
## even if the fixed part doesn't
## when we allow random effect to differ between treatments in non Gaussian GLMM, 
### we can't based our inferences 100% on fx only, at least not naively
### correct based on arithmetic means

### naive use of model coefs gives geometric means of patches


## some note on how mean metapop can differ from sum of mean(pop) due to synchrony
## let's imagine two patches x y each sampled 1000 times with a lognormal biomass
x=rlnorm(1000,0,0.5);y=rlnorm(1000,0,0.5)
##here is the distribution of the sum if the patches are independent
fitdistrplus::fitdist(((x)+(y)),"lnorm")

## let's use sort() to mimic what happens if they are + correlated
fitdistrplus::fitdist((sort(x)+sort(y)),"lnorm")

## with rev if they are - correlated

fitdistrplus::fitdist((sort(x)+rev(sort(y))),"lnorm")

##the mean decreases when patches become positively correlated




## Figure 2

Figs 2 to 4 are each made of 4 panels, one per key variable (metapopulation size, $\alpha,\beta,\gamma$ variabilities). We create each panel separately then use the `patchwork` package to merge them and to apply styles common to all panels.

Figure 2 contains posteriors for the non-shuffled (control) metapopulations:
  
  ```{r fig2}

tab2 <- tab_results %>% 
  filter(SHUFFLE == "NO")

p2a <- tab2 %>% 
  ggplot()+
  stat_halfeye(aes(x=metapopNmean,y=factor(LENGTH),fill=factor(LENGTH)),
               point_interval = median_hdi,normalize="xy") +
  labs(x = "Metapopulation mean density")+
  scale_y_discrete("Length (cm)")

p2b <- tab2 %>% 
  ggplot()+
  stat_halfeye(data= . %>% filter(alpha<=20), ##artificial adjustment for plotting, remove beyond 99% quantile for readability
               aes(x=alpha,y=factor(LENGTH),fill=factor(LENGTH)),
               point_interval = NULL,normalize="xy") +
  stat_pointinterval(aes(x=alpha,y=factor(LENGTH)),
                     point_interval = median_hdi)+
  scale_x_continuous("Alpha variability")+
  scale_y_discrete("")

p2c <- tab2 %>% 
  ggplot()+
  stat_halfeye(aes(x=beta1,y=factor(LENGTH),fill=factor(LENGTH)),
               point_interval = median_hdi,normalize="xy")+
  scale_x_continuous("Beta variability")+
  scale_y_discrete("Length (cm)")

p2d <- tab2 %>% 
  ggplot()+
  stat_halfeye(data=. %>% filter(gamma<=15), ##artificial adjustment for plotting, remove beyond 99% quantile for readability
               aes(x=gamma,y=factor(LENGTH),fill=factor(LENGTH)),
               point_interval = NULL,normalize="xy") +
  stat_pointinterval(aes(x=gamma,y=factor(LENGTH)),
                     point_interval = median_hdi)+
  scale_x_continuous("Gamma variability") +
  scale_y_discrete("")

(p2a | p2b)/ (p2c | p2d) + plot_annotation(tag_levels = 'A') & 
  scale_fill_manual(values = paletteLENGTH) &
  theme_bw() & 
  theme(legend.position = "none", text=element_text(size = 20))
```

## Figure 3

Figure 3 contains posteriors for the shuffled metapopulations:
  
  ```{r fig3}

tab3 <- tab_results %>% 
  filter(SHUFFLE == "R")

p3a <- tab3 %>% 
  ggplot()+
  stat_halfeye(aes(x=metapopNmean,y=factor(LENGTH),fill=factor(LENGTH)),
               point_interval = median_hdi,normalize="xy") +
  scale_x_continuous("Metapopulation mean density")+
  scale_y_discrete("Length (cm)")

p3b <- tab3 %>% 
  ggplot()+
  stat_halfeye(data=. %>% filter(alpha<=20), ##artificial adjustment for plotting, remove beyond 99% quantile for readability
               aes(x=alpha,y=factor(LENGTH),fill=factor(LENGTH)),
               point_interval = NULL,normalize="xy") +
  stat_pointinterval(aes(x=alpha,y=factor(LENGTH)),
                     point_interval = median_hdi)+
  scale_x_continuous("Alpha variability")+
  scale_y_discrete("")

p3c <- tab3 %>% 
  ggplot()+
  stat_halfeye(aes(x=beta1,y=factor(LENGTH),fill=factor(LENGTH)),
               point_interval = median_hdi,normalize="xy")+
  scale_x_continuous("Beta variability")+
  scale_y_discrete("Length (cm)")

p3d <- tab3 %>% 
  ggplot()+
  stat_halfeye(data= . %>% filter(gamma<=15), ##artificial adjustment for plotting, remove beyond 99% quantile for readability
               aes(x=gamma,y=factor(LENGTH),fill=factor(LENGTH)),
               point_interval = NULL,normalize="xy") +
  stat_pointinterval(aes(x=gamma,y=factor(LENGTH)),
                     point_interval = median_hdi)+
  scale_x_continuous("Gamma variability") +
  scale_y_discrete("")

(p3a | p3b)/ (p3c | p3d) + plot_annotation(tag_levels = 'A') & 
  scale_fill_manual(values = paletteLENGTH) & 
  theme_bw() & 
  theme(legend.position = "none", text=element_text(size = 20))

```

## Figure 4

Figure 4 contains posteriors of the comparison between shuffled and non shuffled treatments. There are two options here, either use the difference between treatments (chunk `fig4`), or the ratio (`4bis`) between treatments. The ratio may be more appropriate since we are dealing with multiplicative processes:
  
  ### Figure 4: additive comparison
  
  ```{r fig4}

p4a <- tab_results %>% 
  group_by(LENGTH) %>% 
  compare_levels(variable = "metapopNmean", by = SHUFFLE) %>%
  ggplot() +
  stat_halfeye(aes(y = factor(LENGTH), x = metapopNmean, fill=factor(LENGTH)),
               point_interval = median_hdi,normalize="xy") +
  geom_vline(xintercept = 0) + 
  scale_x_continuous("Metapopulation mean density") +
  scale_y_discrete("Length (cm)")

p4b <- tab_results %>% 
  group_by(LENGTH) %>% 
  compare_levels(variable = "alpha", by = SHUFFLE) %>%
  ggplot() +
  stat_halfeye(data=. %>% filter(abs(alpha)<=15), ##artificial adjustment for plotting, remove beyond 99% quantile for readability
               aes(x=alpha,y=factor(LENGTH),fill=factor(LENGTH)),
               point_interval = NULL,normalize="xy") +
  stat_pointinterval(aes(x=alpha,y=factor(LENGTH)),
                     point_interval = median_hdi)+
  geom_vline(xintercept = 0) + 
  scale_x_continuous("Alpha variability") +
  scale_y_discrete("")

p4c <- tab_results %>% 
  group_by(LENGTH) %>% 
  compare_levels(variable = "beta1", by = SHUFFLE) %>%
  ggplot() +
  stat_halfeye(aes(y = factor(LENGTH), x = beta1,fill=factor(LENGTH)),
               point_interval = median_hdi,normalize="xy") +
  geom_vline(xintercept = 0) +
  scale_x_continuous("Beta variability") +
  scale_y_discrete("Length (cm)")

p4d <- tab_results %>% 
  group_by(LENGTH) %>% 
  compare_levels(variable = "gamma", by = SHUFFLE) %>%
  ggplot() +
  stat_halfeye(data=. %>% filter(abs(gamma)<=5), ##artificial adjustment for plotting, remove beyond 99% quantile for readability
               aes(x=gamma,y=factor(LENGTH),fill=factor(LENGTH)),
               point_interval = NULL,normalize="xy") +
  stat_pointinterval(aes(x=gamma,y=factor(LENGTH)),
                     point_interval = median_hdi)+
  geom_vline(xintercept = 0) +
  scale_x_continuous("Gamma variability") +
  scale_y_discrete("")


(p4a | p4b)/ (p4c | p4d) + plot_annotation(title = "Comparisons, shuffled - control", tag_levels = 'A') & 
  theme_bw() & 
  scale_fill_manual(values = paletteLENGTH) & 
  theme(legend.position = "none", text=element_text(size = 20))
```

### Figure 4bis: multiplicative comparison

```{r fig4bis}
p4a <- tab_results %>% 
  group_by(LENGTH) %>% 
  compare_levels(variable = "metapopNmean", by = SHUFFLE, fun= `/`) %>%
  ggplot() +
  stat_halfeye(aes(y = factor(LENGTH), x = metapopNmean, fill=factor(LENGTH)),
               point_interval = median_hdi,normalize="xy") +
  geom_vline(xintercept = 1) + 
  scale_x_continuous("Metapopulation mean density") +
  scale_y_discrete("Length (cm)")

p4b <- tab_results %>% 
  group_by(LENGTH) %>% 
  compare_levels(variable = "alpha", by = SHUFFLE, fun= `/`) %>%
  ggplot() +
  stat_halfeye(data=. %>% filter(alpha<=10), ##artificial adjustment for plotting, remove beyond 99% quantile for readability
               aes(x=alpha,y=factor(LENGTH),fill=factor(LENGTH)),
               point_interval = NULL,normalize="xy") +
  stat_pointinterval(aes(x=alpha,y=factor(LENGTH)),
                     point_interval = median_hdi)+
  geom_vline(xintercept = 1) + 
  scale_x_continuous("Alpha variability") +
  scale_y_discrete("")

p4c <- tab_results %>% 
  group_by(LENGTH) %>% 
  compare_levels(variable = "beta1", by = SHUFFLE, fun= `/`) %>%
  ggplot() +
  stat_halfeye(aes(y = factor(LENGTH), x = beta1,fill=factor(LENGTH)),
               point_interval = median_hdi,normalize="xy") +
  geom_vline(xintercept = 1) +
  scale_x_continuous("Beta variability") +
  scale_y_discrete("Length (cm)")

p4d <- tab_results %>% 
  group_by(LENGTH) %>% 
  compare_levels(variable = "gamma", by = SHUFFLE, fun= `/`) %>%
  ggplot() +
  stat_halfeye(data=. %>% filter(gamma<=5), ##artificial adjustment for plotting, remove beyond 99% quantile for readability
               aes(x=gamma,y=factor(LENGTH),fill=factor(LENGTH)),
               point_interval = NULL,normalize="xy") +
  stat_pointinterval(aes(x=gamma,y=factor(LENGTH)),
                     point_interval = median_hdi)+
  geom_vline(xintercept = 1) +
  scale_x_continuous("Gamma variability") +
  scale_y_discrete("")


(p4a | p4b)/ (p4c | p4d) + plot_annotation(title = "Comparisons, ratio shuffled / control", tag_levels = 'A') & 
  theme_bw() & 
  scale_fill_manual(values = paletteLENGTH) & 
  theme(legend.position = "none", text=element_text(size = 20))
```

## Figure 5

Figure 5 shows what happens at the local patch scale, depending on the local level of patch connectedness:
  
  ```{r fig5}
tab_results %>% 
  select(Nmean_corners,Nmean_center,Nmean_sides, LENGTH, SHUFFLE,.iteration) %>%
  pivot_longer(cols=c(Nmean_corners,Nmean_center,Nmean_sides),
               values_to="N_local") %>% 
  mutate(name = fct_recode(factor(name), corner="Nmean_corners",edge="Nmean_sides",center="Nmean_center")) %>% 
  mutate(`local connectedness` = fct_relevel(name,"corner",after=Inf)) %>% 
  mutate(SHUFFLE = fct_recode(factor(SHUFFLE),control="NO",randomized="R")) %>% 
  ggplot()+
  stat_halfeye(aes(x=N_local,y=factor(LENGTH),fill=`local connectedness`, alpha=0.5))+
  scale_x_continuous("mean local population density")+
  scale_fill_manual(values=paletteLOCAL)+
  scale_y_discrete("Length (cm)")+
  scale_alpha(guide="none")+
  facet_grid(rows=vars(SHUFFLE))+
  theme_bw() & 
  theme(legend.position = "bottom", text=element_text(size = 20))

```

# More on inferences

The function `compare_levels()` used to make Figure 4 can be used to obtain actual numbers and intervals for many comparisons. Only one example displayed here for simplicity, see the help of both functions for more details, and explore:
  
  ```{r compare}
## an example
## calculating differences between treatment using compare_levels
compare_levels(tab_results, variable = "metapopNmean", by = TREATMENT) %>%
  mean_hdi()


```
