

### we first extract all the fixed effects:
FIXEF <- fixef(mod, summary = FALSE) %>% ### all model fixed effects, with one column = one "variable" * treatment combination
  as_tibble() %>%
  mutate(.iteration = 1:dim(.)[1]) %>% ### we put the implicit iteration # into an explicit column
  pivot_longer(-.iteration) %>%        ### we pivot so we can split the "variable*treatment" into patch ID and treatment variables
  mutate(variable=str_extract(name,"P11|P12|P13|P21|P22|P23|P31|P32|P33|METAPOPSUM"),
         LENGTH=as.numeric(str_extract(name,"4|8|16")),
         SHUFFLE=str_extract(name,"NO$|R$")) %>%  # $ : end of string
  mutate(TREATMENT = interaction(LENGTH,SHUFFLE)) %>% 
  select(-name) %>%
  pivot_wider(values_from=value,names_from=variable) %>%  ###we repivot back so that 1 column per patch
  group_by(.iteration, LENGTH, SHUFFLE, TREATMENT) %>%  ### we group and nest by iteration (with treatment variables along for the ride)
  nest(fixef_patch = c(P11:P33)) %>% 
  mutate(fixef_patch = map(.x=fixef_patch, .f = ~.x %>% unlist())) %>%   #convert from tibble to vector
  rename(fixef_metapop = METAPOPSUM)

### we then extract the among-metapop deviations for the patch by patch model:
RANEF_REPLICATE<- ranef(mod, summary = FALSE)$METAPOP_ID %>% 
  array_tree(margin=c(1,2)) %>% ## convert array to list, enlisting by metapop nested in iteration
  tibble() %>% ## we "flatten" the list, iteration become rows of new tibble
  rename(., ranef_rpl_patch=".") %>% 
  mutate(.iteration = 1:dim(.)[1]) %>% ### we put the implicit iteration # into an explicit column
  unnest_longer(ranef_rpl_patch, indices_to = "METAPOP_ID") %>% 
  mutate(LENGTH=as.numeric(str_extract(METAPOP_ID,"^4|^8|^16")),  ## ^: start of string
         SHUFFLE=str_extract(METAPOP_ID,"NO|R")) %>% 
  mutate(TREATMENT = interaction(LENGTH,SHUFFLE))

### we do the same for the metapop sum model:
RANEF_REPLICATE2<- ranef(mod, summary = FALSE)$METAPOP_ID2 %>%
  array_tree(margin=c(1)) %>% ## convert array to list, enlisting by metapop nested in iteration
  tibble() %>% ## we "flatten" the list, iteration become rows of new tibble
  rename(., ranef_rpl_metapop=".") %>% 
  mutate(ranef_rpl_metapop = map(.x = ranef_rpl_metapop,
                                     .f = ~ t(.x) %>% as.data.frame() %>% c())) %>% 
  mutate(.iteration = 1:dim(.)[1]) %>% ### we put the implicit iteration # into an explicit column
  unnest_longer(ranef_rpl_metapop, indices_to = "METAPOP_ID") %>% 
  mutate(LENGTH=as.numeric(str_extract(METAPOP_ID,"^4|^8|^16")),  ## ^: start of string
         SHUFFLE=str_extract(METAPOP_ID,"NO|R")) %>% 
  mutate(TREATMENT = interaction(LENGTH,SHUFFLE))

### now we need to extract the temporal variance-covariance matrices:
GLOBAL_VCV_time <- VarCorr(mod, summary = FALSE)$WEEK2$cov %>% 
      array_tree(margin=1) %>% 
      tibble() %>% 
  rename(., vcv_global_latent=".") %>% 
  ### this is the VCV matrix including all the patches 
  ### so including the zero covs between the patchs from one metapop and the patch from another
  mutate(.iteration = 1:dim(.)[1])

VCV_rpl_patch <- VarCorr(mod, summary = FALSE)$METAPOP_ID$cov %>% 
  array_tree(margin=1) %>% 
  tibble() %>% 
  rename(., vcv_rpl_latent=".") %>% 
  mutate(.iteration = 1:dim(.)[1])



VAR_time_metapop <- VarCorr(mod,summary=FALSE)$WEEK3$sd %>%
  as.data.frame() %>% 
  mutate(.iteration = 1:dim(.)[1]) %>% 
  pivot_longer(-.iteration) %>% 
  mutate(latent_time_var_metapop = value^2) %>% 
  mutate(METAPOP_ID = str_remove(name,"METAPOPSUM_Intercept:METAPOP_ID2")) %>% 
  select(METAPOP_ID,latent_time_var_metapop, .iteration) %>% 
  mutate(LENGTH=as.numeric(str_extract(METAPOP_ID,"^4|^8|^16")),  ## ^: start of string
         SHUFFLE=str_extract(METAPOP_ID,"NO|R")) %>% 
  mutate(TREATMENT = interaction(LENGTH,SHUFFLE)) 

VAR_rpl_metapop <- VarCorr(mod,summary=FALSE)$METAPOP_ID2$sd %>% 
  as.data.frame() %>% 
  mutate(.iteration = 1:dim(.)[1]) %>% 
  mutate(latent_rpl_var_metapop = METAPOPSUM_Intercept^2) %>% 
  select(.iteration,latent_rpl_var_metapop)

### and we merge all these objects into one table:
### this will in addition help us split the global VCV matrix into one VCV matrix per metapop:
tab <- RANEF_REPLICATE %>% 
  left_join(RANEF_REPLICATE2) %>% 
  left_join(VAR_time_metapop) %>%  
  left_join(VAR_rpl_metapop) %>% 
  left_join(FIXEF) %>% 
  left_join(VCV_rpl_patch) %>% 
  inner_join(GLOBAL_VCV_time, by = ".iteration") %>% 
  mutate(vcv_metapop_latent = map2(
    .x = vcv_global_latent, .y = METAPOP_ID,
    .f = ~.x %>% .[str_detect(colnames(.x),.y),str_detect(colnames(.x),.y)]
  )) %>% 
  mutate(latent_intercepts = map2( ### I need to change the variables names
    .x = fixef_patch, .y = ranef_rpl_patch, ##the latent metapop level intercept is the sum of the global mean + the metapop level ranef
    .f= function(.x,.y){.x + .y}
  ))

###let's reorder the existing columns a bit, just to see better

tab <- tab %>% 
  select(c(.iteration:TREATMENT, METAPOP_ID, 
           fixef_patch, ranef_rpl_patch, latent_intercepts, vcv_rpl_latent, vcv_metapop_latent,
           fixef_metapop, ranef_rpl_metapop, latent_rpl_var_metapop, latent_time_var_metapop)
           )

### we can now start to calculate quantities of interest, replicate by replicate

### first, the predicted mean patch sizes:

tab <- tab %>% 
  mutate(total_vcv = map2(.x = vcv_rpl_latent, .y = vcv_metapop_latent,
                          .f = function(.x,.y){.x+.y})) %>% 
  mutate(patch_avgs = map2( ## patch by patch average
    .x = latent_intercepts, .y = total_vcv,
    .f = ~ exp(.x + diag(.y)/2)  ## analytic form for the Poisson model, see annex villemereuil
  )) %>% 
  ### we then average within metapops:
  mutate(mean_patch = map(.x = patch_avgs, .f = ~mean(.x))) %>% 
  mutate(mean_patch_corners = map(.x = patch_avgs, .f = ~mean(.x[c(1,3,7,9)]))) %>% 
  mutate(mean_patch_center = map(.x = patch_avgs, .f = ~mean(.x[c(5)]))) %>% 
  mutate(mean_patch_sides = map(.x = patch_avgs, .f = ~mean(.x[c(2,4,6,8)]))) %>% 
  unnest(c(mean_patch, mean_patch_corners,mean_patch_center,mean_patch_sides)) 

### then, the metapop-level predictions

tab<- tab %>% 
  mutate(mean_metapop = exp(fixef_metapop + (latent_rpl_var_metapop + latent_time_var_metapop)/2))

### the next step is the estimation of the observed scale temporal VCV to estimate alpha, beta and gamma

### in theory, for that (and the mean prediction above btw), we should use the QGmvparams function
### but runtime and memory use increase massively with number of variables (here 9), become unusable at about 5 variables
### by chance all of our variables are from the same distribution, and an easy one, so we can do things manually:

tab <- tab %>% 
  mutate(psiM = map2(.x = latent_intercepts, .y = vcv_metapop_latent,
                     .f = function(.x,.y){
                       psiM <-NA
                       for (j in 1:9) {
                         psiM[j] <- QGpsi(.x[j], .y[j, j], 
                                          d.link.inv = function(x) { exp(x)})
                       }
                       diag(psiM)
                     }  
  )) %>% 
  mutate(vcv_metapop_obs = map2(.x = psiM, .y= vcv_metapop_latent,
                                .f = ~ (.x %*% .y %*% t(.x))
  ))

tab <- tab %>% 
  mutate(
    alpha = map2(.x = vcv_metapop_obs, .y = patch_avgs, .f = ~.x %>% alpha_wang_loreau(varcorr=., means = .y)),
    gamma = map2(.x = vcv_metapop_obs, .y = patch_avgs, .f = ~.x %>% gamma_wang_loreau(varcorr=., means = .y))
  ) %>% 
  unnest(cols=c(alpha,gamma)) %>% 
  mutate(
    beta1 = alpha / gamma, beta2 = alpha - gamma
  )

tab %>% group_by(.iteration,LENGTH, SHUFFLE) %>% summarise(mean = mean(mean_patch)) %>% 
  ggplot()+stat_halfeye(aes(x=mean,y=factor(LENGTH)))+ scale_x_log10()+ facet_wrap(~SHUFFLE)



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


fits %>% 
  group_by(LENGTH, SHUFFLE, .iteration) %>% 
  mutate(TOTALN = P11 + P12 + P13 + P21 + P22 + P23 + P31 + P32 + P33) %>% 
  summarise(meanTOTAL = mean(TOTALN)) %>% 
  ggplot()+
  stat_halfeye(aes(x=meanTOTAL/9,y=factor(as.numeric(LENGTH))))+
  facet_wrap(~SHUFFLE)
### THIS IS HOW YOU FALL BACK ON THE MANUSCRIPT ESTIMATES

### correct based on arithmetic means

### naive use of model coefs gives geometric means of patches
