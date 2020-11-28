
## columns prefixed by P_ refer to the patch-level model
## by M_ to the metapop-level model
## rpl something related to the among replicates variation
##time to the temporal variance

### we first extract all the fixed effects:
P_M_fixef <- fixef(mod, summary = FALSE) %>% ### all model fixed effects, with one column = one "variable" * treatment combination
  as_tibble() %>%
  mutate(.iteration = 1:dim(.)[1]) %>% ### we put the implicit iteration # into an explicit column
  pivot_longer(-.iteration) %>%        ### we pivot so we can split the "variable*treatment" into patch ID and treatment variables
  mutate(variable=str_extract(name,"P11|P12|P13|P21|P22|P23|P31|P32|P33|METAPOPSUM"),
         LENGTH=as.numeric(str_extract(name,"4|8|16")),
         SHUFFLE=str_extract(name,"NO$|R$")) %>%  # $ : end of string
  select(-name) %>%
  pivot_wider(values_from=value,names_from=variable) %>%  ###we repivot back so that 1 column per patch
  group_by(.iteration, LENGTH, SHUFFLE) %>%  ### we group and nest by iteration (with treatment variables along for the ride)
  nest(P_fixef = c(P11:P33)) %>% 
  mutate(P_fixef = map(.x=P_fixef, .f = ~.x %>% unlist())) %>%   #convert list column from tibbles to vectors for faster calculations
  rename(M_fixef = METAPOPSUM)

### we then extract the among-metapop deviations for the patch by patch model (needed for estimating alpha beta gamma):
P_ranef_rpl<- ranef(mod, summary = FALSE)$METAPOP_ID %>% 
  array_tree(margin=c(1,2)) %>% ## convert array to list, enlisting by metapop nested in iteration
  tibble() %>% ## we "flatten" the list, iteration become rows of new tibble
  rename(., P_ranef_rpl=".") %>% 
  mutate(.iteration = 1:dim(.)[1]) %>% ### we put the implicit iteration # into an explicit column
  unnest_longer(P_ranef_rpl, indices_to = "METAPOP_ID") %>% 
  mutate(LENGTH=as.numeric(str_extract(METAPOP_ID,"^4|^8|^16")),  ## ^: start of string
         SHUFFLE=str_extract(METAPOP_ID,"NO|R"))

### we do the same for the metapop sum model:
M_ranef_rpl<- ranef(mod, summary = FALSE)$METAPOP_ID2 %>%
  array_tree(margin=c(1)) %>% ## convert array to list, enlisting by metapop nested in iteration
  tibble() %>% ## we "flatten" the list, iteration become rows of new tibble
  rename(., M_ranef_rpl=".") %>% 
  mutate(M_ranef_rpl = map(.x = M_ranef_rpl,
                                     .f = ~ t(.x) %>% as.data.frame() %>% c())) %>% 
  mutate(.iteration = 1:dim(.)[1]) %>% ### we put the implicit iteration # into an explicit column
  unnest_longer(M_ranef_rpl, indices_to = "METAPOP_ID") %>% 
  mutate(LENGTH=as.numeric(str_extract(METAPOP_ID,"^4|^8|^16")),  ## ^: start of string
         SHUFFLE=str_extract(METAPOP_ID,"NO|R"))

### now we need to extract the temporal variance-covariance matrices:
P_vcv_time_global <- VarCorr(mod, summary = FALSE)$WEEK2$cov %>% 
      array_tree(margin=1) %>% 
      tibble() %>% 
  rename(., P_vcv_time_latent_global=".") %>% 
  ### this is the VCV matrix including all the patches 
  ### so including the zero covs between the patchs from one metapop and the patch from another
  mutate(.iteration = 1:dim(.)[1])

P_vcv_rpl <- VarCorr(mod, summary = FALSE)$METAPOP_ID$cov %>% 
  array_tree(margin=1) %>% 
  tibble() %>% 
  rename(., P_vcv_rpl_latent=".") %>% 
  mutate(.iteration = 1:dim(.)[1])


M_var_time <- VarCorr(mod,summary=FALSE)$WEEK3$sd %>%
  as.data.frame() %>% 
  mutate(.iteration = 1:dim(.)[1]) %>% 
  pivot_longer(-.iteration) %>% 
  mutate(M_var_time_latent = value^2) %>% 
  mutate(METAPOP_ID = str_remove(name,"METAPOPSUM_Intercept:METAPOP_ID2")) %>% 
  select(METAPOP_ID, M_var_time_latent, .iteration) %>% 
  mutate(LENGTH=as.numeric(str_extract(METAPOP_ID,"^4|^8|^16")),  ## ^: start of string
         SHUFFLE=str_extract(METAPOP_ID,"NO|R"))

M_var_rpl <- VarCorr(mod,summary=FALSE)$METAPOP_ID2$sd %>% 
  as.data.frame() %>% 
  mutate(.iteration = 1:dim(.)[1]) %>% 
  mutate(M_var_rpl_latent = METAPOPSUM_Intercept^2) %>% 
  select(.iteration,M_var_rpl_latent)

### and we merge all these objects into one table:
### this will in addition help us split the global VCV matrix into one VCV matrix per metapop:
tab <- P_ranef_rpl %>% 
  left_join(M_ranef_rpl) %>% 
  left_join(M_var_time) %>%  
  left_join(M_var_rpl) %>% 
  left_join(P_M_fixef) %>% 
  left_join(P_vcv_rpl) %>% 
  inner_join(P_vcv_time_global, by = ".iteration") %>% 
  mutate(P_vcv_time_latent = map2(
    .x = P_vcv_time_latent_global, .y = METAPOP_ID,
    .f = ~.x %>% .[str_detect(colnames(.x),.y),str_detect(colnames(.x),.y)]
  )) %>% 
  mutate(P_latent_intercept = map2(
    .x = P_fixef, .y = P_ranef_rpl,
    .f= function(.x,.y){.x + .y}
  ))

###let's reorder the existing columns a bit, just to see better (and remove the few we don't need anymore)

tab <- tab %>% 
  select(c(.iteration,METAPOP_ID,LENGTH, SHUFFLE, 
           P_fixef, P_ranef_rpl, P_latent_intercept, P_vcv_time_latent, P_vcv_rpl_latent,
           M_fixef, M_ranef_rpl, M_var_time_latent, M_var_rpl_latent)
           )

### we can now start to calculate quantities of interest, replicate by replicate

### first, the predicted mean patch sizes:

tab <- tab %>% 
  mutate(P_vcv_total_latent = map2(.x = P_vcv_time_latent, .y = P_vcv_rpl_latent,
                          .f = function(.x,.y){.x+.y})) %>% 
  mutate(P_pred = map2( ## patch by patch average
    .x = P_fixef, .y = P_vcv_time_latent, # or vcv_total?
    .f = ~ exp(.x + diag(.y)/2)  ## analytic form for the Poisson model, see annex villemereuil
  )) %>% 
  ### we then average within metapops:
  mutate(P_mean_all = map(.x = P_pred, .f = ~mean(.x))) %>% 
  mutate(P_mean_corners = map(.x = P_pred, .f = ~mean(.x[c(1,3,7,9)]))) %>% 
  mutate(P_mean_center = map(.x = P_pred, .f = ~mean(.x[c(5)]))) %>% 
  mutate(P_mean_sides = map(.x = P_pred, .f = ~mean(.x[c(2,4,6,8)]))) %>% 
  unnest(c(P_mean_all, P_mean_corners,P_mean_center,P_mean_sides)) 

### then, the metapop-level predictions

tab<- tab %>% 
  mutate(M_pred = exp(M_fixef + (M_var_rpl_latent + M_var_time_latent)/2))

### the next step is the estimation of the observed scale temporal VCV to estimate alpha, beta and gamma

### in theory, for that (and the mean prediction above btw), we should use the QGmvparams function
### but runtime and memory use increase massively with number of variables (here 9), become unusable at about 5 variables
### by chance all of our variables are from the same distribution, and an easy one, so we can do things manually:

tab <- tab %>% 
  mutate(P_psi = map2(.x = P_latent_intercept, .y = P_vcv_time_latent,
                     .f = function(.x,.y){
                       psi <-NA
                       for (j in 1:9) {
                         psi[j] <- QGpsi(.x[j], .y[j, j], 
                                          d.link.inv = function(x) { exp(x)})
                       }
                       diag(psi)
                     }  
  )) %>% 
  mutate(P_vcv_time_obs = map2(.x = P_psi, .y= P_vcv_time_latent,
                                .f = ~ (.x %*% .y %*% t(.x))
  ))

tab <- tab %>% 
  mutate(   ##very important !! go up to check the right components are in the P_pred above
    ## if not, duplicate the column and include the right ones
    P_alpha = map2(.x = P_vcv_time_obs, .y = P_pred, .f = ~.x %>% alpha_wang_loreau(varcorr=., means = .y)),
    P_gamma = map2(.x = P_vcv_time_obs, .y = P_pred, .f = ~.x %>% gamma_wang_loreau(varcorr=., means = .y))
  ) %>% 
  unnest(cols=c(P_alpha,P_gamma)) %>% 
  mutate(
    P_beta1 = P_alpha / P_gamma, P_beta2 = P_alpha - P_gamma
  )

tab %>% group_by(.iteration,LENGTH, SHUFFLE) %>% summarise(mean = mean(M_pred)) %>% 
  ggplot()+stat_halfeye(aes(x=mean,y=factor(LENGTH)))+ facet_wrap(~SHUFFLE)



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
