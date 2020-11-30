
## columns prefixed by P_ refer to the patch by patch model
## by M_ to the metapop-level model
## rpl something related to the among-replicates level of variation
## time to the within-replicate (temporal) level of variation

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
         SHUFFLE=str_extract(METAPOP_ID,"NO|R")) %>% 
  mutate(TREATMENT=interaction(LENGTH,SHUFFLE))


### now we need to extract the temporal variance-covariance matrices:
memory.limit(40000) #### if needed, a little memory boost to handle the full vcv matrix in final model
P_vcv_time_global <- VarCorr(mod, summary = FALSE)$WEEKrpl$cov

P_vcv_time_latent<-P_ranef_rpl %>% 
  select(.iteration,METAPOP_ID) %>% 
  mutate(P_vcv_time_latent=map2(.x=.iteration,.y=METAPOP_ID,
                                .f=function(iter=.x,metapop=.y,source=P_vcv_time_global){
                                  include <- colnames(source)
                                  include <- include[str_detect(include, pattern = metapop)]
                                  return(source[iter,include,include])
                                }))

rm(P_vcv_time_global);gc() ##we remove the big global VCV from memory since we don't need it anymore

P_vcv_rpl <- VarCorr(mod, summary = FALSE)$METAPOP_ID$sd %>% 
  array_tree(margin=1) %>% 
  tibble() %>% 
  rename(., P_vcv_rpl_latent=".") %>% 
  mutate(.iteration = 1:dim(.)[1]) %>% 
  mutate(P_vcv_rpl_latent=map(.x=P_vcv_rpl_latent,
                       .f=function(.x){diag(.x^2)}))



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


M_var_time <- VarCorr(mod,summary=FALSE)$WEEKrpl2$sd %>%
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
  left_join(P_vcv_rpl) %>% 
  left_join(P_vcv_time_latent) %>% 
  left_join(M_ranef_rpl) %>% 
  left_join(M_var_rpl) %>% 
  left_join(M_var_time) %>%  
  left_join(P_M_fixef) %>% 
  mutate(P_latent_intercept = map2(
    .x = P_fixef, .y = P_ranef_rpl,
    .f= function(.x,.y){.x + .y}
  ))



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

tab <- tab %>% 
  mutate(P_vcv_total_latent = map2(.x = P_vcv_time_latent, .y = P_vcv_rpl_latent,
                          .f = function(.x,.y){.x+.y})) %>% 
  mutate(P_pred = map2( ## patch by patch average
    .x = P_fixef, .y = P_vcv_time_latent, # or vcv_total? or vcv_time_latent, or vcv_time_latent_mean...?
    .f = ~ exp(.x + diag(.y)/2)  ## analytic form for the Poisson model, see annex villemereuil
  )) %>% 
  mutate(P_pred_for_var = map2( ## patch by patch mean, not averaged over anything else, needed for alpha/beta/gamma vars
    .x = P_latent_intercept, .y = P_vcv_time_latent,
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
  mutate(M_pred = exp(M_fixef + (M_var_rpl_latent + M_var_time_latent)/2)) %>% 
  mutate(M_latent_intercept = M_fixef + M_ranef_rpl) %>% 
  mutate(M_alpha = map2(.x = M_latent_intercept, .y = M_var_time_latent,
                        .f = function(.x,.y){
                          M_obscale = QGparams(mu=.x, var.a = .y, var.p = .y, model="Poisson.log", verbose=FALSE)
                          CV_L = sqrt(M_obscale$var.a.obs)/M_obscale$mean.obs
                          return(CV_L^2)
                          })) %>% 
  unnest(M_alpha)

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
    P_alpha = map2(.x = P_vcv_time_obs, .y = P_pred_for_var, .f = ~.x %>% alpha_wang_loreau(varcorr=., means = .y)),
    P_gamma = map2(.x = P_vcv_time_obs, .y = P_pred_for_var, .f = ~.x %>% gamma_wang_loreau(varcorr=., means = .y))
  ) %>% 
  unnest(cols=c(P_alpha,P_gamma)) %>% 
  mutate(
    P_beta1 = P_alpha / P_gamma, P_beta2 = P_alpha - P_gamma
  )


obssummary <- data %>% 
  pivot_longer(P11:P33) %>% 
  group_by(SHUFFLE,LENGTH,METAPOP_ID,name) %>% ##what's the most relevat/correct avging for display?
  summarise(P_mean=mean(value), P_se=plotrix::std.error(value),P_median = median(value),
            M_mean=mean(METAPOPSUM),M_se=plotrix::std.error(METAPOPSUM))

tab %>% group_by(.iteration,LENGTH, SHUFFLE) %>% summarise(mean = mean(P_mean_all)) %>% 
  ggplot()+
  stat_halfeye(aes(x=mean,y=log(LENGTH,base=4)+0.05),
               orientation="horizontal",.width=c(0.01,0.95))+
  geom_pointinterval(
    data = obssummary, 
    aes(x=P_mean, xmin=P_mean-P_se,xmax=P_mean+P_se,y=log(LENGTH,base=4)-0.05), 
    col="grey40",size = 1.5, alpha=0.5, position=position_jitter(height=0.1,width=NULL))+
  scale_x_continuous("mean patch population size (adult females)")+
  scale_y_continuous("bridge length (cm)", breaks=log(c(4,8,16),base=4), labels=c(4,8,16))+
  facet_wrap(~SHUFFLE)+
  cowplot::theme_half_open(11) +
  cowplot::background_grid(colour.major = "grey95", colour.minor = "grey95")


tab %>% group_by(.iteration,LENGTH, SHUFFLE) %>% summarise(mean = mean(P_alpha)) %>% 
  ggplot()+
  stat_halfeye(aes(x=mean,y=log(LENGTH,base=4)),
               orientation="horizontal",.width=c(0.01,0.95))+
  scale_x_continuous(expression(paste("mean ", alpha, " variability")))+
  scale_y_continuous("bridge length (cm)", breaks=log(c(4,8,16),base=4), labels=c(4,8,16))+
  facet_wrap(~SHUFFLE)+
  cowplot::theme_half_open(11) +
  cowplot::background_grid(colour.major = "grey95", colour.minor = "grey95")

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



