

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




obssummary <- data %>% 
  pivot_longer(P11:P33) %>% 
  group_by(SHUFFLE,LENGTH,METAPOP_ID,name) %>% ##what's the most relevat/correct avging for display?
  summarise(P_mean=mean(value), P_se=plotrix::std.error(value),P_median = median(value),
            M_mean=mean(METAPOPSUM),M_se=plotrix::std.error(METAPOPSUM))

P_summary %>% group_by(.iteration,LENGTH, SHUFFLE) %>% summarise(mean = mean(P_mean_all)) %>% 
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


M_summary %>% group_by(.iteration,LENGTH, SHUFFLE) %>% summarise(mean = mean(M_mean)) %>% 
  ggplot()+
  stat_halfeye(aes(x=mean,y=log(LENGTH,base=4)+0.05),
               orientation="horizontal",.width=c(0.01,0.95))+
  geom_pointinterval(
    data = obssummary %>% select(SHUFFLE, LENGTH, M_mean,M_se) %>% distinct(), 
    aes(x=M_mean, xmin=M_mean-M_se,xmax=M_mean+M_se,y=log(LENGTH,base=4)-0.05), 
    col="grey40",size = 1.5, alpha=0.5, position=position_jitter(height=0.1,width=NULL))+
  scale_x_continuous("mean metapopulation size (adult females)")+
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



