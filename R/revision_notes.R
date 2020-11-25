
### notes on revamping the post-processing function using purrr

### make it clearer, more flexible, and (i think?) less memory hungry

FIXEF <- fixef(mod, summary = FALSE) %>% 
  as_tibble() %>% 
  mutate(.iteration = 1:dim(.)[1]) %>% 
  pivot_longer(-.iteration) %>% 
  mutate(patch=substr(name,1,3),    ## !!! TO ADAPT IF METAPOP SUM END UP IN MODEL
         LENGTH=str_extract(name,"4|8|16"),
         SHUFFLE=str_extract(name,"NO$|R$")) %>%  # $ : end of string
  mutate(TREATMENT = interaction(LENGTH,SHUFFLE)) %>% 
  select(-name) %>% 
  pivot_wider(values_from=value,names_from=patch) %>% 
  group_by(.iteration, LENGTH, SHUFFLE, TREATMENT) %>% 
  nest() %>% 
  rename(fixef=data)

RANEF_METAPOP<- ranef(mod, summary = FALSE)$METAPOP_ID %>% 
  array_tree(margin=c(1,2)) %>% 
  tibble() %>% 
  rename(., ranef_metapop=".") %>% 
  mutate(.iteration = 1:dim(.)[1]) %>% 
  unnest_longer(ranef_metapop, indices_to = "METAPOP_ID") %>% 
  mutate(ranef_metapop = map(.x=ranef_metapop, 
                             .f = ~.x %>% bind_rows()  ##change from vector to tibble with one col = one patch
                             )) %>%  
  mutate(LENGTH=str_extract(METAPOP_ID,"^4|^8|^16"),  ## ^: start of string
         SHUFFLE=str_extract(METAPOP_ID,"NO|R")) %>% 
  mutate(TREATMENT = interaction(LENGTH,SHUFFLE))

GLOBAL_VCV_time <- VarCorr(mod, summary = FALSE)$WEEK2$cov %>% 
      array_tree(margin=1) %>% 
      tibble() %>% 
  rename(., vcv_global_latent=".") %>% 
  ### the VCV including all the patches 
  ### so including the zero covs between the patchs from one metapop and the patch from another
  mutate(.iteration = 1:dim(.)[1])


test <- inner_join(RANEF_METAPOP,GLOBAL_VCV_time, by = ".iteration") %>% 
  mutate(vcv_metapop_latent = map2(
    .x = vcv_global, .y = METAPOP_ID,
    .f = ~.x %>% .[str_detect(colnames(.x),.y),str_detect(colnames(.x),.y)]
  )) %>% 
  left_join(FIXEF) %>% 
  mutate(latent_metapop_intercept = map2(
    .x = fixef, .y = ranef_metapop, ##the latent metapop level intercept is the sum of the global mean + the metapop level ranef
    .f=~ as_tibble(.x+.y))
  )



