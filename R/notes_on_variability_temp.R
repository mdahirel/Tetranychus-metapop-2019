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
