library(tidyr)
library(ggplot2)
df <- sampled$samples$theta_mu[,8500:sampled$samples$idx]
tmp<-as.data.frame(t(df))
tmp2 <- pivot_longer(tmp, cols = everything(),names_to = "Parameter", values_to = "value")
tmp2$Iteration = rep(seq(8500:sampled$samples$idx),each = sampled$n_pars)
#ggplot(tmp2, aes(x=Iteration, y=value))+geom_line()+facet_wrap(~Parameter, nrow = sampled$n_pars)+theme_bw()
ggplot(tmp2, aes(x=Iteration, y=value))+geom_line()+facet_wrap(~Parameter, scales = "free")+theme_bw()
#ggplot(tmp2, aes(x=Iteration, y=value, colour = Parameter))+geom_line()+theme_bw()

df <- sampled$samples$theta_mu[,500:sampled$samples$idx]
tmp<-as.data.frame(t(df))
tmp2 <- pivot_longer(tmp, cols = everything(),names_to = "Parameter", values_to = "value")
tmp2$Iteration = rep(seq(500:sampled$samples$idx),each = sampled$n_pars)
#ggplot(tmp2, aes(x=Iteration, y=value))+geom_line()+facet_wrap(~Parameter, nrow = sampled$n_pars)+theme_bw()
ggplot(tmp2, aes(x=Iteration, y=value))+geom_line()+facet_wrap(~Parameter, scales = "free")+theme_bw()
#ggplot(tmp2, aes(x=Iteration, y=value, colour = Parameter))+geom_line()+theme_bw()


## for incompleted data
df <- adapted$samples$theta_mu[,13000:adapted$samples$idx]
tmp<-as.data.frame(t(df))
tmp2 <- pivot_longer(tmp, cols = everything(),names_to = "Parameter", values_to = "value")
tmp2$Iteration = rep(seq(13000:adapted$samples$idx),each = adapted$n_pars)
ggplot(tmp2, aes(x=Iteration, y=value))+geom_line()+facet_wrap(~Parameter, scales = "free")+theme_bw()
