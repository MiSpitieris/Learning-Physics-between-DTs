# experiments with toy model--different u-param
## uncomment to install packages
# install.packages("SAVE")
# install.packages("ggplot2")
# install.packages("rstan")
library(SAVE)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = 3)

M = function(u,x) 5*exp(-u*x) # misspecified model
R = function(u,x,b) 3.5*exp(-u*x)+b # model we simulate data

sd_noise = 0.3
data("synthfield") # data from the rpackage SAVE
X_loc = unique(synthfield$x) 

# simulate data for different u val and add iid noise
u_val = seq(0.8,1.7,by=0.1)
dl=list()
set.seed(123)
offsets=runif(length(u_val),0.5,5) # sample offsets in [0.5,5]
# input locations with 3 relpicates
xobs = c(unique(synthfield$x),unique(synthfield$x),unique(synthfield$x)); N=length(xobs);
X_mat = matrix(NA, nrow = length(u_val), ncol = length(xobs))
set.seed(0)
# sample random input locations 
dev=c(runif(length(u_val), 0.1,0.25))
for(i in 1:nrow(X_mat)){
  X_mat[i,] = xobs+dev[i]
}
# Simulate data from the true model R  
for(i in seq_along(u_val)){
  set.seed(0)
  y=R(u_val[i], xobs, offsets[i])+rnorm(N,0,sd=sd_noise); 
  dl[[i]] = list(x= X_mat[i,], y = y, N = N, x_pred=X_mat[i,], N_pred=N)
}
y = matrix(NA, nrow = length(u_val), ncol = length(xobs))
for(i in 1:nrow(y)) y[i,] = dl[[i]]$y

id = seq_along(u_val)
Ns = length(id)
data_population = list(x= X_mat, y = y, N = N, Ns=Ns, id=id)

#------------------------------------------------------------------------------------------
# Fit without accounting for discrepancy for each invidual separately (no-without delta))
lu_no=list()
for(i in seq_along(u_val)){
  fit_no_without_delta = stan(
    file='STAN/toy/toy_nodelta.stan', # without delta
    data=dl[[i]],
    chains=3,
    iter=2*1000,
    seed=123
  )
  lu_no[[i]] = extract(fit_no_without_delta)$u
}
stan_trace(fit_no_without_delta)

#------------------------------------------------------------------------------------------
# Accounting for discrepancy for each invidual separately (no-with delta)
lu=list()
for(i in seq_along(u_val)){
  fit_no_with_delta = stan(
    file='STAN/toy/toy_delta.stan', # with delta
    data=dl[[i]],
    chains=3,
    iter=2*1000,
    seed=123
  )
  lu[[i]] = extract(fit_no_with_delta)$u
}
stan_trace(fit_no_with_delta)

#------------------------------------------------------------------------------------------
# shared u common delta model
fit_yes_common_delta = stan(file='STAN/toy/toy_common_delta.stan',
                            data=data_population,
                            chains=3,
                            iter=2*1000,
                            seed=123
)
names(fit_yes_common_delta)
stan_trace(fit_yes_common_delta, pars = "u")
ex_ycd=extract(fit_yes_common_delta)

#------------------------------------------------------------------------------------------
# shared u and delta model
fit_yes_shared_delta = stan(file='STAN/toy/toy_shared_delta.stan',
                            data=data_population,
                            chains=3,
                            iter=2*1000,
                            seed=123
)
names(fit_yes_shared_delta)
stan_trace(fit_yes_shared_delta, pars = "u")
ex_ysd=extract(fit_yes_shared_delta)

nodelta_cis=m_nwd=matrix(NA,length(u_val),2)
fn_s = function(x) c(quantile(x, probs = c(0.025,0.975)))

for(i in seq_along(u_val)){
  nodelta_cis[i,] = fn_s(lu_no[[i]])
  m_nwd[i,] = fn_s(lu[[i]])
}

ysd_cis=data.frame(t(apply(ex_ysd$u,2,quantile,probs=c(0.025,0.975))))
yes_common_delta=data.frame(t(apply(ex_ycd$u,2,quantile,probs=c(0.025,0.975))))

m_nwd=data.frame(m_nwd)
colnames(m_nwd) = colnames(ysd_cis)
nodelta_cis=data.frame(nodelta_cis)
colnames(nodelta_cis) = colnames(ysd_cis)
m_nwd$u_true = u_val
m_nwd$sharing= "no-with delta"
nodelta_cis$u_true = u_val
nodelta_cis$sharing= "no-without delta"
yes_common_delta$u_true=u_val
yes_common_delta$sharing="yes/common delta"
ysd_cis$u_true=u_val
ysd_cis$sharing="yes/shared delta"

df_CIs = rbind(m_nwd,yes_common_delta,ysd_cis,nodelta_cis)
colnames(df_CIs)[1:2] = c("lower","upper")
df_CIs$id = rep(id,4)

(toy_cis=ggplot(df_CIs,aes(x = as.factor(id), y = u_true, ymin = lower, ymax = upper, color=sharing))+ 
    geom_errorbar(width = 0.3, size=0.9) +
    theme_classic()+ 
    geom_point(size = 1.5, color="black")+
    ylab("u")+xlab("ID")+
    annotate("text", x=1.5, y=4.2, label= "95% CIs", size=6)+
    theme(axis.title.y = element_text(size = rel(1.5)))+
    scale_color_manual(
      breaks=c('no-without delta', 'no-with delta', "yes/common delta", "yes/shared delta"),
      values=c("magenta","firebrick","dodgerblue4", "forestgreen")))
ggsave("figs/toy_cis.pdf", plot = toy_cis, width = 18, height = 8, units = "cm")
