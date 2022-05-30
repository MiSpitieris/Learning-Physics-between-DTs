# experiments with toy model--different u-param
## uncomment to install packages
# install.packages("SAVE")
# install.packages("ggplot2")
# install.packages("rstan")
library(SAVE)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-1)

M = function(u,x) 5*exp(-u*x)
R = function(u,x,b) 3.5*exp(-u*x)+b

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
# Predictions at 20 locations x_pred (both interpolation and extrapolation)
Ns = length(u_val) # total number of individuals
id = seq_along(u_val) # individual ids
N_pred=20
x_pred_vec=seq(0.1,5,length.out = N_pred);
x_pred = matrix(rep(x_pred_vec,Ns), nrow = Ns, ncol = N_pred, byrow = TRUE)
# add i.i.d. N(0,0.3^2) noise
for(i in seq_along(u_val)){
  set.seed(0)
  y=R(u_val[i], xobs, offsets[i])+rnorm(N,0,sd=sd_noise); 
  dl[[i]] = list(x= X_mat[i,], y = y, N = N, x_pred=x_pred_vec, N_pred=N_pred)
}
y = matrix(NA, nrow = length(u_val), ncol = length(xobs))
for(i in 1:nrow(y)) y[i,] = dl[[i]]$y


# population data
data_population = list(x= X_mat, y = y, N = N, Ns=Ns, id=id, N_pred = N_pred, x_pred=x_pred)

#------------------------------------------------------------------------------------------
# Fit without accounting for discrepancy for each invidual separately (no-without delta))
lu_no=lpred_no=list()
for(i in seq_along(u_val)){
  fit_no_without_delta = stan(
    file='STAN/toy/toy_nodelta_pred.stan', # without delta
    data=dl[[i]],
    chains=3,
    iter=2*1000,
    seed=123
  )
  lu_no[[i]] = extract(fit_no_without_delta)$u
  smr_nnd=summary(fit_no_without_delta)$summary
  ind=grep("y_pred",rownames(smr_nnd))
  lpred_no[[i]]=smr_nnd[ind,c("mean", "2.5%", "97.5%")]
}
stan_trace(fit_no_without_delta)
pred_nnd=data.frame(do.call(rbind, lpred_no))
pred_nnd$sharing = "no-without delta"

#------------------------------------------------------------------------------------------
# Accounting for discrepancy for each invidual separately (no-with delta)
lu=lpred=list()
for(i in seq_along(u_val)){
  fit_no_with_delta = stan(
    file='STAN/toy/toy_delta_pred.stan', # with delta
    data=dl[[i]],
    chains=3,
    iter=2*1000,
    seed=123
  )
  lu[[i]] = extract(fit_no_with_delta)$u
  smr_nwd=summary(fit_no_with_delta)$summary
  ind=grep("y_pred",rownames(smr_nwd))
  lpred[[i]]=smr_nwd[ind,c("mean", "2.5%", "97.5%")]
}
stan_trace(fit_no_with_delta)
pred_nwd=data.frame(do.call(rbind, lpred))
pred_nwd$sharing = "no-with delta"

#------------------------------------------------------------------------------------------
# shared u common delta model
fit_yes_common_delta = stan(file='STAN/toy/toy_common_delta_pred.stan',
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
fit_yes_shared_delta = stan(file='STAN/toy/toy_shared_delta_pred.stan',
                            data=data_population,
                            chains=3,
                            iter=2*1000,
                            seed=123
)
smr_ysd= summary(fit_yes_shared_delta)
stan_trace(fit_yes_shared_delta, pars = "u")
ex_ycd=extract(fit_yes_shared_delta)
# Extract means and quantiles
indy=grep("y_", rownames(smr_ysd$summary))
pred_ysd=data.frame(smr_ysd$summary[indy,c("mean", "2.5%", "97.5%")])
pred_ysd$sharing = "yes/shared delta"
smr_ycd= summary(fit_yes_common_delta)
indy=grep("y_", rownames(smr_ycd$summary))
pred_ycd=data.frame(smr_ycd$summary[indy,c("mean", "2.5%", "97.5%")])
pred_ycd$sharing = "yes/common delta"

df_pred = rbind( pred_nwd, pred_ysd, pred_nnd,pred_ycd)

df_pred$x = rep(as.vector(t(x_pred)),4)
colnames(df_pred)[2:3] = c("lower", "upper")

ID = paste("ID = ", id)
bs = paste("| b =", round(offsets,2))
ID=paste(ID, bs)
df_obs = data.frame(x=as.vector(data_population$x), y=as.vector(data_population$y), id=rep(ID,N))
df_pred$id = rep(rep(ID,each=N_pred),4) 

y_true=matrix(NA,nrow = nrow(y), ncol=N_pred)
for (i in 1:Ns) {
  y_true[i,]=R(u_val[i], sort(x_pred[i,]), offsets[i])
}
df_true = data.frame(x=as.vector(t(x_pred)), y=as.vector(t(y_true)), id=rep(ID, each=N_pred))
data_fig_toy_pred = list(df_true=df_true, df_obs=df_obs, df_pred=df_pred, data_population=data_population)
# saveRDS(data_fig_toy_pred, file = "Data/data_fig_toy_pred.rds")
pl_toy_pred=ggplot()+
  geom_vline(xintercept = max(X_mat), linetype=2, colour="grey68", size=0.7)+
  geom_line(data=df_true, aes(x=x, y=y, linetype="true"))+
  geom_point(data=df_obs, aes(x=x, y=y, shape="observed"))+
  geom_line(data=df_pred, aes(x=x, y=mean, color= sharing), size=1)+
  geom_ribbon(data = df_pred, aes(x=x,ymin=lower,ymax=upper, fill=sharing),alpha=0.28)+
  facet_wrap(~id, nrow = 2)+
  theme_bw()+
  theme(legend.position = "bottom", legend.title = element_blank())+
  scale_color_manual(
    breaks=c('no-without delta','no-with delta', "yes/common delta", "yes/shared delta"),
    values=c("magenta","firebrick","dodgerblue4", "forestgreen"))+
  scale_fill_manual(
    breaks=c('no-without delta','no-with delta', "yes/common delta", "yes/shared delta"),
    values=c("magenta","firebrick","dodgerblue4", "forestgreen"))
pl_toy_pred
ggsave("Figs/toy_pred.pdf", plot = pl_toy_pred, width = 20, height = 12, units = "cm")

### create posterior densities plot
ex_ycd=extract(fit_yes_common_delta)
ex_ysd=extract(fit_yes_shared_delta)

nnd_u = nwd_u = c()
for(i in 1:Ns){
  nwd_u = c(nwd_u,lu[[i]])
  nnd_u = c(nnd_u,lu_no[[i]])
}

n_sam = nrow(ex_ycd$u)
post_nnd = data.frame(u = nnd_u, id = rep(id,each=n_sam), sharing='no-without delta')
post_nwd = data.frame(u = nwd_u, id = rep(id,each=n_sam), sharing='no-with delta')
post_ycd = data.frame(u = as.vector(ex_ycd$u), id = rep(id,each=n_sam), sharing='yes/common delta')
post_ysd = data.frame(u = as.vector(ex_ysd$u), id = rep(id, each = n_sam), sharing='yes/shared delta')
post_u_df = rbind( post_nnd,post_nwd, post_ycd, post_ysd)
u_true_df= data.frame(u = u_val, id = as.factor(id))
data_post_toy = list(post_u_df = post_u_df, u_true_df = u_true_df)
# saveRDS(data_post_toy, file = "Data/data_post_toy.rds")
(post_toy=ggplot() + 
  geom_violin(data = post_u_df, aes(x = as.factor(id), y=u, color = sharing))+
  theme_bw()+ 
  theme(axis.title.y = element_text(size = rel(1.5)), legend.position = "none")+
  geom_point(data = u_true_df, aes(x = id, y = u, shape = "true"))+
  facet_wrap(sharing~.)+
  xlab("ID")+
  scale_color_manual(
    breaks=c('no-without delta', 'no-with delta', "yes/common delta", "yes/shared delta"),
    values=c("magenta","firebrick","dodgerblue4", "forestgreen")))
ggsave("Figs/post_toy.pdf", plot = post_toy, width = 18, height = 12, units = "cm")
