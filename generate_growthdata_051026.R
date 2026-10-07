# simulated growth-only parasite datasets, based on 

library(ggplot2)
library(lubridate)
library(dplyr)
library(ggpubr)
library(MASS)



output_name_str <- "fitCOV_2_MVNbscse_2108_gam6"
start_date <- today()
now <- Sys.time()
start_time <- format(now, "%H_%M_%S")

# combine these to make a folder name for this run
folder_name <- paste0(output_name_str, "_", start_date, "_", start_time)
dir.create(folder_name)

PD_initialisation =  function(ipl, 
                              iplMu, 
                              iplSigma,
                              nstages) {
  pAge =  sapply(X =  seq(1,nstages), mean =  iplMu, sd =  iplSigma, FUN =  dnorm)
  
  x0 =  round(ipl * (pAge)/sum(pAge))
  
  return(x0)}


simu_growth_1pt <- function(ipl, 
                            ipl_mu, 
                            ipl_sd,
                            pmf, 
                            tseq_50,
                            tseq_gamma,
                            soft_tmax=F, # binary model
                            tmax_gamma=100,
                            tmax_50=40,
                            sample_times # needs to be in h
){
  
  # max possible stage (h) should be far clear of the tmax_50
  #nStages = round(min(60, tmax_50+30))
  nStages=40
  if(soft_tmax){nStages <- 60}
  
  # create empty matrix to store the entire run (hourly, by age)
  tot_matrix <- matrix(NA, nrow=nStages, ncol= max(sample_times))
  
  # initialise the first entry at t0
  #print(length(PD_initialisation(ipl=ipl, iplMu=ipl_mu, iplSigma = ipl_sd, nstages = nStages)))
  tot_matrix[,1] <- PD_initialisation(ipl=ipl, iplMu=ipl_mu, iplSigma = ipl_sd, nstages = nStages)/5000 # need /5000 to go from count to permL
  
  # set up the shapes of both 'fost' transitions
  seq_curve <- ((1:nStages)^tseq_gamma)/(tseq_50^tseq_gamma+(1:nStages)^tseq_gamma)
   tmax_curve <- ((1:nStages)^tmax_gamma)/(tmax_50^tmax_gamma+(1:nStages)^tmax_gamma)
  
  if(soft_tmax){print("softtmax")
  for(i in 2:max(sample_times)){
    tot_matrix[2:nStages,i] <-  tot_matrix[(1:(nStages-1)),i-1] * (1-tmax_curve[1:(nStages-1)])
    tot_matrix[1,i] <-  pmf*sum(tot_matrix[,i-1] * (tmax_curve))
  }} else{

  for(i in 2:max(sample_times)){
    tot_matrix[2:nStages,i] <-  tot_matrix[(1:(nStages-1)),i-1] 
    tot_matrix[1,i] <-  pmf*tot_matrix[nStages,i-1] 
  }
    }
  
  tot_cir <- matrix(1-seq_curve, nrow=1) %*% tot_matrix[,sample_times]
  return(tot_cir)
  
}

simu_growth <- function(ipl, 
                        ipl_mu, 
                        ipl_sd,
                        pmf, 
                        tseq_50 = 26,
                        tseq_gamma = 10,
                        nPt = 30,
                        sample_times = 24*c(3,4,5.5,6,6.5,7,7.6,8,8.5,9),
                        spread_ipl =0.2,
                        spread_ipl_mu =0.4,
                        spread_ipl_sd =0.2,
                        spread_pmf =0.2,
                        spread_tseq_50 =0.2,
                        spread_tseq_gamma =0.2,
                        soft_tmax=F,
                        tmax_50=40,
                        tmax_gamma=100,
                        seed=123,
                        MVN=F
){
  set.seed(seed)
  param <- data.frame(parameter = c("ipl",   "ipl_mu",   "ipl_sd",  "pmf", 
                                    "tseq_50","tseq_gamma" ),
                      location = c(ipl, ipl_mu, ipl_sd, pmf, tseq_50, tseq_gamma),
                      spread = c(spread_ipl, spread_ipl_mu, spread_ipl_sd, 
                                 spread_pmf, 
                                 spread_tseq_50, spread_tseq_gamma))
  
  
  if(MVN){
    # note for this file we are rmeoving tmax_gamma and tmax_50
    
    #First set up distribution for eta, patient variations from transformed pop avg
    omega = diag(c(spread_ipl, spread_ipl_mu, spread_ipl_sd, 
                   spread_pmf, 
                   spread_tseq_50, spread_tseq_gamma))

    #create correlation matrix from L:
    L = matrix(c((runif(36,0,1))), nrow=6)
    
    #make L a lower triangular matrix by removing all values above the diagonal
    L[1,2:6] <- L[2,3:6] <- L[3,4:6] <- L[4,5:6] <- L[5,6] <- 0
    
    #Eta will be MVN w/ mean 0 and covariance Sigma:
    Sigma = omega%*%((L)%*%t(L))%*%omega
    print(Sigma)
    
   # eta = mvrnorm(nPt,rep(0,6),Sigma)
    eta = mvrnorm(nPt,rep(0,6),cov_ext_mod3)
    # Input desired population average values
    pop_avg = c(ipl, ipl_mu, ipl_sd, pmf, tseq_50, tseq_gamma)
    
    a = 3*pop_avg #upper bounds double estimates
    b = 0.25*pop_avg #lower bounds half estimates
    
    #If we want each row to be 1 patient, repeat the vector 8 times
    a_matrix = matrix(rep(a,nPt), nrow=nPt, byrow=T)
    b_matrix = matrix(rep(b,nPt), nrow=nPt, byrow=T)
    
    # These calculations are essentially element-wise (cant divide vector by vector)
    phi_i <- matrix(rep((log((pop_avg - a)/(b-pop_avg))),nPt),nrow=nPt,byrow=T)  + eta
    
    # This is the inverse of the transformation in report section 2.3
    samp_paraMVN = (b_matrix*exp(phi_i)+a_matrix)/(exp(phi_i)+1)
    print( samp_paraMVN )
    samples <- as.vector(t(samp_paraMVN)) # t so we have each patient's full set of parameters then the next. 
    # print(samp_paraMVN)
    # print(samples)
  }else{
  
  samples <- as.vector(t(sapply(1:nrow(param), function(i){ param$location[i]*exp(rnorm(nPt, 0, param$spread[i]))})))
 
  }
  
    sampled <- data.frame(id = rep(1:nPt, each=nrow(param)),
                        param = rep(param$parameter,nPt),
                        val = samples)
  # bit inefficient here
  para_cir <- sapply(1:nPt, function(i){#print(sampled[sampled$id==i,])
    simu_growth_1pt(ipl=sampled$val[sampled$id==i & sampled$param=="ipl"], 
                    ipl_mu=sampled$val[sampled$id==i & sampled$param=="ipl_mu"], 
                    ipl_sd=sampled$val[sampled$id==i & sampled$param=="ipl_sd"], 
                    pmf=sampled$val[sampled$id==i & sampled$param=="pmf"], 
                    soft_tmax=soft_tmax,
                    tmax_50=tmax_50, 
                    tmax_gamma=tmax_gamma,
                  tseq_50=sampled$val[sampled$id==i & sampled$param=="tseq_50"], 
                    tseq_gamma=sampled$val[sampled$id==i & sampled$param=="tseq_gamma"], 
                    sample_times = sample_times)})
  
  simu <- data.frame(id = rep(1:nPt, each=length(sample_times)),
                     time = rep(sample_times, nPt),
                     para_cir = as.vector((para_cir)))
  
  return(list(param = param, 
              sampled = sampled,
              simu = simu ))
  
  
}

lrg_combo <- data.frame(scenario = NA, id = NA, time = NA, para_cir = NA,
                        ipl = NA,
                        ipl_mu = NA,
                        ipl_sd = NA,
                        pmf = NA,
                        tseq_50 = NA,
                        tseq_gamma = NA,
                        nPt = NA,
                        spread_ipl = NA,
                        spread_ipl_mu = NA,
                        spread_ipl_sd  = NA,
                        spread_pmf = NA,
                        spread_tseq_50  = NA,
                        spread_tseq_gamma  = NA )

setup <- expand.grid(ipl = c(100),
                     ipl_mu = c(1, 5, 10), 
                     ipl_sd = c(1, 5), 
                     pmf = c(10, 20), 
                     tseq_50 = c(22),
                     tseq_gamma = c(6),
                     nPt = 100,
                     spread_ipl = 0.1,
                     spread_ipl_mu =0.1,
                     spread_ipl_sd  = 0.5,
                     spread_pmf = 0.2,
                     spread_tseq_50  = 0.1,
                     spread_tseq_gamma  = 1,
                     soft_tmax=F, tmax_gamma=100, tmax_50=40)

#cov_ext_mod3 <- readRDS("cov_ext_mod3_v2.rds") # this is from a fit to model 3 of the data, L*L' *diag(omega)


setup$scenario = 1:nrow(setup)
write.csv(setup, paste0(folder_name,"/setup.csv"))
# for a FIXED universal sample times
sample_times_data = 24*c(3,4,4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8, 8.5,9)
#sample_times = 24*c(4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.5, 7, 7.5, 8,9)
sample_times =sample_times_data

for(i in 1:nrow(setup)){
  #for(i in 69){
  a <- simu_growth(ipl = setup$ipl[i],
                   ipl_mu = setup$ipl_mu[i],
                   ipl_sd = setup$ipl_sd[i],
                   pmf = setup$pmf[i],
                   tseq_50 = setup$tseq_50[i],
                   tseq_gamma = setup$tseq_gamma[i],
                   nPt= setup$nPt[i],
                   soft_tmax=setup$soft_tmax[i],
                   tmax_gamma = setup$tmax_gamma[i],
                   tmax_50=setup$tmax_50[i],
                   sample_times = sample_times,
                   spread_ipl = setup$spread_ipl[i],
                   spread_ipl_mu = setup$spread_ipl_mu[i],
                   spread_ipl_sd = setup$spread_ipl_sd[i],
                   spread_pmf = setup$spread_pmf[i],
                   spread_tseq_50= setup$spread_tseq_50[i],
                   spread_tseq_gamma = setup$spread_tseq_gamma[i],
                   seed=123)
  
  n_newrow <- setup$nPt[i]*length(sample_times)
  #this_sc =  data.frame(rep(setup[i,], n_newrow))
  this_sc = setup %>% slice(rep(i, n_newrow))
  lrg_combo <- rbind(na.omit(lrg_combo), cbind(as.data.frame(a[["simu"]]), as.data.frame(this_sc)))
}


write.csv( lrg_combo,paste0(folder_name,"/growth_data_simulated", today(), "_",format(Sys.time(), "%H_%M_%S"), "_numsc_", nrow(setup), ".csv"),row.names=F)
simu <- ggplot(lrg_combo[lrg_combo$id==1,], aes(x=time/24, y=log10(para_cir), col=as.factor(tseq_50),
                                                 linetype=as.factor(tseq_gamma), group=interaction(scenario,id))) + geom_point() + 
    geom_line() + facet_grid(ipl_mu~ipl_sd) + theme_light() + ggtitle(paste0("Patient 1 across ", nrow(setup), " simulation scenarios. "))


title <- paste0(colnames(setup)[c(8:13)],": ",setup[1,c(c(8:13))], ". ")
title <- paste(title[1], title[2], title[3], title[4], title[5], title[6])


noise <- 0.3
# lets add observaion noise!
set.seed(123)
lrg_combo$para_cir_0noise <- lrg_combo$para_cir
lrg_combo$para_cir_noise <- lrg_combo$para_cir*(rnorm(nrow(lrg_combo), 1, noise))
lrg_combo_trim <- lrg_combo[lrg_combo$time%in%c(sample_times),]


lrg_combo_trim_s <- lrg_combo_trim %>% filter(time%in%c(24*c(4, 5, 5.5, 6, 6.5, 7, 7.5, 8))) %>% group_by(time) %>% summarise(num = length(para_cir))
ggplot(lrg_combo_trim[lrg_combo_trim$time%in%c(24*c(4, 5, 5.5, 6, 6.5, 7, 7.5, 8)),], aes(x=time/24,
                                   y=log10(para_cir_noise), 
                                   group=id, 
                                   col=as.factor(id))) + ylim(0,5)+
                    geom_line(lwd=0.2) + theme_light() +xlim(3,9)+ 
                    theme(legend.position = "none") + labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia")+
                    theme(text=element_text(size=20)) + facet_wrap(~scenario)


ggplot(lrg_combo_trim[lrg_combo_trim$time%in%c(24*c(4, 5, 5.5, 6, 6.5, 7, 7.5, 8)),], aes(x=time/24,
                                                                                          y=log10(para_cir_noise), 
                                                                                          group=interaction(id,pmf), 
                                                                                          col=as.factor(pmf))) + ylim(0,5)+
  geom_line(lwd=0.2) + theme_light() +xlim(3,9)+ 
  theme(legend.position = "none") + labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia")+
  theme(text=element_text(size=20)) + facet_grid(ipl_sd~ipl_mu)



ggplot(lrg_combo_trim[lrg_combo_trim$time%in%c(24*c(4, 5, 5.5, 6, 6.5, 7, 7.5, 8)),], aes(x=time/24,
                                                                                          y=log10(para_cir_noise), 
                                                                                          group=interaction(id,pmf), 
                                                                                          col=as.factor(pmf))) + ylim(0,5)+
  geom_line(lwd=0.2) + theme_light() +xlim(3,9)+ 
  theme(legend.position = "none") + labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia")+
  theme(text=element_text(size=20)) + facet_grid(ipl_sd~ipl_mu)








ggsave(paste0(folder_name,"/growth_data_simulated_noise_dataplot_", today(), "_",format(Sys.time(), "%H_%M_%S"), "_numsc_", nrow(setup), ".png"),
       width=10, height=12)



setwd("/data_ext_4TB/PMF_est/aug_25_2026")
data <- read.csv("para.csv")
head(data)
data_full <- data
data <- data # tets firts 400 obs only?
# ggplot(data, aes(x=Time,
#                  y=log10.av.ppml, 
#                  group=Inoculum.Coh.Subject.ID, 
#                  col= Inoculum.Coh.Subject.ID)) + geom_line() + theme_light() +
#   theme(legend.position = "none")

data <- data[order(data$Inoculum.Coh.Subject.ID),]
data <- data[!is.na(data$log10.av.ppml),] # dropped 1

data$id <- as.numeric(factor(data$Inoculum.Coh.Subject.ID, labels=1:length(unique(data$Inoculum.Coh.Subject.ID))))
data$log10perml <- data$log10.av.ppml
#data$Time <- data$Time - 1 # we adjust the time back to consider no 1st cycle??

# list.files(pattern="growth_data_simu")
# lrg_combo <- read.csv("growth_data_simulated2026-07-24_10_58_05_numsc_288.csv")
prop_inside_dat <- data %>% group_by(Time) %>% summarise(max = max(log10perml),
                                                         min = min(log10perml),
                                                         q25 = quantile(log10perml, 0.25),
                                                         q75 = quantile(log10perml, 0.75),
                                                         q50 = quantile(log10perml, 0.5),
                                                         q2.5 = quantile(log10perml, 0.025),
                                                         q97.5 = quantile(log10perml, 0.975) )

dat_quantiles <- lrg_combo_trim %>% group_by(time) %>% summarise(max = max(log10(para_cir)),
                                                         min = min(log10(para_cir)),
                                                         q25 = quantile(log10(para_cir), 0.25),
                                                         q75 = quantile(log10(para_cir), 0.75),
                                                         q50 = quantile(log10(para_cir), 0.5),
                                                         q2.5 = quantile(log10(para_cir), 0.025),
                                                         q97.5 = quantile(log10(para_cir), 0.975) )


simu2 <- ggplot(lrg_combo_trim,  aes(x=time/24, y=log10(para_cir))) + geom_point(cex=0.6) + 
    geom_line(lwd=0.3, col="darkgreen", aes(group=id)) + facet_wrap(~scenario) + theme_light() + ylim(0,5) + xlim(c(3, 9))+
    labs(title="Testing variability, blue is true data min-max range:", subtitle=paste0(title))+
    geom_segment(data=prop_inside_dat, aes(x=Time, y=max, yend=min), col="blue", lwd=2, alpha=0.4) +
    geom_segment(data=prop_inside_dat,aes(x=Time-0.2, xend=Time+0.2, y=max), col="blue", lwd=1, alpha=0.5) +
    geom_segment(data=prop_inside_dat,aes(x=Time-0.2, xend=Time+0.2, y=min), col="blue", lwd=1, alpha=0.5) 
# ggsave(paste0("varblty_test_060826/rowth_data_simulated", today(), "_",format(Sys.time(), "%H_%M_%S"), "_numsc_", nrow(setup), ".png"),simu2,
#        width=15, height=10)
(simu3 <- ggplot(lrg_combo_trim,  aes(x=time/24, y=log10(para_cir_noise))) + geom_point(cex=0.6) + 
    geom_line(lwd=0.3, col="darkgreen", aes(group=id), alpha=0.3) + facet_wrap(~scenario) + theme_light() + ylim(-0.5,5.5) + xlim(c(3, 9))+
    labs(title="Testing variability, data with noise, blue is true data min-max range:", subtitle=paste0(title))+
   
    geom_segment(data=prop_inside_dat, aes(x=Time, y=max, yend=min), col="blue", lwd=3, alpha=0.4) +
    geom_segment(data=prop_inside_dat,aes(x=Time-0.2, xend=Time+0.2, y=max), col="blue", lwd=1, alpha=0.2) +
    geom_segment(data=prop_inside_dat,aes(x=Time-0.2, xend=Time+0.2, y=min), col="blue", lwd=1, alpha=0.2) +
    geom_segment(data=prop_inside_dat, aes(x=Time, y=q25, yend=q75), col="darkblue", lwd=2, alpha=0.8) +
    geom_segment(data=prop_inside_dat,aes(x=Time-0.1, xend=Time+0.1, y=q25), col="darkblue", lwd=1, alpha=0.6) +
    geom_segment(data=prop_inside_dat,aes(x=Time-0.1, xend=Time+0.1, y=q75), col="darkblue", lwd=1, alpha=0.6)+
    geom_segment(data=dat_quantiles, aes(x=(time+2)/24, y=max, yend=min), col="forestgreen", lwd=3, alpha=0.4) +
    geom_segment(data=dat_quantiles,aes(x=(time+2)/24-0.2, xend=(time+2)/24+0.2, y=max), col="forestgreen", lwd=1, alpha=0.2) +
    geom_segment(data=dat_quantiles,aes(x=(time+2)/24-0.2, xend=(time+2)/24+0.2, y=min), col="forestgreen", lwd=1, alpha=0.2) +
    geom_segment(data=dat_quantiles, aes(x=(time+2)/24, y=q25, yend=q75), col="darkgreen", lwd=2, alpha=0.8) +
    geom_segment(data=dat_quantiles,aes(x=(time+2)/24-0.1, xend=(time+2)/24+0.1, y=q25), col="darkgreen", lwd=1, alpha=0.6) +
    geom_segment(data=dat_quantiles,aes(x=(time+2)/24-0.1, xend=(time+2)/24+0.1, y=q75), col="darkgreen", lwd=1, alpha=0.6)+
    geom_point(data=dat_quantiles , lwd=0.3, col="darkgreen", aes(x=time/24, y=q50), cex=8, alpha=0.7, pch=15) +
    geom_point(data=prop_inside_dat, aes(x=Time, y=q50), col="darkblue",cex=8, alpha=0.7, pch=17))
# ggsave(paste0("varblty_test_060826/", noise,"noise_growth_data_simulated", today(), "_",format(Sys.time(), "%H_%M_%S"), "_numsc_", nrow(setup), ".png"),simu3,
#        width=15, height=10)

(compare <- ggarrange(ggplot(data, aes(x=Time,
                           y=log10.av.ppml, 
                           group=Inoculum.Coh.Subject.ID, 
                           col= Inoculum.Coh.Subject.ID)) + geom_line() + theme_light() +ylim(0,5) + xlim(c(3, 9))+
            theme(legend.position = "none"), simu2, simu3, nrow=3))

ggsave(paste0(folder_name,"/growth_data_simulated_noise_", today(), "_",format(Sys.time(), "%H_%M_%S"), "_numsc_", nrow(setup), ".png"),compare,  width=12, height=18)

invlogit <- function(x){
  return(1/(1+exp(-x)) )
}
