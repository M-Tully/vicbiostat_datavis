# VicBiostat data visualisation workshop 08/10/2026

library(dplyr)
library(ggplot2)
library(gganimate) #info: https://gganimate.com/
library(colorBlindness)#info: https://cran.r-project.org/web/packages/colorBlindness/vignettes/colorBlindness.html

# colour palette: https://coolors.co/
col1 <- "#264653"
col2 <- "#287271"
col3 <- "#2a9d8f"
col4 <- "#8ab17d"
col5 <- "#e9c46a"
col6 <- "#f4a261"
col7 <- "#e76f51"
col8 <- "#be1e2d"  
col9 <- "#9e0059"

# set default ggplot theme
theme_set(theme_light())

# First some simulated parasitaemia data, based on this paper by Wockner et al; https://pubmed.ncbi.nlm.nih.gov/31679015/

setwd("C:/Users/tully/OneDrive - The University of Melbourne/Documents/Presentations/VicBiostat datavis/")

# import data
data <- read.csv( "simu_growth.csv")

# look at the breakdown of unique scenarios/setups
unique(data[,c("pmf", "ipl_mu", "ipl_sd", "scenario")])
 
# ensure these 3 variables are factors
data$ipl_mu_lab <- factor(data$ipl_mu, levels=c(unique(data$ipl_mu)), labels=paste0("mean age = ",unique(data$ipl_mu)))
data$ipl_sd_lab <- factor(data$ipl_sd, levels=c(unique(data$ipl_sd)), labels=paste0("sd age = ",unique(data$ipl_sd)))
data$pmf_lab <- factor(data$pmf, levels=c(unique(data$pmf)), labels=paste0("growth rate = ",unique(data$pmf)))


# hours to days
data$time_d <- data$time/24

# first let's just look at the data, 
ggplot(data, aes(x=time_d, 
                 y=log10(para_cir_noise),
                 group=interaction(id,scenario),
                 col=as.factor(scenario))) + 
  ylim(0,5)+ xlim(3,9)+ 
  geom_line(lwd=0.2) + 
  theme(legend.position = "none") + 
  #  annotation_logticks(sides = "l", outside=T)  +
  #  coord_cartesian(clip = "off")+ # this is req for outside ticks to appear
  labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia")+
  theme(text=element_text(size=20))


# we can facet by scenario
ggplot(data, aes(x=time_d, 
                 y=log10(para_cir_noise),
                 group=id,
                 col=as.factor(id))) + 
  ylim(0,5)+ xlim(3,9)+ 
  geom_line(lwd=0.2) + 
  theme(legend.position = "none") + 
#  annotation_logticks(sides = "l", outside=T)  +
#  coord_cartesian(clip = "off")+ # this is req for outside ticks to appear
  labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia")+
  theme(text=element_text(size=20)) + facet_wrap(~scenario)

# how do we determine the optimal allocation of our 3 variables to aesthetic mappings?

# an obvious starting point is x, y, and colour
(opt1 <- ggplot(data, aes(x=time_d, 
                 y=log10(para_cir_noise),
                 group=interaction(id, pmf_lab),
                 col=pmf_lab)) + 
  ylim(0,5)+ xlim(3,9)+ 
  geom_line(lwd=0.2) + 
  labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia, option 1")+
  theme(text=element_text(size=20)) + facet_grid(ipl_mu_lab~ipl_sd_lab))

(opt2 <- ggplot(data, aes(x=time_d, 
                          y=log10(para_cir_noise),
                          group=interaction(id, ipl_mu_lab),
                          col=ipl_mu_lab)) + 
    ylim(0,5)+ xlim(3,9)+ 
    geom_line(lwd=0.2) + 
    labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia, option 2")+
    theme(text=element_text(size=20)) + facet_grid(pmf_lab~ipl_sd_lab))

(opt3 <- ggplot(data, aes(x=time_d, 
                          y=log10(para_cir_noise),
                          group=interaction(id,ipl_sd_lab),
                          col=ipl_sd_lab)) + 
    ylim(0,5)+ xlim(3,9)+ 
    geom_line(lwd=0.2) + 
    labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia, option 3")+
    theme(text=element_text(size=20)) + facet_grid(ipl_mu_lab~pmf_lab))

ggarrange(opt1+theme(legend.position="bottom"), 
          opt2+theme(legend.position="bottom"), 
          opt3+theme(legend.position="bottom"),  nrow=1)

data_sum <- data %>% group_by(scenario, 
                              pmf_lab,
                              ipl_sd_lab,
                              ipl_mu_lab,
                              time_d) %>% summarise(med = median(para_cir_noise),
                                                              low = quantile(para_cir_noise, 0.025),
                                                              upp = quantile(para_cir_noise, 0.975),
                                                              low2 = quantile(para_cir_noise, 0.25),
                                                              upp2 = quantile(para_cir_noise, 0.75)
                                                              )

(opt1_sum <- ggplot(data_sum, aes(x=time_d, 
                          y=log10(med),
                          group=interaction(pmf_lab),
                          col=pmf_lab)) + 
    ylim(0,5)+ xlim(3,9)+ 
    geom_line(lwd=2) + 
    labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia, option 1")+
    theme(text=element_text(size=20)) + facet_grid(ipl_mu_lab~ipl_sd_lab))

(opt2_sum <- ggplot(data_sum, aes(x=time_d, 
                          y=log10(med),
                          group=interaction(ipl_mu_lab),
                          col=ipl_mu_lab)) + 
    ylim(0,5)+ xlim(3,9)+ 
    geom_line(lwd=2) + 
    labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia, option 2")+
    theme(text=element_text(size=20)) + facet_grid(pmf_lab~ipl_sd_lab))

(opt3_sum <- ggplot(data_sum, aes(x=time_d, 
                          y=log10(med),
                          group=interaction(ipl_sd_lab),
                          col=ipl_sd_lab)) + 
    ylim(0,5)+ xlim(3,9)+ 
    geom_line(lwd=2) + 
    labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia, option 3")+
    theme(text=element_text(size=20)) + facet_grid(ipl_mu_lab~pmf_lab))

ggarrange(opt1_sum+theme(legend.position="bottom"), 
          opt2_sum+theme(legend.position="bottom"), 
          opt3_sum+theme(legend.position="bottom"),  nrow=1)


# how about with range ribbons?
(opt1_rib <- ggplot(data_sum, aes(x=time_d, 
                                  y=log10(med),
                                  group=interaction(pmf_lab),
                                  col=pmf_lab)) + 
    coord_cartesian(ylim=c(0,5), xlim=c(3,9))+ 
    geom_line(lwd=2) + 
    geom_ribbon(aes(x=time_d, 
                    ymax = log10(upp), 
                    ymin=log10(low),
                    group=pmf_lab, 
                    fill=pmf_lab), alpha=0.1) +
    labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia, option 1")+
    theme(text=element_text(size=20)) + facet_grid(ipl_mu_lab~ipl_sd_lab))

(opt2_rib <- ggplot(data_sum, aes(x=time_d, 
                                  y=log10(med),
                                  group=interaction(ipl_mu_lab),
                                  col=ipl_mu_lab)) + 
    coord_cartesian(ylim=c(0,5), xlim=c(3,9))+ 
    geom_line(lwd=2) + 
    geom_ribbon(aes(x=time_d, 
                    ymax = log10(upp), 
                    ymin=log10(low),
                    group=ipl_mu_lab, 
                    fill=ipl_mu_lab), alpha=0.1) +
    labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia, option 2")+
    theme(text=element_text(size=20)) + facet_grid(pmf_lab~ipl_sd_lab))

(opt3_rib <- ggplot(data_sum, aes(x=time_d, 
                                  y=log10(med),
                                  group=interaction(ipl_sd_lab),
                                  col=ipl_sd_lab)) + 
    coord_cartesian(ylim=c(0,5), xlim=c(3,9))+ 
    geom_line(lwd=2) + 
    geom_ribbon(aes(x=time_d, 
                    ymax = log10(upp), 
                    ymin=log10(low),
                    group=ipl_sd_lab, 
                    fill=ipl_sd_lab), alpha=0.1) +
    labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia, option 3")+
    theme(text=element_text(size=20)) + facet_grid(ipl_mu_lab~pmf_lab))

ggarrange(opt1_rib+theme(legend.position="bottom"), 
          opt2_rib+theme(legend.position="bottom"), 
          opt3_rib+theme(legend.position="bottom"),  nrow=1)



# double ribbons?
(opt3_sum <- ggplot(data_sum, aes(x=time_d, 
                                  y=log10(med),
                                  col=ipl_sd_lab)) + 
    coord_cartesian(ylim=c(0,5), xlim=c(3,9))+ 
    
    geom_ribbon(aes(x=time_d, 
                    ymax = log10(upp), 
                    ymin=log10(low),
                    group=ipl_sd_lab,
                    col=NULL,
                    fill=ipl_sd_lab), alpha=0.1) +
    geom_ribbon(aes(x=time_d, 
                    ymax = log10(upp2), 
                    ymin=log10(low2),
                    group=ipl_sd_lab, 
                    col=NULL,
                    fill=ipl_sd_lab), alpha=0.4) +
    geom_line(lwd=0.4) + 
    labs(x="Days since inoculation", y="log10(Parasites/mL)", title="Simulated Parasitaemia, option 3")+
    theme(text=element_text(size=20)) + facet_grid(ipl_mu_lab~pmf_lab))

#tidy final version
data_sum$ipl_sd_lab2 <- factor(data_sum$ipl_sd_lab, levels = unique(data_sum$ipl_sd_lab),
                              labels = c(expression(sigma[ipl]==1),
                                         expression(sigma[ipl]==5) ))

data_sum$pmf_lab2 <- factor(data_sum$pmf_lab, levels = unique(data_sum$pmf_lab),
                              labels = c(expression(PMF==10),
                                         expression(PMF==20)))

# lets clean the favourite up a bit!
(opt2_final <- ggplot(data_sum, aes(x=time_d, 
                                  y=log10(med),
                                  group=interaction(ipl_mu_lab),
                                  col=ipl_mu_lab)) + 
    annotation_logticks(sides = "l", outside=F)  +
    coord_cartesian(xlim=c(3,9), ylim=c(0,5))+ 
    scale_y_continuous(breaks=c(0, 2, 4), 
                       labels=c(expression(1),
                                expression(100),
                                expression(10^4))) +
    geom_line(lwd=1.4) + 
    geom_ribbon(aes(x=time_d, 
                    ymax =log10(upp), 
                    ymin=log10(low),
                    group=ipl_mu_lab, 
                    fill=ipl_mu_lab), alpha=0.1) +
    labs(x="Days since inoculation", y="parasites/mL", 
         title="Simulated Parasitaemia, clean version", 
         fill=expression(mu[ipl]), 
         col=expression(mu[ipl]))+
    scale_color_manual(values=c(col4, col2, col8))+
    scale_fill_manual(values=c(col4, col2, col8  ))+
    theme(text=element_text(size=20), strip.text.y = element_text(family="serif")) + 
    facet_grid(ipl_sd_lab2~pmf_lab2, labeller=label_parsed))

# check colour blindness
cvdPlot(opt2_final)


data$ipl_sd_lab2 <- factor(data$ipl_sd_lab, levels = unique(data$ipl_sd_lab),
                               labels = c(expression(sigma[ipl]==1),
                                          expression(sigma[ipl]==5) ))

data$pmf_lab2 <- factor(data$pmf_lab, levels = unique(data$pmf_lab),
                            labels = c(expression(PMF==10),
                                       expression(PMF==20)))

# how about stat_smooth instead?
(opt2_final2 <- ggplot() + 
    annotation_logticks(sides = "l", outside=F)  +
    coord_cartesian(xlim=c(3,9), ylim=c(0,5))+ 
    scale_y_continuous(breaks=c(0, 2, 4), 
                       labels=c(expression(1),
                                expression(100),
                                expression(10^4))) +
    geom_ribbon(data=data_sum, aes(x=time_d, 
                    ymax =log10(upp), 
                    ymin=log10(low),
                    col=ipl_mu_lab, 
                    linetype=ipl_mu_lab,
                    fill=ipl_mu_lab), alpha=0.1) +
    stat_smooth(data=data, aes(x=time_d, y=log10(para_cir_noise), 
                               col=ipl_mu_lab, group=ipl_mu_lab),method="loess", span=0.2)+
    labs(x="Days since inoculation", y="parasites/mL", 
         title="Simulated Parasitaemia, clean version, loess curve", 
         fill=expression(mu[ipl]), 
         linetype=expression(mu[ipl]), 
         col=expression(mu[ipl]))+
    scale_color_manual(values=c(col4, col2, col8))+
    scale_fill_manual(values=c(col4, col2, col8  ))+
    theme(text=element_text(size=20), strip.text.y = element_text(family="serif")) + 
    facet_grid(ipl_sd_lab2~pmf_lab2, labeller=label_parsed))




# animate via reveal 
opt2_final + geom_point(cex=3) + transition_reveal(time_d)

anim <- opt2_final + geom_point(cex=3) + transition_reveal(time_d)
anim_gif <- animate(anim, height = 700, width =800, renderer =gifski_renderer(), fps=7, nframes=200)
anim_save("example.gif", animation = anim_gif, 
          path="C:/Users/tully/OneDrive - The University of Melbourne/Documents/Presentations/VicBiostat datavis/")

# online gif maker: https://ezgif.com/maker

