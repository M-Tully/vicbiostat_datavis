# VicBiostat data visualisation workshop 08/10/2026

library(dplyr)
library(ggplot2)
library(gganimate) #info: https://gganimate.com/
library(colorBlindness)#info: https://cran.r-project.org/web/packages/colorBlindness/vignettes/colorBlindness.html
require(patchwork)
require(cowplot)

# PK data made available by: Abd-Rahman, A.N., Kaschek, D., Kümmel, A. et al. 
# Characterizing the pharmacological interaction of the antimalarial combination 
# artefenomel-piperaquine in healthy volunteers with induced blood-stage Plasmodium 
# falciparum to predict efficacy in patients with malaria. BMC Med 22, 563 (2024). 
# https://doi.org/10.1186/s12916-024-03787-0

# The actual study that originally collected this data is:
# McCarthy JS, Baker M, O'Rourke P, Marquart L, Griffin P, Hooft van Huijsduijnen R, Möhrle JJ. 
# Efficacy of OZ439 (artefenomel) against early Plasmodium falciparum blood-stage malaria infection 
# in healthy volunteers. J Antimicrob Chemother. 2016 Sep;71(9):2620-7. doi: 10.1093/jac/dkw174.
# Epub 2016 Jun 5. PMID: 27272721; PMCID: PMC4992851.
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4992851/

setwd("C:/Users/tully/OneDrive - The University of Melbourne/Documents/Presentations/VicBiostat datavis/artefenomel-ppq")
pk <- read.csv("12916_2024_3787_MOESM5_ESM (2).csv") # VIS PK
pk <- pk[pk$STUDY=="QP12C10"&pk$NAME=="OZ439:::Plasma:::Concentration",] # we only want the mono-therapy Artefenomel data

pk$conc <- pk$VALUE*1000 #converting from ug/mL to ng/mL

# to summarise by time-point (not exact time of sample) we need to round 'TIME'
# some will be to nearest 12h, but the first ~50 hours are richer sampling
#cutoff <- 50
#pk$time[pk$TIME<cutoff] <- round(pk$TIME[pk$TIME<cutoff])
#pk$time[pk$TIME>cutoff] <- 12*round(pk$TIME[pk$TIME>cutoff]/12)
pk$time <- pk$TIME

#colour by dose level; log10 scale
ggplot(pk, aes(x=time,y=(conc), col=DOSELEVEL1, group=ID)) + 
  geom_line()+
  scale_y_log10(guide="axis_logticks") + 
  labs(x="Time (h)", y="Artefenomel plasma concentation (ng/mL)") + 
  cowplot::theme_half_open()+
  theme(text=element_text(size=15), 
        legend.title = element_blank(), 
        legend.key.size = unit(2, "lines"))

pk$dose <- as.factor(pk$DOSELEVEL1)
# tidy up factor, perhaps look at early samples more closely
ggplot(pk, aes(x=time,y=(conc), col=dose, group=ID)) + 
  geom_line(lwd=0.8, alpha=0.3)+
  scale_y_log10(guide="axis_logticks") + 
  labs(x="Time (h)", y="Artefenomel plasma concentation (ng/mL)") + 
  cowplot::theme_half_open()+
  theme(text=element_text(size=15), 
        legend.title = element_blank(), 
        legend.key.size = unit(2, "lines")) #+xlim(0,100)

# facet wrap
ggplot(pk, aes(x=time,y=(conc), col=dose, group=ID)) + 
  geom_line(lwd=0.1, alpha=0.5)+
  scale_y_log10(guide="axis_logticks") + 
  labs(x="Time (h)", y="Artefenomel plasma concentation (ng/mL)") + 
  cowplot::theme_half_open()+
  theme(text=element_text(size=15), 
        legend.title = element_blank(), 
        legend.key.size = unit(2, "lines")) +facet_wrap(~dose, nrow=3)

# facet wrap + grey example
ggplot() + 
  geom_line(data=pk[,c("time", "conc", "ID")], 
            aes(x=time,y=(conc), group=ID), col="grey",
            lwd=0.1, alpha=0.5)+
  geom_line(data=pk, 
            aes(x=time,y=(conc), col=dose, group=ID),
            lwd=0.1, alpha=0.5)+
  scale_y_log10(guide="axis_logticks") + 
  labs(x="Time (h)", y="Artefenomel plasma concentation (ng/mL)") + 
  cowplot::theme_half_open()+
  theme(text=element_text(size=15), 
        legend.title = element_blank(), 
        legend.key.size = unit(2, "lines")) +facet_wrap(~dose, nrow=3)

# add in stat_smooth
ggplot() + 
  geom_line(data=pk[,c("time", "conc", "ID")], 
            aes(x=time,y=(conc), group=ID), col="grey",
            lwd=0.1, alpha=0.5)+
  geom_line(data=pk, 
            aes(x=time,y=(conc), col=dose, group=ID),
            lwd=0.1, alpha=0.5)+
  stat_smooth(data= pk#[pk$time>5,]
              , aes(x=time,y=(conc), col=dose),se=F#, span=0.1
              )+
  scale_y_log10(guide="axis_logticks") + 
  labs(x="Time (h)", y="Artefenomel plasma concentation (ng/mL)") + 
  cowplot::theme_half_open()+
  theme(text=element_text(size=15), 
        legend.title = element_blank(), 
        legend.key.size = unit(2, "lines")) +facet_wrap(~dose, nrow=3)

# only stat_smooth
ggplot(pk#[pk$time>5,]
       ) + 
  stat_smooth( aes(x=time, y=conc, col=dose, group=dose), se=F, span=0.5)+
  scale_y_log10(guide="axis_logticks") + 
  labs(x="Time (h)", y="Artefenomel plasma concentation (ng/mL)") + 
  cowplot::theme_half_open()+
  theme(text=element_text(size=15), 
        legend.title = element_blank(), 
        legend.key.size = unit(2, "lines"))


# calculate medians for these new groups observation times;
pk_sum <- pk %>% group_by(DOSELEVEL1, time) %>% filter(time>0) %>% mutate(mean_conc = mean(conc),
                                                       med_conc = median(conc),
                                                       low_conc = quantile(conc, 0.025),
                                                       upp_conc = quantile(conc, 0.975),
                                                       dose = factor(DOSELEVEL1, levels=c(100,200,500), 
                                                                     labels=c("100 mg", "200 mg", "500 mg"))) %>%
  dplyr::select(mean_conc, med_conc,upp_conc, low_conc, dose) %>% distinct()


ggplot(pk_sum, aes(x=time, group=dose, y=(med_conc), col=dose)) + 
  geom_line(lwd=1)+
  geom_point(aes(pch=dose), cex=4)+
  scale_y_log10(guide="axis_logticks") + 
  scale_colour_manual(values=c("black",  "grey80","grey50")) +
  scale_shape_manual(values=c(16,15,17)) +
  labs(x="Time (h)", y="Artefenomel plasma concentation (ng/mL)") + 
  xlim(0,500) +
  cowplot::theme_half_open()+
  theme(text=element_text(size=15), 
        legend.title = element_blank(), 
        legend.key.size = unit(2, "lines"))

#issue with this one;
# need to group time to the neares observation point
cutoff <- 46
pk$time[pk$TIME<cutoff] <- round(pk$TIME[pk$TIME<cutoff])
pk$time[pk$TIME>cutoff] <- 12*round(pk$TIME[pk$TIME>cutoff]/12)

# manually change last point because the tmax is 312 instad of 450 on the plot
pk_sum$time[pk_sum$time>300] <- 450
# also remove the one at 288 because they seem to have censored it?
pk_sum <- pk_sum %>% filter(time != 288)

# repeat, (now with fake last point)
(pub_fig <- ggplot(pk_sum, aes(x=time, group=dose, y=(med_conc), col=dose)) + 
  geom_line(lwd=1)+
  geom_point(aes(pch=dose), cex=4)+
  scale_y_log10(guide="axis_logticks") + 
  scale_colour_manual(values=c("black",  "grey80","grey50")) +
  scale_shape_manual(values=c(16,15,17)) +
  labs(x="Time (h)", y="Artefenomel plasma concentation\n(ng/mL)") + 
  scale_x_continuous(limits = c(0, 500), breaks = scales::breaks_width(100)) +
  cowplot::theme_half_open()+
  theme(text=element_text(size=15), 
        legend.title = element_blank(), 
        legend.key.size = unit(2, "lines")))#+ transition_reveal(time)

p <- pub_fig

p | p

p / p

p_combined <- 
  ((p + xlim(c(0, 200))) | p ) + 
  patchwork::plot_layout(
    guides = "collect", widths = c(4, 1))

ggsave(
  plot = p_combined,
  device = svglite::svglite,
  here::here("Part 3 - Inkscape demo/plot_out.svg"),
  width = 18/1.5, height = 6/1.5)

ggsave(
  plot = p_combined,
  here::here("Part 3 - Inkscape demo/plot_out.pdf"),
  width = 18/1.5, height = 6/1.5)
#-----------------------------------------------------------



pk$dose <- factor(pk$DOSELEVEL1, levels=c(100,200,500), 
                  labels=c("100 mg", "200 mg", "500 mg"))

# repeat, with background
ggplot(pk_sum, aes(x=time, group=dose, y=(med_conc), col=dose)) + 
  geom_line(lwd=1)+
  geom_line(data=pk, 
            aes(x=time,y=(conc), col=dose, group=ID),
            lwd=0.1, alpha=0.5)+
  geom_point(aes(pch=dose), cex=4)+
  scale_y_log10(guide="axis_logticks") + 
# scale_colour_manual(values=c("black",  "grey80","grey50")) +
  scale_shape_manual(values=c(16,15,17)) +
  labs(x="Time (h)", y="Artefenomel plasma concentation (ng/mL)") + 
  xlim(0,500) +
  cowplot::theme_half_open()+
  theme(text=element_text(size=15), 
        legend.title = element_blank(), 
        legend.key.size = unit(2, "lines"))

# ribbons?
ggplot(pk_sum, aes(x=time, y=(med_conc), col=dose)) + 
  geom_line(lwd=1)+
  geom_ribbon(data=pk_sum, aes(x=time, ymin=low_conc, ymax=upp_conc, fill=dose), alpha=0.2)+
  geom_point(aes(pch=dose), cex=3)+
  scale_y_log10() + 
  coord_cartesian(ylim=c(1,10^3))+
 # scale_colour_manual(values=c("black",  "grey80","grey50")) +
  scale_shape_manual(values=c(16,15,17)) +
  labs(x="Time (h)", y="Artefenomel plasma concentation (ng/mL)") + 
  xlim(0,500) +
  cowplot::theme_half_open()+
  theme(text=element_text(size=15), 
        legend.title = element_blank(), 
        legend.key.size = unit(2, "lines"))#+ transition_reveal(time)



# animate via reveal 
pub_fig + transition_reveal(time)

anim <- opt2_final + geom_point(cex=3) + transition_reveal(time_d)
anim_gif <- animate(anim, height = 700, width =800, renderer =gifski_renderer(), fps=7, nframes=200)
anim_save("example.gif", animation = anim_gif, 
          path="C:/Users/tully/OneDrive - The University of Melbourne/Documents/Presentations/VicBiostat datavis/")



fig_density <- ggplot(pk[pk$time>0&!(pk$time %in%c(288, 312)),]) + 
  geom_density(aes(x=log10(conc), fill=dose),
           # lwd=0.1, 
           alpha=0.5)+
  labs(title='Hours from dose:{closest_state}', 
       x="Artefenomel plasma concentation (ng/mL)") + 
  cowplot::theme_half_open()+
  scale_fill_manual(values=c(col1,  col3,col7)) +
  theme(text=element_text(size=20), 
        legend.title = element_blank(), 
        legend.key.size = unit(2, "lines"),
        ) +transition_states(time, state_length=3)
anim_gif <- animate(fig_density, height = 700, width =800, renderer =gifski_renderer(), fps=8, nframes=200)
anim_save("pk_density.gif", animation = anim_gif, 
          path="C:/Users/tully/OneDrive - The University of Melbourne/Documents/Presentations/VicBiostat datavis/")

fig_median <- ggplot() + 
  geom_line(data=pk_sum[pk_sum$time>0,], aes(x=time, group=dose, y=(med_conc), col=dose), lwd=1)+
  geom_point(data=pk_sum[pk_sum$time>0,], aes(x=time, group=dose, y=(med_conc), col=dose, pch=dose), cex=4)+
  scale_y_log10(guide="axis_logticks") + 
  geom_vline(data=data.frame(x=unique(pk_sum$time[pk_sum$time>=0&pk_sum$time<150])), aes(xintercept=x)) +
  scale_colour_manual(values=c(col1,  col3,col7)) +
  scale_shape_manual(values=c(16,15,17)) +
  labs(title='Hours from dose:{closest_state}', 
       x="Time (h)", y="Artefenomel plasma concentation (ng/mL)") + 
  xlim(0,150) +
  cowplot::theme_half_open()+
  theme(text=element_text(size=20), 
        legend.title = element_blank(), 
        legend.key.size = unit(2, "lines"))+transition_states(x)
anim_gif2 <- animate(fig_median, height = 700, width =800, renderer =gifski_renderer(), fps=8, nframes=200)
anim_save("pk_medians.gif", animation = anim_gif2, 
          path="C:/Users/tully/OneDrive - The University of Melbourne/Documents/Presentations/VicBiostat datavis/")




