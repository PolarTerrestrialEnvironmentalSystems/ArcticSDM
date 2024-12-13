

library(tidyverse)

library(sf)

sf_use_s2(FALSE)

library(stars)

library(wesanderson)

library(ggtext)

library(ggpubr)


setwd("//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/")
load("sdm_occurance_array_metadata_25km.rda")

species <- read_csv("species_sigma1.csv")


load("sig1Array.rda")

sig1 <- spArray

load("sig5Array.rda")

sig5 <- spArray


setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/hypotheses_results/plot1 data")
load("NWdist.rda")
load("New.rda")
load("NewSm.rda")

setwd("C:/Users/roschw001/Documents/R/SDM/SDM/ArcticSDMserver/ArcticSDM/hypotheses_results/plot1 data")

newSm <- read.csv("newSm_distcons.csv")
NWsum <- read.csv("NWsum_distcons.csv")

### Emerging habitat ###

########################


# New <- lapply(1:dim(sig1)[4], function(s) {
#   
#   parallel::mclapply(1:dim(sig1)[1], function(i) {
#     
#     # cumsum(apply(t(apply(sig1[i,,,s], 1, diff)), 2, function(z) sum(z>0)))
#     
#     apply(t(apply(sig1[i,,,s], 1, function(x) x[-1] > x[1])), 2, sum)
#     
#   }, mc.cores = 1) %>% Reduce("rbind",.) %>% as_tibble() %>%
#     
#     setNames(c("2026", "2056", "2086")) %>%
#     
#     pivot_longer(cols = names(.), values_to = 'cells') %>%
#     
#     mutate(year = as.numeric(name),
#            
#            scenario = s) %>% dplyr::select(scenario, year, cells)
#   
# }) %>% Reduce("rbind",.) %>% mutate(disp = "sig1") %>%
#   
#   bind_rows(
#     
#     lapply(1:dim(sig5)[4], function(s) {
#       
#       parallel::mclapply(1:dim(sig5)[1], function(i) {
#         
#         # cumsum(apply(t(apply(sig5[i,,,s], 1, diff)), 2, function(z) sum(z>0)))
#         
#         apply(t(apply(sig5[i,,,s], 1, function(x) x[-1] > x[1])), 2, sum)
#         
#       }, mc.cores = 1) %>% Reduce("rbind",.) %>% as_tibble() %>%
#         
#         setNames(c("2026", "2056", "2086")) %>%
#         
#         pivot_longer(cols = names(.), values_to = 'cells') %>%
#         
#         mutate(year = as.numeric(name),
#                
#                scenario = s) %>% dplyr::select(scenario, year, cells)
#       
#     }) %>% Reduce("rbind",.) %>% mutate(disp = "sig5")) %>% mutate(area = cells * 25 *25)
# 
# 
# 
# 
# 
# newSm <- New %>% filter(area>0) %>% group_by(scenario, year, disp) %>%
#   
#   summarise(median = median(area, na.rm = T),
#             
#             mean = mean(area, na.rm = T),
#             
#             sd = sd(area, na.rm = T),
#             
#             lower_50 = quantile(area, probs = 0.25, na.rm = T),
#             
#             lower_95 = quantile(area, probs = 0.05, na.rm = T),
#             
#             upper_50 = quantile(area, probs = 0.75, na.rm = T),
#             
#             upper_95 = quantile(area, probs = 0.975, na.rm = T)
#             
#   ) %>% mutate(x_axis = year + (as.numeric(as.character(scenario))-2)*5)
# 

#####################
#New--> newSm


#### newSm ####


newSm <- newSm %>%
  mutate(x_axis_old = x_axis)


newSm <- newSm %>%
  mutate(yearnew = case_when(
  year == 2021 ~ 2040,
  year == 2051 ~ 2070,
  year == 2086 ~ 2100
  
))

pl1 <- ggplot(newSm %>% mutate(disp = factor(disp, levels=c("sig5","sig1"))), 
              
              aes(x = x_axis, y = median, color = as.factor(scenario), fill = as.factor(scenario))) +
  
  geom_bar(stat="identity", color = NA) +
  
  geom_segment(mapping = aes(x = x_axis, xend = x_axis,
                             
                             y = lower_50, yend = upper_50), color = "grey40") +
  
 # scale_fill_manual(values = wes_palette("Darjeeling1")[2:4], labels = c("SSP 1.26", "SSP 3.70", "SSP 5.85")) +
  scale_fill_manual(values = hcl.colors(6, "heat")[3:1], labels = c("SSP 1.26", "SSP 3.70", "SSP 5.85")) +
 
  # Customizing x-axis breaks and labels
  scale_x_continuous(breaks = c(2026, 2056, 2086), 
                     labels = c("2040", "2070", "2100")) + 
 #  
 # scale_x_continuous(breaks = c(2026, 2056, 2086)) +
 #  
 #  scale_x_continuous(breaks = c(2040, 2070, 2100)) +
  
  facet_wrap(~disp, labeller = labeller(disp = c("sig5" = "unconstrained","sig1" = "dispersal constrained"))) +
  
  guides(fill = guide_legend(title = "Scenarios")) +
  
  labs(x = "", y = "Emerging habitat [\u0394 1980, km <sup>2</sup>]") +
  
  theme_light() +
  
  theme(strip.text = element_text(size = 15),
        
        axis.text = element_text(size = 12),
        
        axis.title = element_text(size = 16),
        
        axis.title.y = element_markdown(),
        
        legend.position = c(0.085, 0.75),#0.095
        
        legend.key.size = unit(0.6, 'cm'),
        
        legend.text = element_text(size=11),
        
        legend.title = element_text(size=15),
        
        strip.background = element_rect(fill = "white", color = "black"),
        
        strip.text.x = element_text(color = "black"),
        
        panel.border = element_rect(colour = "black", fill=NA))


print(pl1)



##### Ndist ####

ggplot(NWdist, aes(x = as.numeric(name), y = dist, group = scenario, color = as.factor(scenario))) +
  
  geom_smooth(method = "lm") +
  
  scale_color_manual(values = wes_palette("Darjeeling1")[2:4], labels = c("SSP 1.26", "SSP 3.70", "SSP 5.85")) +
  
  facet_wrap(~disp)





NWsum <- NWdist %>% group_by(scenario, year, disp) %>%
  
  summarise(median = median(dist, na.rm = T),
            
            mean = mean(dist, na.rm = T),
            
            sd = sd(dist, na.rm = T),
            
            lower_50 = quantile(dist, probs = 0.25, na.rm = T),
            
            lower_95 = quantile(dist, probs = 0.05, na.rm = T),
            
            upper_50 = quantile(dist, probs = 0.75, na.rm = T),
            
            upper_95 = quantile(dist, probs = 0.975, na.rm = T)
            
  ) %>% mutate(x_axis = year + (as.numeric(as.character(scenario))-2)*5)

#NWsum$year
NWsum <- NWsum %>%
  mutate(yearnew = case_when(
    year == 2026 ~ 2040,
    year == 2056 ~ 2070,
    year == 2086 ~ 2100
    
  ))

#### p2 ####
# pl2 <- ggplot(NWsum %>% mutate(disp = factor(disp, levels=c("sig5","sig1"))), 
#               
#               aes(x = x_axis, y = median, color = as.factor(scenario), fill = as.factor(scenario))) +
#   # 
#   # geom_bar(stat="identity", color = NA) +
#   # 
#   # geom_segment(mapping = aes(x = x_axis, xend = x_axis,
#   #                            
#   #                            y = lower_50, yend = upper_50), color = "grey40") +
#   
#   geom_point(size = 5, show.legend = F, alpha = 0.5) +
#   
#   geom_segment(mapping = aes(x = x_axis, xend = x_axis,
#                              
#                              y = lower_50, yend = upper_50), show.legend = F, alpha = 0.5) +
#   
#   # scale_fill_manual(values = wes_palette("Darjeeling1")[2:4], labels = c("SSP 1.26", "SSP 3.70", "SSP 5.85")) +
#   scale_color_manual(values = hcl.colors(6, "heat")[3:1], labels = c("SSP 1.26", "SSP 3.70", "SSP 5.85")) +
#   
#   
#   scale_x_continuous(breaks = c(2026, 2056, 2086)) +
#   
#   facet_wrap(~disp, labeller = labeller(disp = c("sig5" = "unconstrained","sig1" = "dispersal constrained"))) +
#   
#   guides(fill = guide_legend(title = "Scenarios")) +
#   
#   labs(x = "", y = "Distance to north [km]") +
#   
#   theme_light() +
#   
#   theme(strip.text = element_text(size = 15),
#         
#         axis.text = element_text(size = 12),
#         
#         axis.title = element_text(size = 16),
#         
#         axis.title.y = element_markdown(),
#         
#         legend.position = c(0.085, 0.75),#0.095
#         
#         legend.key.size = unit(0.6, 'cm'),
#         
#         legend.text = element_text(size=11),
#         
#         legend.title = element_text(size=15),
#         
#         strip.background = element_rect(fill = "white", color = "black"),
#         
#         strip.text.x = element_text(color = "black"),
#         
#         panel.border = element_rect(colour = "black", fill=NA))
# 
# 
# print(pl2)


### Northward shift ###

#######################


# NW <- lapply(1:dim(sig1)[4], function(s) {
#   
#   parallel::mclapply(1:dim(sig1)[1], function(i) {
#     
#     t(apply(sig1[i,,,s], 1, diff)) %>% as_tibble() %>% setNames(c("2026", "2056", "2086")) %>%
#       
#       mutate(id = 1:nrow(.)) %>% pivot_longer(cols = -id) %>%
#       
#       filter(value>0)
#     
#   }, mc.cores = 1) %>% Reduce("rbind", .) %>% mutate(scenario = s)
#   
# }) %>% Reduce("rbind", .) %>% mutate(disp = "sig1") %>% bind_rows(
#   
#   lapply(1:dim(sig5)[4], function(s) {
#     
#     parallel::mclapply(1:dim(sig5)[1], function(i) {
#       
#       t(apply(sig5[i,,,s], 1, diff)) %>% as_tibble() %>% setNames(c("2026", "2056", "2086")) %>%
#         
#         mutate(id = 1:nrow(.)) %>% pivot_longer(cols = -id) %>%
#         
#         filter(value>0)
#       
#     }, mc.cores = 1) %>% Reduce("rbind", .) %>% mutate(scenario = s)
#     
#   }) %>% Reduce("rbind", .) %>% mutate(disp = "sig5")
#   
# )
# 
# 
# distGrid <- spArrayMeta$dim2 %>% 
#   
#   mutate(id = 1:nrow(.), dist = as.numeric(st_distance(., st_point(c(0,0)) %>% st_sfc(crs = st_crs(spArrayMeta$dim2))))/1000) %>%
#   
#   dplyr::select(id, dist) %>% st_drop_geometry()
# 
# 
# NWdist <- NW %>% left_join(distGrid, by = "id") %>% mutate(year = as.numeric(name))


pl2 <- ggplot(NWsum %>% mutate(disp = factor(disp, levels=c("sig5","sig1"))),
              
              aes(x = x_axis, y = median, color = as.factor(scenario), fill = as.factor(scenario))) +
  
  geom_point(size = 5, show.legend = F, alpha = 1) +
  
  geom_segment(mapping = aes(x = x_axis, xend = x_axis,
                             
                             y = lower_50, yend = upper_50), show.legend = F, alpha = 1) +
  
  # geom_smooth(data = NWdist %>% mutate(x_axis = year + (as.numeric(as.character(scenario))-2)*5), 
  # 
  #            mapping = aes(x = x_axis, y = dist, group = scenario, color = as.factor(scenario)), method = "lm",
  # 
  #             show.legend = F) +
  # 
 # scale_color_manual(values = wes_palette("Darjeeling1")[2:4], labels = c("SSP 1.26", "SSP 3.70", "SSP 5.85")) +
  scale_color_manual(values = hcl.colors(6, "heat")[3:1], labels = c("SSP 1.26", "SSP 3.70", "SSP 5.85")) +
  
  
  # Customizing x-axis breaks and labels
  scale_x_continuous(breaks = c(2026, 2056, 2086), 
                     labels = c("2040", "2070", "2100")) + 
  
  #scale_x_continuous(breaks = c(2026, 2056, 2086)) +
  
  facet_wrap(~disp, labeller = labeller(disp = c("sig5" = "","sig1" = ""))) +
  
  guides(fill = guide_legend(title = "Scenarios")) +
  
  xlab("") + ylab(bquote("Distance to north [km]")) +
  theme_light() +
  
  theme(strip.text = element_text(size = 15),
        
        axis.text = element_text(size = 12),
        
        axis.title = element_text(size = 16),
      
        legend.position = c(0.085, 0.75),#0.095
        
        legend.key.size = unit(0.6, 'cm'),
        
        legend.text = element_text(size=11),
        
        legend.title = element_text(size=15),
        
        strip.background = element_blank(),
        
        strip.text.x = element_blank(),
        
        panel.border = element_rect(colour = "black", fill=NA)
        
  )

#### p2-median ####

m = NWsum$median[1]
#m = 3042774  
# 
# pl2 <- ggplot(NWsum %>% mutate(disp = factor(disp, levels=c("sig5","sig1"))),
#               
#               aes(x = x_axis, y = median-m, color = as.factor(scenario), fill = as.factor(scenario))) +
#   
#   geom_point(size = 5, show.legend = F, alpha = 0.5) +
#   
#   geom_segment(mapping = aes(x = x_axis, xend = x_axis,
#                              
#                              y = lower_50 - m, yend = upper_50 - m), show.legend = F, alpha = 0.5) +
#   
#   
  
#pl2 <- ggplot(NWsum %>% mutate(disp = factor(disp, levels=c("sig5","sig1"))),

sig5 <- subset(NWsum, NWsum$disp == "sig5")
              
p2 <- ggplot(sig5 %>% mutate(disp = factor(disp, levels= "sig5")),
              
              aes(x = x_axis, y = median, color = as.factor(scenario), fill = as.factor(scenario))) +
  
  geom_point(size = 5, show.legend = F, alpha = 1) +
  
  geom_segment(mapping = aes(x = x_axis, xend = x_axis,
                             
                             y = lower_50, yend = upper_50), show.legend = F, alpha = 0.5) +
  
  # geom_smooth(data = NWdist %>% mutate(x_axis = year + (as.numeric(as.character(scenario))-2)*5), 
  # 
  #            mapping = aes(x = x_axis, y = dist, group = scenario, color = as.factor(scenario)), method = "lm",
  # 
  #             show.legend = F) +
  # 
  # scale_color_manual(values = wes_palette("Darjeeling1")[2:4], labels = c("SSP 1.26", "SSP 3.70", "SSP 5.85")) +
  scale_color_manual(values = hcl.colors(6, "heat")[3:1], labels = c("SSP 1.26", "SSP 3.70", "SSP 5.85")) +
  
  
  # Customizing x-axis breaks and labels
  scale_x_continuous(breaks = c(2026, 2056, 2086), 
                     labels = c("2040", "2070", "2100")) + 
  
  #scale_x_continuous(breaks = c(2026, 2056, 2086)) +
  
  #facet_wrap(~disp, labeller = labeller(disp = c("sig5" = "","sig1" = ""))) +
  facet_wrap(~disp, labeller = labeller(disp = c("sig5" = ""))) +
  
  guides(fill = guide_legend(title = "Scenarios")) +
  
  xlab("") + ylab(bquote("Distance to north [km]")) +
  theme_light() +
  
  theme(strip.text = element_text(size = 15),
        
        axis.text = element_text(size = 12),
        
        axis.title = element_text(size = 16),
        
        legend.position = c(0.085, 0.75),#0.095
        
        legend.key.size = unit(0.6, 'cm'),
        
        legend.text = element_text(size=11),
        
        legend.title = element_text(size=15),
        
        strip.background = element_blank(),
        
        strip.text.x = element_blank(),
        
        panel.border = element_rect(colour = "black", fill=NA)
        )

print(p2)

pl3 <- p2 + geom_line(data = sig5, 
          aes(x = x_axis, y = median, color = as.factor(scenario)), 
          size = 0.2, linetype = "dashed", show.legend = FALSE)  # Line color according to scen

print(pl3)
ggarrange(pl1, pl3, nrow = 2)




# pl2 <- ggplot(newSm %>% mutate(disp = factor(disp, levels = c("sig5", "sig1"))), 
#               aes(x = x_axis, y = median, color = as.factor(scenario), fill = as.factor(scenario))) +
#   
#   
#   
#   geom_point(size = 5, show.legend = F, alpha = 0.5) +
#   
#   geom_segment(mapping = aes(x = x_axis, xend = x_axis,
#                              y = lower_50, yend = upper_50), alpha = 0.5) +
#   
#   scale_color_manual(values = hcl.colors(6, "heat")[3:1], labels = c("SSP 126", "SSP 370", "SSP 585")) +
#   
#   scale_fill_manual(values = hcl.colors(6, "heat")[3:1], labels = c("SSP 126", "SSP 370", "SSP 585"), guide = guide_legend(title = "Scenarios")) + # Ensure fill has a legend
#   
#   scale_x_continuous(breaks = c(2026, 2056, 2086)) +
#   
#   facet_wrap(~disp, labeller = labeller(disp = c("sig5" = "unconstrained", "sig1" = "dispersal constrained"))) +
#   
#  # guides(color = TRUE) +  # Remove color legend
#   
#   labs(x = "", y = "Emerging habitat [\u0394 1980, km <sup>2</sup>]") +
#   
#   theme_light() +
#   
#   theme(strip.text = element_text(size = 15),
#         axis.text = element_text(size = 12),
#         axis.title = element_text(size = 16),
#         axis.title.y = element_markdown(),
#         legend.position = c(0.195, 0.80),
#         legend.key.size = unit(0.6, 'cm'),
#         legend.text = element_text(size=11),
#         legend.title = element_text(size=15),
#         strip.background = element_rect(fill="white", color="black"),
#         strip.text.x = element_text(color="black"),
#         panel.border = element_rect(colour="black", fill=NA))
# 
# print(pl1)

#### slope ####
##### dist #####
sig1 <- NWdist  %>%  filter(disp == "sig1") %>%
  filter(scenario == "1")

cor.test(sig1$year, sig1$dist)
lm(dist ~ year, data = sig1)

sig1 <- NWdist  %>%  filter(disp == "sig1") %>%
  filter(scenario == "2")

cor.test(sig1$year, sig1$dist)
lm(dist ~ year, data = sig1)

sig1 <- NWdist  %>%  filter(disp == "sig1") %>%
  filter(scenario == "3")

cor.test(sig1$year, sig1$dist)
lm(dist ~ year, data = sig1)

#sig5


library(arm)
nsim <- 1000

sig5 <- NWdist  %>%  filter(disp == "sig5") %>%
  filter(scenario == "3")
model <- lm(dist ~ year, data = sig5)
bsim <- sim(model, n.sim = nsim)
apply(coef(bsim), 2, quantile, prob=c(0.025, 0.975))
summary(model)$coefficients[2]

sig1 <- NWdist  %>%  filter(disp == "sig1") %>%
  filter(scenario == "3")
model <- lm(dist ~ year, data = sig1)
bsim <- sim(model, n.sim = nsim)
apply(coef(bsim), 2, quantile, prob=c(0.025, 0.975))
summary(model)$coefficients[2]



##### get data area #####
sig5 <- newSm %>%  filter(disp == "sig5")  %>%
  filter(year == "3")
round(sig5$median)

mean(round(sig5$median))
mean(sig1$median)

#######################

model <- lm(median ~ year, data = sig5)
bsim <- sim(model, n.sim = nsim)
apply(coef(bsim), 2, quantile, prob=c(0.025, 0.975))
summary(model)$coefficients[2]

sig1 <- NWdist  %>%  filter(disp == "sig1") %>%
  filter(scenario == "3")
model <- lm(dist ~ year, data = sig1)
bsim <- sim(model, n.sim = nsim)
apply(coef(bsim), 2, quantile, prob=c(0.025, 0.975))
summary(model)$coefficients[2]

########
sig1 <- NWsum  %>%  filter(disp == "sig1") %>%
  filter(scenario == "1")

cor.test(sig1$year, sig1$dist)
lm(dist ~ year, data = sig1)

sig1 <- NWdist  %>%  filter(disp == "sig1") %>%
  filter(scenario == "2")

cor.test(sig1$year, sig1$dist)
lm(dist ~ year, data = sig1)

sig1 <- NWdist  %>%  filter(disp == "sig1") %>%
  filter(scenario == "3")

cor.test(sig1$year, sig1$dist)
lm(dist ~ year, data = sig1)

#sig5

sig5 <- NWdist  %>%  filter(disp == "sig5") %>%
  filter(scenario == "1")

cor.test(sig5$year, sig5$dist)
lm(dist ~ year, data = sig5)

sig5 <- NWdist  %>%  filter(disp == "sig5") %>%
  filter(scenario == "2")

cor.test(sig5$year, sig5$dist)
lm(dist ~ year, data = sig5)

sig5 <- NWdist  %>%  filter(disp == "sig5") %>%
  filter(scenario == "3")

cor.test(sig5$year, sig5$dist)
lm(dist ~ year, data = sig5)


sig5 <- subset(newSm, newSm$disp == "sig5")

sig1 <- subset(NWsum, NWsum$disp == "sig1")
sig5 <- subset(NWsum, NWsum$disp == "sig5")





#############

sig1 <- subset(newSm, newSm$disp == "sig1")
sig5 <- subset(newSm, newSm$disp == "sig5")

sig1 <- subset(NWsum, NWsum$disp == "sig1")
sig5 <- subset(NWsum, NWsum$disp == "sig5")

cor.test(sig1$year, sig1$median)
cor.test(sig5$year, sig5$median)
cor.test(sig1$scenario, sig1$median)
cor.test(sig5$scenario, sig5$median)

3122/15642  
# Fit linear model for sig1$median vs sig1$year
lm(median ~ year, data = sig1)

# Fit linear model for sig5$median vs sig5$year
lm(median ~ year, data = sig5)

# Fit linear model for sig1$median vs sig1$scenario
lm(median ~ scenario, data = sig1)

# Fit linear model for sig5$median vs sig5$scenario
lm(median ~ scenario, data = sig5)

sig1 <- subset(newSm, newSm$disp == "sig1")
sig5 <- subset(newSm, newSm$disp == "sig5")



sig1$ratio <- sig1$median/sig5$median
mean(sig1$ratio) 
packageVersion("rgbif")
