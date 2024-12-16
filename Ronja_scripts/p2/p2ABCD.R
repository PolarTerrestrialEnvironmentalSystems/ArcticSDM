#### plot 2 ####
#packages
library(tidyverse)
library(sf)
sf_use_s2(FALSE)
library(stars)
library(wesanderson)
library(ggtext)
library(ggpubr)

path <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/data_paper/"

load(glue::glue("{path}/sdm_occurance_metadata_25km_dispersal_99.rda"))

species <- read_csv(glue::glue("{path}/species_sigma1.csv"))


load(glue::glue("{path}/sigma1Array.rda"))
sig1 <- spArray

load(glue::glue("{path}/distconsArray.rda")) #sig5
sig5 <- spArray

### Emerging habitat ###


New <- lapply(1:dim(sig1)[4], function(s) {

  parallel::mclapply(1:dim(sig1)[1], function(i) {

    # cumsum(apply(t(apply(sig1[i,,,s], 1, diff)), 2, function(z) sum(z>0)))

    apply(t(apply(sig1[i,,,s], 1, function(x) x[-1] > x[1])), 2, sum)

  }, mc.cores = 1) %>% Reduce("rbind",.) %>% as_tibble() %>%

    setNames(c("2026", "2056", "2086")) %>%

    pivot_longer(cols = names(.), values_to = 'cells') %>%

    mutate(year = as.numeric(name),

           scenario = s) %>% dplyr::select(scenario, year, cells)

}) %>% Reduce("rbind",.) %>% mutate(disp = "sig1") %>%

  bind_rows(

    lapply(1:dim(sig5)[4], function(s) {

      parallel::mclapply(1:dim(sig5)[1], function(i) {

        apply(t(apply(sig5[i,,,s], 1, function(x) x[-1] > x[1])), 2, sum)

      }, mc.cores = 1) %>% Reduce("rbind",.) %>% as_tibble() %>%

        setNames(c("2026", "2056", "2086")) %>%

        pivot_longer(cols = names(.), values_to = 'cells') %>%

        mutate(year = as.numeric(name),

               scenario = s) %>% dplyr::select(scenario, year, cells)

    }) %>% Reduce("rbind",.) %>% mutate(disp = "sig5")) %>% mutate(area = cells * 25 *25)





newSm <- New %>% filter(area>0) %>% group_by(scenario, year, disp) %>%

  summarise(median = median(area, na.rm = T),

            mean = mean(area, na.rm = T),

            sd = sd(area, na.rm = T),

            lower_50 = quantile(area, probs = 0.25, na.rm = T),

            lower_95 = quantile(area, probs = 0.05, na.rm = T),

            upper_50 = quantile(area, probs = 0.75, na.rm = T),

            upper_95 = quantile(area, probs = 0.975, na.rm = T)

  ) %>% mutate(x_axis = year + (as.numeric(as.character(scenario))-2)*5)



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
  
  scale_fill_manual(values = hcl.colors(6, "heat")[3:1], labels = c("ssp126", "ssp370", "ssp585")) +
 
  # Customizing x-axis breaks and labels
  scale_x_continuous(breaks = c(2026, 2056, 2086), 
                     labels = c("2040", "2070", "2100")) + 

  
  facet_wrap(~disp, labeller = labeller(disp = c("sig5" = "A) unconstrained ","sig1" = "B) constrained"))) +
  
  guides(fill = guide_legend(title = "Scenarios")) +
  
  labs(x = "", y = "Emerging habitat [\u0394 2010, km <sup>2</sup>]") +
  
  theme_light() +
  
  theme(#strip.text = element_text(size = 15),
        
        axis.text = element_text(size = 12),
        
        axis.title = element_text(size = 16),
        
        axis.title.y = element_markdown(),
        
        legend.position = c(0.085, 0.75),#0.095
        
        legend.key.size = unit(0.6, 'cm'),
        
        legend.text = element_text(size=11),
        
        legend.title = element_text(size=15),
        
        strip.background = element_rect(fill = "white", color = "black"),
      #  
       strip.text.x = element_text(color = "black", size=14, hjust=0), # Align text to the left
        
        panel.border = element_rect(colour = "black", fill=NA))


print(pl1)


### Northward shift ###

NW <- lapply(1:dim(sig1)[4], function(s) {
  
  parallel::mclapply(1:dim(sig1)[1], function(i) {
    
    t(apply(sig1[i,,,s], 1, diff)) %>% as_tibble() %>% setNames(c("2026", "2056", "2086")) %>%
      
      mutate(id = 1:nrow(.)) %>% pivot_longer(cols = -id) %>%
      
      filter(value>0)
    
  }, mc.cores = 1) %>% Reduce("rbind", .) %>% mutate(scenario = s)
  
}) %>% Reduce("rbind", .) %>% mutate(disp = "sig1") %>% bind_rows(
  
  lapply(1:dim(sig5)[4], function(s) {
    
    parallel::mclapply(1:dim(sig5)[1], function(i) {
      
      t(apply(sig5[i,,,s], 1, diff)) %>% as_tibble() %>% setNames(c("2026", "2056", "2086")) %>%
        
        mutate(id = 1:nrow(.)) %>% pivot_longer(cols = -id) %>%
        
        filter(value>0)
      
    }, mc.cores = 1) %>% Reduce("rbind", .) %>% mutate(scenario = s)
    
  }) %>% Reduce("rbind", .) %>% mutate(disp = "sig5")
  
)


distGrid <- spArrayMeta$dim2 %>%
  
  mutate(id = 1:nrow(.), dist = as.numeric(st_distance(., st_point(c(0,0)) %>% st_sfc(crs = st_crs(spArrayMeta$dim2))))/1000) %>%
  
  dplyr::select(id, dist) %>% st_drop_geometry()


NWdist <- NW %>% left_join(distGrid, by = "id") %>% mutate(year = as.numeric(name))

##### Ndist ####

NWsum <- NWdist %>% group_by(scenario, year, disp) %>%
  
  summarise(median = median(dist, na.rm = T),
            
            mean = mean(dist, na.rm = T),
            
            sd = sd(dist, na.rm = T),
            
            lower_50 = quantile(dist, probs = 0.25, na.rm = T),
            
            lower_95 = quantile(dist, probs = 0.05, na.rm = T),
            
            upper_50 = quantile(dist, probs = 0.75, na.rm = T),
            
            upper_95 = quantile(dist, probs = 0.975, na.rm = T)
            
  ) %>% mutate(x_axis = year + (as.numeric(as.character(scenario))-2)*5)


NWsum <- NWsum %>%
  mutate(yearnew = case_when(
    year == 2026 ~ 2040,
    year == 2056 ~ 2070,
    year == 2086 ~ 2100
    
  ))


pl2 <- ggplot(NWsum %>% mutate(disp = factor(disp, levels=c("sig5","sig1"))),
              
              aes(x = x_axis, y = median, color = as.factor(scenario), fill = as.factor(scenario))) +
  
  geom_point(size = 5, show.legend = F, alpha = 1) +
  
  geom_segment(mapping = aes(x = x_axis, xend = x_axis,
                             
                             y = lower_50, yend = upper_50), show.legend = F, alpha = 1) +
 
  scale_color_manual(values = hcl.colors(6, "heat")[3:1], labels = c("ssp126", "ssp370", "ssp585")) +
  
  
  # Customizing x-axis breaks and labels
  scale_x_continuous(breaks = c(2026, 2056, 2086), 
                     labels = c("2040", "2070", "2100")) + 
  
  facet_wrap(~disp, labeller = labeller(disp = c("sig5" = "C) unconstrained","sig1" = "D) constrained"))) +
  
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

        panel.border = element_rect(colour = "black", fill=NA),
     strip.background = element_rect(fill = "white", color = "black"),
     #  
     strip.text.x = element_text(color = "black", size=14, hjust=0), # Align text to the left
  )

print(pl2)


pl3 <- pl2 + geom_line(data = NWsum %>% mutate(disp = factor(disp, levels=c("sig5","sig1"))),
          aes(x = x_axis, y = median, color = as.factor(scenario)), 
          size = 0.2, linetype = "dashed", show.legend = FALSE)  # Line color according to scen

print(pl3)
ggarrange(pl1, pl3, nrow = 2)

#### estimations ####

##### share #####
s1 <- subset(newSm, newSm$disp == "sig1")
s5 <- subset(newSm, newSm$disp == "sig5")

s1$ratio <- s1$median/s5$median
mean(s1$ratio)
median(s1$ratio)

##### get max area #####
sig5 <- newSm %>%  filter(disp == "sig5")  %>%
  filter(yearnew == "2100") %>%
  filter(scenario == 3)
sig5$median

##### slope #####

library(arm)
nsim <- 1000

sig <- NWdist  %>%  filter(disp == "sig1") %>% #choose sig1, sig5
  filter(scenario == "1") #choose 1,2,3

model <- lm(dist ~ year, data = sig)
bsim <- sim(model, n.sim = nsim)
apply(coef(bsim), 2, quantile, prob=c(0.025, 0.975))
summary(model)$coefficients[2]

