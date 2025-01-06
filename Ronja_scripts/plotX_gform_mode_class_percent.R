
#packages
library(stars)
library(tidyr)
library(terra)
library(sf)
sf_use_s2(FALSE)
library(dplyr)
library(ggplot2)

#### general ####
data <- "//smb.isipd.dmawi.de/projects/p_ecohealth/projects/Arctic_SDM/data_paper"

traits <- read.csv(glue::glue("{data}/migration_meters_all.csv"))

###########

load(glue::glue("{data}/predArray.rda")) #array
load(glue::glue("{data}/grid_5km.rda")) #grid
load(glue::glue("{data}/spTable.rda")) #spTable

#### res generate h0/result ####
result <- as.data.frame(matrix(nrow=nrow(spTable), ncol=11))
result[,1] <- spTable[,1]
colnames(result) = c("species", "pres", "y2s2", "y2s3", "y2s4", 
                     "y3s1", "y3s2", "y3s3",
                     "y4s1","y4s2","y4s3")


spArray <- predArray[,,,,2]

# Precompute the reference values
reference_values <- spArray[, , 1, 1]

# Calculate the sum for the second column
result[, 2] <- rowSums(reference_values)

# Create a function to calculate percentages with ifelse
calculate_newcells <- function(sp_data, reference, result_column, result_matrix) {
  comparison <- sp_data > reference
  sum_values <- rowSums(reference)
  result_matrix[, result_column] <- rowSums(comparison)
  result_matrix
}

# Calculate percentages for y = 2, 3, 4
for (y in 2:4) {
  start_col <- 3 * y - 4
  result <- calculate_newcells(spArray[, , y, 1], reference_values, start_col + 1, result)
  result <- calculate_newcells(spArray[, , y, 2], reference_values, start_col + 2, result)
  result <- calculate_newcells(spArray[, , y, 3], reference_values, start_col + 3, result)
}

result <- result[-670,]

write.csv2(result, glue::glue("{data}/sp_pres_newcells.csv"), row.names = FALSE)


#### unres generate h0/result ####
result <- as.data.frame(matrix(nrow=nrow(spTable), ncol=11))
result[,1] <- spTable[,1]
colnames(result) = c("species", "pres", "y2s2", "y2s3", "y2s4", 
                     "y3s1", "y3s2", "y3s3",
                     "y4s1","y4s2","y4s3")


spArray <- predArray[,,,,1]

# Precompute the reference values
#reference_values <- spArray[, , 1, 1] take res

# Calculate the sum for the second column
result[, 2] <- rowSums(reference_values)

# Create a function to calculate percentages with ifelse
calculate_newcells <- function(sp_data, reference, result_column, result_matrix) {
  comparison <- sp_data > reference
  sum_values <- rowSums(reference)
  result_matrix[, result_column] <- rowSums(comparison)
  result_matrix
}

# Calculate percentages for y = 2, 3, 4
for (y in 2:4) {
  start_col <- 3 * y - 4
  result <- calculate_newcells(spArray[, , y, 1], reference_values, start_col + 1, result)
  result <- calculate_newcells(spArray[, , y, 2], reference_values, start_col + 2, result)
  result <- calculate_newcells(spArray[, , y, 3], reference_values, start_col + 3, result)
}


write.csv2(result, glue::glue("{data}/sp_pres_newcells_unrespres.csv"), row.names = FALSE)

write.csv2(result2, glue::glue("{data}/sp_pres_newcells_unres.csv"), row.names = FALSE)

#### p1 ####
#h0 <- result
h0 <- read.csv(glue::glue("{data}/sp_pres_newcells.csv"), sep=";")

h0t <- h0 %>%
  left_join(traits, by = "species")

h0t$class[134] <- 6
h0t$mode[134] <- "Dyszoochory"
h0t$gform[134] <- "Shrub"



unres <- read.csv(glue::glue("{data}/sp_pres_newcells_unrespres.csv"), sep=";")
#unres <- read.csv(glue::glue("{data}/sp_pres_newcells_unres.csv"), sep=";")

unres <- unres[-670,]
h0t$unres <- unres$y4s3
h0t <- subset(h0t, h0t$unres > 0)
h0t$percent <- h0t$y4s3/h0t$unres

summary(h0t$pres)
summary(unres$pres)

  
summary(predArray[1,,1,1,1])
summary(predArray[1,,1,1,2])
##### mode ####

result <- h0t 

# result <- as.data.frame(matrix(nrow = 1175, ncol=3))
# result[,1] <- h0t$species
# result[,2] <- h0t$mode
# result[,3] <- h0t$y4s3
# result[,4] <- h0t$percent
# 
# colnames(result) <- c("species", "mode", "new_cells", "percent")

result <- result %>%
  mutate(mode = case_when(
    mode == "Local non-specific dispersal" ~ "local",
    TRUE ~ as.character(mode)  # Keep the original value as a string if none of the above conditions match
  ))

result1 <-  result

result1a <- subset(result1, result1$new_cells > 0)
sum(result1a$mode =="Myrmecochory")
sum(result1$mode =="Myrmecochory")

result1a <- subset(result1, result1$new_cells > 0)
sum(result1a$mode =="local")
sum(result1$mode =="local")

result1 <- subset(result1, result1$mode !="local")

summary(result$percent)

# Create a boxplot
p1 <- ggplot(result1, aes(x = reorder(mode, -percent, FUN = median), y = percent)) +
  geom_boxplot(fill = "#0072B2") +
  labs(#title = expression(paste(Delta, " 2100 constrained")),  # Update labels    
       x = "dispersal mode",
       y = "new cells ratio") +
  theme_minimal()+
  ggtitle("B)")+
  theme(margin(t = 10, r = 10, b = 60, l = 10))+
  theme(   plot.title = element_text(size = 16),
           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

print(p1)

##### class ####
# 
# result <- as.data.frame(matrix(nrow = 1175, ncol=3))
# result[,1] <- h0t$species
# result[,2] <- h0t$class
# result[,3] <- h0t$y4s3
# 
# colnames(result) <- c("species", "class", "new_cells")

result <- result %>%
  mutate(class = case_when(
    class == "Local non-specific dispersal" ~ "local",
    TRUE ~ as.character(class)  # Keep the original value as a string if none of the above conditions match
  ))

result1m <-  result

result1a <- subset(result1m, result1m$new_cells > 0)
sum(result1a$class == "2")
sum(result1m$class == "2")

result1a <- subset(result1m, result1$new_cells > 0)
sum(result1a$mode =="local")
sum(result1$mode =="local")

result1m$classn <- as.numeric(result1m$class)
result1m <-  subset(result1m, result1m$classn > 2)

p1m <- ggplot(result1m, aes(x = reorder(class, -percent, FUN = median), y = percent)) +
  geom_boxplot(fill = "#D95F02") +
  labs(title = "",
       x = "dispersal distance class",
       y = "new cells ratio") +
  theme_minimal()+
  ggtitle("C)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18),
        plot.title = element_text(size = 18))
      
print(p1m)
##### gform ####

# result<- as.data.frame(matrix(nrow = 1175, ncol=3))
# result[,1] <- h0t$species
# result[,2] <- h0t$gform
# result[,3] <- h0t$y4s3
# 
# colnames(result) <- c("species", "gform", "new_cells")
# unique(result$gform)

result <- result %>%
  mutate(gform = case_when(
    gform == "Submerged/floating aquatic" ~ "aquatic",
    TRUE ~ as.character(gform)  # Keep the original value as a string if none of the above conditions match
  ))


result1g <- result 
#summary(as.factor(result1g$gform))

# Load ggplot2 library
library(ggplot2)

# Create a boxplot
p1g <- ggplot(result1g, aes(x = reorder(gform, -percent, FUN = median), y = percent)) +
  geom_boxplot(fill = "#009E73") +
  labs(title = "",
       x = "growth form",
       y = "new cells ratio") +
  theme_minimal()+
  ggtitle("A)")+
  theme(plot.title = element_text(size = 16),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14),
        axis.title.x = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 18))

print(p1g)


library(cowplot)

newcells <- plot_grid(p1g, p1, p1m, nrow=3)
print(newcells)


