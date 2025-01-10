#### plot 4 combined ####

library(cowplot)

BC <- plot_grid(B, C, ncol = 2, nrow=1, rel_widths = c(2, 1))

DE <- plot_grid(D, E, ncol = 2, nrow= 1, rel_widths = c(2,1))

BCD <- plot_grid(B, C, D, ncol = 3, nrow=1, rel_widths = c(1,2,2))
#ABCDE <- plot_grid(A, BC, DE, ncol= 1, nrow=3)
ABCD <- plot_grid(A, BCD, ncol= 1, nrow=2)

#print(ABCDE)
print(ABCD)
