
#' This script contains configuration for how the ggplot2 files look like.
#' 
#' References: 
#' - https://rpubs.com/mclaire19/ggplot2-custom-themes


require(ggplot2)
require(ggsci)

# define theme 
theme_ggpr <- function(){ 
  font <- "Arial"  
  # theme_bw(base_size = 14) %+replace%    
  theme_bw(base_size = 14) %+replace%    
    theme(legend.background = element_rect(fill=alpha('white', 0.6), color = NA),
          panel.grid.major = element_line(size = 0.3),  
          panel.grid.minor = element_blank()) 
}
theme_set(theme_ggpr())

# numerical experiments
my_pal <- pal_futurama()(12)
my_pal_numexp <- my_pal[c(4,6)]
legend_font_size_numexp   <- 14
facetrid_font_size_numexp <- 14
title_font_size_numexp    <- 14
gg_base_size_numexp <- 12

message("The file config_figures.R with ggplot2 theme was read.")
