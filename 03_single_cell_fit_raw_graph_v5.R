################################################################################
# Output of graphs for single cell fit analysis
################################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(tibble)

#Get rid of all currently defined objects
rm(list = ls())

#Set the work directory to where the file is currently located -----------------

workdirectory <- getSrcDirectory(function(x) {x})
setwd(workdirectory)

################################################################################
# Read input data
################################################################################

# Read in combined and fitted single-cell raw data -----------------------------

setwd("../02_01_single_cell_fit_oldschool")

alldata <- as_tibble(read.table(file = "single_cell_fit.csv",
                                sep = ",",
                                header = TRUE,
                                stringsAsFactors = FALSE
                     )
)

# Change data type of wrongly imported variables -------------------------------

alldata <- alldata %>% 
  mutate_at(vars(matches("405|488|561|640")
            ),
            as.double
  )

################################################################################
# Creating fit data for later plotting of single cells and corresponding fit
# with and without sd
################################################################################

# Preparation of unique cells --------------------------------------------------

single.cell.coef <- alldata %>%
  select(cell.unique,
         exp.id,
         exp.analyte.unique,
         analyte,
         analyte.trans,
         matches("^epi.(405|488|561|640).int.nls.lin.(k|d)$"),
         matches("^sd.fit.(405|488|561|640)$"),
         matches("tp.postinj.(405|488|561|640).(start|end)$")
  ) %>%
  droplevels() %>% 
  distinct()

# Each cell now gets a large number of x-axis points ---------------------------
    
single.cell.coef.list <- list()

xaxis <- tibble(xaxis = seq(from = 0, to = 25, length.out = 100))

for (x in c(1:length(single.cell.coef$cell.unique))) {
  single.cell.coef.list[[x]] <- cbind(xaxis,
                                      single.cell.coef[x,]
  )
}

# Bind lists together vertically -----------------------------------------------

single.cell.coef.fit <- as_tibble(do.call(rbind,
                                          single.cell.coef.list
                                  )
)

# Now create fit  --------------------------------------------------------------

single.cell.coef.fit <- single.cell.coef.fit %>% 
  mutate(fit405 = 
         epi.405.int.nls.lin.d * exp(epi.405.int.nls.lin.k * xaxis),
         fit488 = 
         epi.488.int.nls.lin.d * exp(epi.488.int.nls.lin.k * xaxis),
         fit561 = 
         epi.561.int.nls.lin.d * exp(epi.561.int.nls.lin.k * xaxis),
         fit640 = 
         epi.640.int.nls.lin.d * exp(epi.640.int.nls.lin.k * xaxis)
  )

# Create standard deviation graph ----------------------------------------------

single.cell.coef.fit <- single.cell.coef.fit %>% 
  mutate(fit405sdplus  = fit405 + 1.96*sd.fit.405,
         fit405sdminus = fit405 - 1.96*sd.fit.405,
    
         fit488sdplus  = fit488 + 1.96*sd.fit.488,
         fit488sdminus = fit488 - 1.96*sd.fit.488,
    
         fit561sdplus  = fit561 + 1.96*sd.fit.561,
         fit561sdminus = fit561 - 1.96*sd.fit.561,
    
         fit640sdplus  = fit640 + 1.96*sd.fit.640,
         fit640sdminus = fit640 - 1.96*sd.fit.640
  )

################################################################################
# Creating fit data for later plotting of single analytes and corresponding fit
# with and without sd
################################################################################

# Preparation of analytes ------------------------------------------------------

single.analyte.coef <- alldata %>%
  select(exp.analyte.unique,
         exp.id,
         analyte,
         exp.analyte.unique,
         analyte.trans,
         matches("^epi.(405|488|561|640).int.nls.lin.k.out.(sd|mean)$"),
         matches("^tp.postinj.(405|488|561|640).(start|end)$")
  ) %>%
  droplevels() %>%
  distinct()

# Each cell now gets a large number of x-axis points ---------------------------
    
single.analyte.coef.list <- list()

xaxis <- tibble(xaxis = seq(from = 0, to = 25, length.out = 100))

for (x in c(1:length(single.analyte.coef$exp.analyte.unique))) {
  single.analyte.coef.list[[x]] <- cbind(xaxis,
                                         single.analyte.coef[x,]
  )
}

# Bind lists together vertically -----------------------------------------------

single.analyte.coef.fit <- as_tibble(do.call(rbind,
                                          single.analyte.coef.list
                                     )
)

# Now create fit  --------------------------------------------------------------

single.analyte.coef.fit <- single.analyte.coef.fit %>% 
  mutate(meank405 = 
         exp(epi.405.int.nls.lin.k.out.mean * xaxis),
         meank488 = 
         exp(epi.488.int.nls.lin.k.out.mean * xaxis),
         meank561 = 
         exp(epi.561.int.nls.lin.k.out.mean * xaxis),
         meank640 = 
         exp(epi.640.int.nls.lin.k.out.mean * xaxis)
  )

# Create standard deviation graph ----------------------------------------------

single.analyte.coef.fit <- single.analyte.coef.fit %>% 
  mutate(meank405plus = 
           exp((epi.405.int.nls.lin.k.out.mean +
                epi.405.int.nls.lin.k.out.sd * 1.96
           ) * xaxis),
         meank405minus = 
           exp((epi.405.int.nls.lin.k.out.mean -
                epi.405.int.nls.lin.k.out.sd * 1.96
           ) * xaxis),
         meank488plus = 
           exp((epi.488.int.nls.lin.k.out.mean +
                epi.488.int.nls.lin.k.out.sd * 1.96
           ) * xaxis),
         meank488minus = 
           exp((epi.488.int.nls.lin.k.out.mean -
                epi.488.int.nls.lin.k.out.sd * 1.96
           ) * xaxis), 
         meank561plus = 
           exp((epi.561.int.nls.lin.k.out.mean +
                epi.561.int.nls.lin.k.out.sd * 1.96
           ) * xaxis),
         meank561minus = 
           exp((epi.561.int.nls.lin.k.out.mean -
                epi.561.int.nls.lin.k.out.sd * 1.96
           ) * xaxis),
         meank640plus = 
           exp((epi.640.int.nls.lin.k.out.mean +
                epi.640.int.nls.lin.k.out.sd * 1.96
           ) * xaxis),
         meank640minus = 
           exp((epi.640.int.nls.lin.k.out.mean -
                epi.640.int.nls.lin.k.out.sd * 1.96
           ) * xaxis)
  )
    
################################################################################
# Plotting of single-cell graphs
################################################################################         
# 
# imagewidth <- 5
# imageheight <- 2.5
# imageres <- 300
# 
# setwd(workdirectory)
# setwd("../03_single_cell_fit_raw_graph/01_single_cell")
# 
# # For each cell do the following -----------------------------------------------
# 
# for (x in levels(fct_inorder(alldata$cell.unique))){
# 
#   tempsubset1 <- alldata %>%
#     filter(cell.unique == x) %>%
#     droplevels()
# 
#   tempsubset2 <- single.cell.coef.fit %>%
#     filter(cell.unique == x) %>%
#     droplevels()
# 
#   tempsubset2.405 <- tempsubset2 %>%
#   filter(xaxis >= tp.postinj.405.start &
#          xaxis <= tp.postinj.405.end
#   )
# 
#   tempsubset2.488 <- tempsubset2 %>%
#   filter(xaxis >= tp.postinj.488.start &
#          xaxis <= tp.postinj.488.end
#   )
# 
#   tempsubset2.561 <- tempsubset2 %>%
#   filter(xaxis >= tp.postinj.561.start &
#          xaxis <= tp.postinj.561.end
#   )
#   tempsubset2.640 <- tempsubset2 %>%
#   filter(xaxis >= tp.postinj.640.start &
#          xaxis <= tp.postinj.640.end
#   )
# 
#   # Check whether the output file already exists
# 
#   pngname <- paste(tempsubset1$exp.analyte.unique[[1]],
#                    "_",
#                    tempsubset1$pos.id[[1]],
#                    "_",
#                    tempsubset1$cell.id[[1]],
#                    ".png",
#                    sep = ""
#   )
# 
#   if(!file.exists(pngname)){
# 
#     png(filename = pngname,
#         type     = "cairo",
#         units    = "cm",
#         width    = imagewidth,
#         height   = imageheight,
#         res      = imageres
#     )
# 
#     print(
# 
#       ggplot(tempsubset1)
# 
#       # Addition of raw data points --------------------------------------------
# 
#       + geom_point(aes(x      = tp.postinj,
#                        y      = epi.405.int,
#                        colour = "darkblue"
#                    ),
#                    shape = "circle",
#                    size  = 0.176389
#       )
#       + geom_point(aes(x      = tp.postinj,
#                        y      = epi.488.int,
#                        colour = "darkgreen"
#                    ),
#                    shape = "triangle",
#                    size  = 0.176389
#       )
#       + geom_point(aes(x      = tp.postinj,
#                        y      = epi.561.int,
#                        colour = "darkorange"
#                    ),
#                    shape = "square",
#                    size  = 0.176389
#       )
#       + geom_point(aes(x      = tp.postinj,
#                        y      = epi.640.int,
#                        colour = "darkred"
#                    ),
#                    shape = "diamond",
#                    size  = 0.176389
#       )
# 
#       #Add NLS fit 405 ---------------------------------------------------------
# 
#       + {if(!length(tempsubset2.405[[1]]) == 0){
#            geom_line(data = tempsubset2.405,
#                      aes(x = xaxis,
#                          y = fit405,
#                          colour = "darkblue"
#                      ),
#                     size = 0.17
#           )
#          }
#       }
#       + {if(!length(tempsubset2.405[[1]]) == 0){
#         geom_line(data = tempsubset2.405,
#                      aes(x = xaxis,
#                          y = fit405sdplus,
#                          colour = "darkblue"
#                      ),
#                     size = 0.17,
#                     linetype = "dashed",
#                     alpha = 1
#           )
#          }
#       }
# 
#       + {if(!length(tempsubset2.405[[1]]) == 0){
#            geom_line(data = tempsubset2.405,
#                      aes(x = xaxis,
#                          y = fit405sdminus,
#                          colour = "darkblue"
#                      ),
#                     size = 0.17,
#                     linetype = "dashed",
#                     alpha = 1
#           )
#          }
#       }
# 
#       #Add NLS fit 488 ---------------------------------------------------------
# 
#       + {if(!length(tempsubset2.488[[1]]) == 0){
#            geom_line(data = tempsubset2.488,
#                      aes(x = xaxis,
#                          y = fit488,
#                          colour = "darkgreen"
#                      ),
#                     size = 0.17
#           )
#          }
#       }
#       + {if(!length(tempsubset2.488[[1]]) == 0){
#            geom_line(data = tempsubset2.488,
#                      aes(x = xaxis,
#                          y = fit488sdplus,
#                          colour = "darkgreen"
#                      ),
#                     size = 0.17,
#                     linetype = "dashed",
#                     alpha = 1
#           )
#          }
#       }
# 
#       + {if(!length(tempsubset2.488[[1]]) == 0){
#            geom_line(data = tempsubset2.488,
#                      aes(x = xaxis,
#                          y = fit488sdminus,
#                          colour = "darkgreen"
#                      ),
#                     size = 0.17,
#                     linetype = "dashed",
#                     alpha = 1
#           )
#          }
#       }
# 
#       #Add NLS fit 561 ---------------------------------------------------------
# 
#       + {if(!length(tempsubset2.561[[1]]) == 0){
#            geom_line(data = tempsubset2.561,
#                      aes(x = xaxis,
#                          y = fit561,
#                          colour = "darkorange"
#                      ),
#                     size = 0.17
#           )
#          }
#       }
#       + {if(!length(tempsubset2.561[[1]]) == 0){
#            geom_line(data = tempsubset2.561,
#                      aes(x = xaxis,
#                          y = fit561sdplus,
#                          colour = "darkorange"
#                      ),
#                     size = 0.17,
#                     linetype = "dashed",
#                     alpha = 1
#           )
#          }
#       }
# 
#       + {if(!length(tempsubset2.561[[1]]) == 0){
#            geom_line(data = tempsubset2.561,
#                      aes(x = xaxis,
#                          y = fit561sdminus,
#                          colour = "darkorange"
#                      ),
#                     size = 0.17,
#                     linetype = "dashed",
#                     alpha = 1
#           )
#          }
#       }
# 
#       #Add NLS fit 640 ---------------------------------------------------------
# 
#       + {if(!length(tempsubset2.640[[1]]) == 0){
#            geom_line(data = tempsubset2.640,
#                      aes(x = xaxis,
#                          y = fit640,
#                          colour = "darkred"
#                      ),
#                     size = 0.17
#           )
#          }
#       }
#       + {if(!length(tempsubset2.640[[1]]) == 0){
#            geom_line(data = tempsubset2.640,
#                      aes(x = xaxis,
#                          y = fit640sdplus,
#                          colour = "darkred"
#                      ),
#                     size = 0.17,
#                     linetype = "dashed",
#                     alpha = 1
#           )
#          }
#       }
# 
#       + {if(!length(tempsubset2.640[[1]]) == 0){
#            geom_line(data = tempsubset2.640,
#                      aes(x = xaxis,
#                          y = fit640sdminus,
#                          colour = "darkred"
#                      ),
#                     size = 0.17,
#                     linetype = "dashed",
#                     alpha = 1
#           )
#          }
#       }
# 
#       # Adding fit coefficients ------------------------------------------------
# 
#       # + geom_text(aes(x     = 0,
#       #                 y     = 0,
#       #                 label = s488
#       #             ),
#       #             size  =  15 / .pt * (5/14),
#       #             vjust = "bottom",
#       #             hjust = "left"
#       # )
# 
#       + theme_bw()#base_size = 3.88056)
#       + labs(#title    = paste(tempsubset1$exp.analyte.unique[[1]],
#              #                 "_",
#              #                 tempsubset1$pos.id[[1]],
#              #                 "_",
#              #                 tempsubset1$cell.id[[1]],
#              #                 ".png",
#              #                  sep = ""
#              #),
#              x        = "Time [h]",
#              y        = "Fluor. [AU]"
#       )
# 
#       + coord_cartesian(xlim = c(0, 15)
#       )
#       #+ scale_y_log10()
#       + scale_colour_identity("Channel",
#                               labels = c("darkblue" = "405 nm",
#                                          "darkgreen" = "488 nm",
#                                          "darkorange" = "561 nm",
#                                          "darkred" = "640 nm"
#                               ),
#                               guide = "legend"
#       )
#       + theme(text = element_text(size = 15 / .pt),
#               axis.line = element_line(colour = "black"),
#               axis.text = element_text(colour = "black"),
#               panel.border = element_rect(colour = "black", size = 0.75),
#               axis.ticks = element_line(colour = "black"),
#               legend.position = "none"
#           )
#     )
# 
#     dev.off()
#   }
# }

################################################################################
# Plotting of RAW data per experiment and analyte
# 
# This section will create one compilation of all single cell graphs for each
# analyte and experiment
################################################################################         

# original
imagewidth <- 26
imageheight <- 11

# alternative for large datasets
imagewidth <- 30
imageheight <- 16

imageres <- 300

setwd(workdirectory)
setwd("../03_single_cell_fit_raw_graph/02_single_cell_per_exp_per_analyte_grouped")

for (x in levels(fct_inorder(alldata$exp.id))){

  tempsubset1 <- alldata %>%
    filter(exp.id == x) %>%
    droplevels()

  tempsubset2 <- single.cell.coef.fit %>%
    filter(exp.id == x) %>%
    droplevels()

  for (y in levels(as.factor(tempsubset1$analyte))){
    
    # Check if image already exists
    pngname <- paste(x, "_", y, ".png", sep = "")
    
    if(!file.exists(pngname)){
      
      tempsubset11 <- tempsubset1 %>%
        filter(analyte == y)
  
      tempsubset22 <- tempsubset2 %>%
        filter(analyte == y)
  
      tempsubset222.405 <- tempsubset22 %>%
        filter(xaxis >= tp.postinj.405.start &
               xaxis <= tp.postinj.405.end
        )
  
      tempsubset222.488 <- tempsubset22 %>%
        filter(xaxis >= tp.postinj.488.start &
               xaxis <= tp.postinj.488.end
        )
  
      tempsubset222.561 <- tempsubset22 %>%
        filter(xaxis >= tp.postinj.561.start &
               xaxis <= tp.postinj.561.end
        )
  
      tempsubset222.640 <- tempsubset22 %>%
        filter(xaxis >= tp.postinj.640.start &
               xaxis <= tp.postinj.640.end
        )
  
      png(filename = pngname,
          type     = "cairo",
          units    = "cm",
          width    = imagewidth,
          height   = imageheight,
          res      = imageres
      )
  
      print(
  
        ggplot(tempsubset11)
  
        # Addition of raw data points --------------------------------------------
  
        + geom_point(aes(x      = tp.postinj,
                         y      = epi.405.int,
                         colour = "darkblue"
                     ),
                     shape = "circle",
                     size  = 0.176389
        )
  
        + geom_point(aes(x      = tp.postinj,
                         y      = epi.488.int,
                         colour = "darkgreen"
                     ),
                     shape = "triangle",
                     size  = 0.176389
        )
  
        + geom_point(aes(x      = tp.postinj,
                         y      = epi.561.int,
                         colour = "darkorange"
                     ),
                     shape = "square",
                     size  = 0.176389
        )
  
        + geom_point(aes(x      = tp.postinj,
                         y      = epi.640.int,
                         colour = "darkred"
                     ),
                     shape = "diamond",
                     size  = 0.176389
        )
  
        #Add NLS fit 405 ---------------------------------------------------------
  
        + {if(!length(tempsubset222.405[[1]]) == 0){
             geom_line(data = tempsubset222.405,
                       aes(x = xaxis,
                           y = fit405,
                           colour = "darkblue"
                       ),
                      size = 0.17
            )
           }
        }
  
        + {if(!length(tempsubset222.405[[1]]) == 0){
          geom_line(data = tempsubset222.405,
                       aes(x = xaxis,
                           y = fit405sdplus,
                           colour = "darkblue"
                       ),
                      size = 0.17,
                      linetype = "dashed",
                      alpha = 1
            )
           }
        }
  
        + {if(!length(tempsubset222.405[[1]]) == 0){
             geom_line(data = tempsubset222.405,
                       aes(x = xaxis,
                           y = fit405sdminus,
                           colour = "darkblue"
                       ),
                      size = 0.17,
                      linetype = "dashed",
                      alpha = 1
            )
           }
        }
  
        #Add NLS fit 488 ---------------------------------------------------------
  
  
        + {if(!length(tempsubset222.488[[1]]) == 0){
             geom_line(data = tempsubset222.488,
                       aes(x = xaxis,
                           y = fit488,
                           colour = "darkgreen"
                       ),
                      size = 0.17
            )
           }
        }
  
        + {if(!length(tempsubset222.488[[1]]) == 0){
             geom_line(data = tempsubset222.488,
                       aes(x = xaxis,
                           y = fit488sdplus,
                           colour = "darkgreen"
                       ),
                      size = 0.17,
                      linetype = "dashed",
                      alpha = 1
            )
           }
        }
  
        + {if(!length(tempsubset222.488[[1]]) == 0){
             geom_line(data = tempsubset222.488,
                       aes(x = xaxis,
                           y = fit488sdminus,
                           colour = "darkgreen"
                       ),
                      size = 0.17,
                      linetype = "dashed",
                      alpha = 1
            )
           }
        }
  
        #Add NLS fit 561 ---------------------------------------------------------
  
  
        + {if(!length(tempsubset222.561[[1]]) == 0){
             geom_line(data = tempsubset222.561,
                       aes(x = xaxis,
                           y = fit561,
                           colour = "darkorange"
                       ),
                      size = 0.17
            )
           }
        }
  
        + {if(!length(tempsubset222.561[[1]]) == 0){
             geom_line(data = tempsubset222.561,
                       aes(x = xaxis,
                           y = fit561sdplus,
                           colour = "darkorange"
                       ),
                      size = 0.17,
                      linetype = "dashed",
                      alpha = 1
            )
           }
        }
  
        + {if(!length(tempsubset222.561[[1]]) == 0){
             geom_line(data = tempsubset222.561,
                       aes(x = xaxis,
                           y = fit561sdminus,
                           colour = "darkorange"
                       ),
                      size = 0.17,
                      linetype = "dashed",
                      alpha = 1
            )
           }
        }
  
        #Add NLS fit 640 ---------------------------------------------------------
  
  
        + {if(!length(tempsubset222.640[[1]]) == 0){
             geom_line(data = tempsubset222.640,
                       aes(x = xaxis,
                           y = fit640,
                           colour = "darkred"
                       ),
                      size = 0.17
            )
           }
        }
  
        + {if(!length(tempsubset222.640[[1]]) == 0){
             geom_line(data = tempsubset222.640,
                       aes(x = xaxis,
                           y = fit640sdplus,
                           colour = "darkred"
                       ),
                      size = 0.17,
                      linetype = "dashed",
                      alpha = 1
            )
           }
        }
  
        + {if(!length(tempsubset222.640[[1]]) == 0){
             geom_line(data = tempsubset222.640,
                       aes(x = xaxis,
                           y = fit640sdminus,
                           colour = "darkred"
                       ),
                      size = 0.17,
                      linetype = "dashed",
                      alpha = 1
            )
           }
        }
        # Pretty settings --------------------------------------------------------
  
        + facet_wrap(~fct_inorder(as.factor(cell.unique)), scale = "free_y")
        + theme_bw(base_size = 3.88056)
        + labs(title    = paste(x, "_", y, sep = ""),
               x        = "Time [h]",
               y        = "Fluor. [AU]"
        )
  
  
        #+ scale_y_continuous(trans = "log")
        + scale_colour_identity("Channel",
                                labels = c("darkblue" = "405 nm",
                                           "darkgreen" = "488 nm",
                                           "darkorange" = "561 nm",
                                           "darkred" = "640 nm"
                                ),
                                guide = "legend"
        )
        + coord_cartesian(xlim = c(0,12),
                          ylim = c(0,NA)
        )
        + theme(strip.background = element_rect(color = "white",
                                                fill="white"
                ),
                text = element_text(size = 15 / .pt),
                axis.line = element_line(colour = "black"),
                axis.text = element_text(colour = "black"),
                panel.border = element_rect(colour = "black", size = 1),
                axis.ticks = element_line(colour = "black"),
                legend.position = "bottom"
        )
      )
  
      dev.off()
    }  
  }
}

################################################################################
# Plotting single analytes including single cell rawdata, mean and 
################################################################################         

imagewidth <- 3
imageheight <- 2.5
imageres <- 300

setwd(workdirectory)
setwd("../03_single_cell_fit_raw_graph/03_single_cell_per_exp_analyte_fit_raw")

# For each cell do the following -----------------------------------------------

for (x in levels(fct_inorder(alldata$exp.analyte.unique))){
  
  tempsubset1 <- alldata %>%
    filter(exp.analyte.unique == x) %>%
    droplevels()
  
  tempsubset2 <- single.analyte.coef.fit %>%
    filter(exp.analyte.unique == x) %>%
    droplevels()
    
  pngname <- paste(tempsubset1$exp.analyte.unique[[1]],
                              ".png",
                              sep = ""
  )
  
  if(!file.exists(pngname)){
  
    png(filename = paste(tempsubset1$exp.analyte.unique[[1]],
                              ".png",
                              sep = ""
        ),
        type     = "cairo",
        units    = "cm",
        width    = imagewidth,
        height   = imageheight,
        res      = imageres
    )
  
      print(
        
        ggplot(tempsubset1)
  
        + geom_line(aes(x = tp.postinj,
                         y = epi.405.int / epi.405.int.nls.lin.d, 
                         colour = "darkblue",
                         group = cell.unique
                     ),
                     size = 0.14,
                     alpha = 0.2
        )
        
        + geom_line(aes(x = tp.postinj,
                         y = epi.488.int / epi.488.int.nls.lin.d,
                         colour = "darkgreen",
                         group = cell.unique
                     ),
                     size = 0.14,
                     alpha = 0.2
        )
        
        + geom_line(aes(x = tp.postinj,
                         y = epi.561.int / epi.561.int.nls.lin.d,
                         colour = "darkorange",
                         group = cell.unique
                     ),
                     size = 0.14,
                     alpha = 0.2
        )
        
        + geom_line(aes(x = tp.postinj,
                         y = epi.640.int / epi.640.int.nls.lin.d,
                         colour = "darkred",
                         group = cell.unique
                     ),
                     size = 0.14,
                     alpha = 0.2
        )
  
        #Add NLS fit 405 ---------------------------------------------------------
        
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank405,
              colour = "black"
          ),
          size = 0.9
        )
        
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank405,
              colour = "blue"
          ),
          size = 0.3
        )
        
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank405minus,
        #       colour = "darkblue"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        # 
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank405plus,
        #       colour = "darkblue"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
  
        #Add NLS fit 488 ---------------------------------------------------------
        
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank488,
              colour = "black"
          ),
          size = 0.9
        )
      
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank488,
              colour = "green"
          ),
          size = 0.3
        )
          
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank488minus,
        #       colour = "darkgreen"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        # 
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank488plus,
        #       colour = "darkgreen"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        # 
        #Add NLS fit 561 ---------------------------------------------------------
        
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank561,
              colour = "black"
          ),
          size = 0.9
        )
     
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank561,
              colour = "darkorange"
          ),
          size = 0.3
        )
       
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank561minus,
        #       colour = "darkorange"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        # 
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank561plus,
        #       colour = "darkorange"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        
        #Add NLS fit 640 ---------------------------------------------------------
        
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank640,
              colour = "black"
          ),
          size = 0.9
        )
        
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank640,
              colour = "red"
          ),
          size = 0.3
        )
          
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank640minus,
        #       colour = "darkred"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        # 
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank640plus,
        #       colour = "darkred"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        
        + theme_bw(base_size = 3.88056)
        + labs(#title    = paste(tempsubset1$exp.analyte.unique[[1]]),
             x        = "Time [h]",
             y        = "Fluor. [AU]"
          )
        + scale_colour_identity("Channel",
                                labels = c("darkblue" = "405 nm",
                                           "darkgreen" = "488 nm",
                                           "darkorange" = "561 nm",
                                           "darkred" = "640 nm"
                                ),
                                guide = "legend"
        )
  
      + coord_cartesian(xlim = c(0,12),
                        ylim = c(0,2),
                        expand = T
      )
      + scale_x_continuous(breaks=c(seq(0,12,2)))
      + theme(strip.background = element_rect(color = "white",
                                              fill="white"
              ),
              text = element_text(size = 15 / .pt),
              axis.line = element_line(colour = "black"),
              axis.text = element_text(colour = "black"),
              panel.border = element_rect(colour = "black", size = 1),
              axis.ticks = element_line(colour = "black"),
              legend.position = "bottom"
      )
      
      + theme(legend.position = "none")
        #+ scale_y_continuous(trans = "log")
      #+ scale_y_log10()
    )
  
    dev.off()
  }
}

################################################################################
# Plotting scaled raw data and statistical cell, now only with different file
# names, in order to multiple injections of the same analyte
################################################################################         

imagewidth <- 3
imageheight <- 2.5
imageres <- 300

setwd(workdirectory)
setwd(paste("../",
            "03_single_cell_fit_raw_graph/",
            "04_single_cell_per_analyte_incl_all_exp_fit_raw",
            sep = ""
      )
)

# For each cell do the following -----------------------------------------------

for (x in levels(fct_inorder(alldata$exp.analyte.unique))){
  
  tempsubset1 <- alldata %>%
    filter(exp.analyte.unique == x) %>%
    droplevels()
  
  tempsubset2 <- single.analyte.coef.fit %>%
    filter(exp.analyte.unique == x) %>%
    droplevels()
    
  pngname <- paste(tempsubset1$analyte.trans[[1]],
                   "_",
                   tempsubset1$exp.id[[1]],
                   ".png",
                   sep = ""
  )
  
  if(!file.exists(pngname)){
  
    png(filename = pngname,
        type     = "cairo",
        units    = "cm",
        width    = imagewidth,
        height   = imageheight,
        res      = imageres
    )
  
      print(
        
        ggplot(tempsubset1)
  
        + geom_line(aes(x = tp.postinj,
                         y = epi.405.int / epi.405.int.nls.lin.d, 
                         colour = "darkblue",
                         group = cell.unique
                     ),
                     size = 0.14,
                     alpha = 0.2
        )
        
        + geom_line(aes(x = tp.postinj,
                         y = epi.488.int / epi.488.int.nls.lin.d,
                         colour = "darkgreen",
                         group = cell.unique
                     ),
                     size = 0.14,
                     alpha = 0.2
        )
        
        + geom_line(aes(x = tp.postinj,
                         y = epi.561.int / epi.561.int.nls.lin.d,
                         colour = "darkorange",
                         group = cell.unique
                     ),
                     size = 0.14,
                     alpha = 0.2
        )
        
        + geom_line(aes(x = tp.postinj,
                         y = epi.640.int / epi.640.int.nls.lin.d,
                         colour = "darkred",
                         group = cell.unique
                     ),
                     size = 0.14,
                     alpha = 0.2
        )
  
        #Add NLS fit 405 ---------------------------------------------------------
        
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank405,
              colour = "black"
          ),
          size = 0.9
        )
        
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank405,
              colour = "blue"
          ),
          size = 0.3
        )
        
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank405minus,
        #       colour = "darkblue"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        # 
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank405plus,
        #       colour = "darkblue"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
  
        #Add NLS fit 488 ---------------------------------------------------------
        
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank488,
              colour = "black"
          ),
          size = 0.9
        )
      
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank488,
              colour = "green"
          ),
          size = 0.3
        )
          
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank488minus,
        #       colour = "darkgreen"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        # 
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank488plus,
        #       colour = "darkgreen"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        # 
        #Add NLS fit 561 ---------------------------------------------------------
        
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank561,
              colour = "black"
          ),
          size = 0.9
        )
     
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank561,
              colour = "darkorange"
          ),
          size = 0.3
        )
       
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank561minus,
        #       colour = "darkorange"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        # 
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank561plus,
        #       colour = "darkorange"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        
        #Add NLS fit 640 ---------------------------------------------------------
        
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank640,
              colour = "black"
          ),
          size = 0.9
        )
        
        + geom_line(data = tempsubset2,
          aes(x = xaxis,
              y = meank640,
              colour = "red"
          ),
          size = 0.3
        )
          
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank640minus,
        #       colour = "darkred"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        # 
        # + geom_line(data = tempsubset2,
        #   aes(x = xaxis,
        #       y = meank640plus,
        #       colour = "darkred"
        #   ),
        #   linetype = "dashed",
        #   alpha = 0.5,
        #   size = 0.17
        # )
        
        + theme_bw(base_size = 3.88056)
        + labs(#title    = paste(tempsubset1$exp.analyte.unique[[1]]),
             x        = "Time [h]",
             y        = "Fluor. [AU]"
          )
        + scale_colour_identity("Channel",
                                labels = c("darkblue" = "405 nm",
                                           "darkgreen" = "488 nm",
                                           "darkorange" = "561 nm",
                                           "darkred" = "640 nm"
                                ),
                                guide = "legend"
        )
  
      + coord_cartesian(xlim = c(0,12),
                        ylim = c(0,2),
                        expand = T
      )
      + scale_x_continuous(breaks=c(seq(0,12,2)))
      + theme(strip.background = element_rect(color = "white",
                                              fill="white"
              ),
              text = element_text(size = 15 / .pt),
              axis.line = element_line(colour = "black"),
              axis.text = element_text(colour = "black"),
              panel.border = element_rect(colour = "black", size = 1),
              axis.ticks = element_line(colour = "black"),
              legend.position = "bottom"
      )
      
      + theme(legend.position = "none")
        #+ scale_y_continuous(trans = "log")
      #+ scale_y_log10()
    )
  
    dev.off()
  }
}

