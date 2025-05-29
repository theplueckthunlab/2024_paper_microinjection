################################################################################
# SINGLE CELL FIT V12
#
# V12: Removed removal of channels based on analyte translation document.
################################################################################

library(tidyverse)
library(readxl)

# Get rid of all currently defined objects -------------------------------------

rm(list = ls())

#Set the work directory to where the file is currently located -----------------

workdirectory <- getSrcDirectory(function(x) {x})
setwd(workdirectory)

max = function(x, na.rm = TRUE) {
  if(all(is.na(x))) return(x[1])
  return(base::max(x, na.rm = na.rm))
}

min = function(x, na.rm = TRUE) {
  if(all(is.na(x))) return(x[1])
  return(base::min(x, na.rm = na.rm))
}

################################################################################
# Read single cell input data
################################################################################

# Read in all single-cell raw data CSV files -----------------------------------

# Find all filenames ending with .csv
allfilenames <- list.files(path = "./../01_raw_data_input",
                           pattern = "*.csv",
                           full.names = T
)

# Initiate list that will then be filled in a for loop
allfiles <- list()

for (x in c(1:length(allfilenames))) {
  allfiles[[x]] <- read_csv(file = allfilenames[[x]],
                            col_types = cols(exp.id      = col_character(),
   	 																			   pos.id      = col_double(),
   	 																				 cell.id     = col_double(),
   	 																			   tp.id       = col_double(),
   	 																			   tp.posix    = col_double(),
   	 																			   tp.postinj  = col_double(),
   	 																			   epi.405.int = col_double(),
   	 																			   epi.488.int = col_double(),
   	 																			   epi.561.int = col_double(),
   	 																			   epi.640.int = col_double(),
   	 																			   analyte     = col_character()
                            )
  )
}

alldata <- bind_rows(allfiles)

# Create identifiers for split-apply-combine operations ------------------------

alldata <- alldata %>%
  
  # For single cells generally (for fitting)
  unite(cell.unique,
        exp.id,
        pos.id,
        cell.id,
        sep = "_",
        remove = F
  ) %>%
  
  # For splitting single experiments and single analytes
  unite(exp.analyte.unique,
        exp.id,
        analyte,
        sep = "_",
        remove = F
  )

################################################################################
# IMPORT OF ANALYTE METADATA
#
# Add column containing global analyte names E.g. 0002_GS-eGFP and 00029_GS-eGFP
# are translated to GS-eGFP, thereby dropping the protein ID.
#
# During some experiments, channels that physically cannot record any data, e.g.
# DAPI when eGFP without injection marker is recorded.
################################################################################

translation.read <- read_excel("conversion_analytes_v17.xlsx")

# Do all analytes have conversion data? ----------------------------------------

allanalytes     <- sort(levels(as.factor(alldata$analyte)))
allanalytetrans <- sort(levels(as.factor(translation.read$analytefactor)))

# This is super important > do all analytes have a common name after
# translation?

nottranslated <- setdiff(allanalytes,
												 allanalytetrans
)

if(length(nottranslated) > 0) {
  stop("The following analytes were not translated. Therefore, no channel ",
       "information available.\n",
       paste(nottranslated,
             "\n",
             sep = ""
       )
  )
}

# Do all conversion names exist as analytes? -----------------------------------

# Why do I have data names in the translation sheet that are not found in the
# list of analytes? One reason might be that I changed the names at some point
# in the data extraction script found in every single-experiment folder.
# > investigate further

notfoundinalldata <- setdiff(allanalytetrans,
														 allanalytes
)

if(length(notfoundinalldata) > 0) {
  stop("The following analytes found in translation sheet to not exist in the ",
       "real data.\n",
       paste(notfoundinalldata,
             "\n",
             sep = ""
       )
  )
}

# Add translations to alldata file ---------------------------------------------

# If the analyte is not translated, use sample name
alldata <- alldata %>%
  mutate(analyte.trans = analyte)

for (x in seq_along(translation.read[[1]])){
	
  # Which rows contain the translation in question?
  matches <- which(alldata$analyte == as.character(translation.read[x,1]))
  alldata[matches, "analyte.trans"] <- translation.read[x,2]
}

# Adding information about analyte fluorescent channels ------------------------     # EDIT
# This section is actually useful, work around experiments with same analyte but 
# different acquisition contexts > e.g. dextran in different channels, or no
# dextran for example

# adding.analyte.fluor <- translation.read %>% 
#   ungroup() %>% 
#   select(c(2:6)) %>% 
#   rename(analyte.trans = analytetranslation) %>% 
#   distinct()
# 
# alldata <- alldata %>%
#   full_join(adding.analyte.fluor, by = "analyte.trans")

# Set channels that are not coming from analyte to NA --------------------------

# alldata <- alldata %>%
#   mutate(epi.405.int = ifelse(y405 == 1,
#                               epi.405.int,
#                               NA
#          ),
#          epi.488.int = ifelse(y488 == 1,
#                               epi.488.int,
#                               NA
#          ),
#          epi.561.int = ifelse(y561 == 1,
#                               epi.561.int,
#                               NA
#          ),
#          epi.640.int = ifelse(y640 == 1,
#                               epi.640.int,
#                               NA
#          )
# )

################################################################################
# Remove datapoints too low from fitting procedure
################################################################################

################################################################################
# Remove certain analytes
################################################################################

# alldata <- alldata %>%
#   filter(!exp.analyte.unique == "2017-12-04_0065_GS-mNG") %>%
    # different settings than later on in the experiments, confocal data at end of experiment available


################################################################################
# When was the first and last timepoint acquired after microinjection
# This is performed for each cell and channel and added to alldata
################################################################################

alldata <- alldata %>% 
  group_by(cell.unique) %>%
  
  # No signal in channel > put NA into channel time column ---------------------
  
  mutate(tp.postinj.405 = ifelse(is.na(epi.405.int),
                                 epi.405.int,
                                 tp.postinj
         ),
         tp.postinj.488 = ifelse(is.na(epi.488.int),
                                 epi.488.int,
                                 tp.postinj
         ),
         tp.postinj.561 = ifelse(is.na(epi.561.int),
                                 epi.561.int,
                                 tp.postinj
         ),
         tp.postinj.640 = ifelse(is.na(epi.640.int),
                                 epi.640.int,
                                 tp.postinj
         )
  ) %>%
  
  # Extract start end from each channel ----------------------------------------

  mutate(tp.postinj.405.start = min(tp.postinj.405, na.rm = T),
         tp.postinj.405.end   = max(tp.postinj.405, na.rm = T),
         tp.postinj.488.start = min(tp.postinj.488, na.rm = T),
         tp.postinj.488.end   = max(tp.postinj.488, na.rm = T),
         tp.postinj.561.start = min(tp.postinj.561, na.rm = T),
         tp.postinj.561.end   = max(tp.postinj.561, na.rm = T),
         tp.postinj.640.start = min(tp.postinj.640, na.rm = T),
         tp.postinj.640.end   = max(tp.postinj.640, na.rm = T)
  ) %>%
  
  # Set to NA if min/max causes infinite values --------------------------------

  mutate(tp.postinj.405.start = ifelse(is.infinite(tp.postinj.405.start),
  	                                     NA,
  	                                     tp.postinj.405.start
         ),
				 tp.postinj.405.end   = ifelse(is.infinite(tp.postinj.405.end),
 	                                     NA,
 	                                     tp.postinj.405.end
         ),
				 tp.postinj.488.start = ifelse(is.infinite(tp.postinj.488.start),
 	                                     NA,
 	                                     tp.postinj.488.start
         ),
				 tp.postinj.488.end   = ifelse(is.infinite(tp.postinj.488.end),
 	                                     NA,
 	                                     tp.postinj.488.end
         ),
				 tp.postinj.561.start = ifelse(is.infinite(tp.postinj.561.start),
 	                                     NA,
 	                                     tp.postinj.561.start
         ), 
				 tp.postinj.561.end   = ifelse(is.infinite(tp.postinj.561.end),
 	                                     NA,
 	                                     tp.postinj.561.end
         ), 
				 tp.postinj.640.start = ifelse(is.infinite(tp.postinj.640.start),
 	                                     NA,
 	                                     tp.postinj.640.start
         ), 
				 tp.postinj.640.end   = ifelse(is.infinite(tp.postinj.640.end),
 	                                     NA,
 	                                     tp.postinj.640.end
         )
  )

################################################################################
# 405 Fitting procedure
################################################################################

# Logarithmic linearization and linear regression ------------------------------

exp.fit.405.lin <- alldata %>%
  filter(is.na(epi.405.int) == FALSE) %>%
  group_by(cell.unique) %>%
  do(epi.405.int.fit.lin.d = 
     exp(as.numeric(lm(log(epi.405.int) ~ tp.postinj, data = .)[[1]][[1]])),
     epi.405.int.fit.lin.k =
     as.numeric(lm(log(epi.405.int) ~ tp.postinj, data = .)[[1]][[2]])
  )

exp.fit.405.lin[[2]] <- as.numeric(exp.fit.405.lin[[2]])
exp.fit.405.lin[[3]] <- as.numeric(exp.fit.405.lin[[3]])

# Now we join these newly found parameters into the alldata table --------------

alldata <- alldata %>%
  full_join(exp.fit.405.lin, by = "cell.unique")

# They serves as starting points for NLS procedure -----------------------------

exp.fit.405.nls <- alldata %>%
  filter(is.na(epi.405.int) == F) %>%
  filter(is.na(epi.405.int.fit.lin.d) == F) %>%
  filter(is.na(epi.405.int.fit.lin.k) == F) %>%
  group_by(cell.unique) %>%
  do(epi.405.int.fit.nls.k = 
     as.numeric(try(summary(nls(formula = epi.405.int ~ ((exp(k*tp.postinj)*C)),
                                data = .,
                                start = list(k = unique(.$epi.405.int.fit.lin.k),
                                             C = unique(.$epi.405.int.fit.lin.d)
                                )
                            )
                    )$parameters[[1]]
                )
     ),
     epi.405.int.fit.nls.d =
     as.numeric(try(summary(nls(formula = epi.405.int ~ (exp(k*tp.postinj)*C),
                                data = .,
                                start = list(k = unique(.$epi.405.int.fit.lin.k),
                                             C = unique(.$epi.405.int.fit.lin.d)
                                )
                            )
                    )$parameters[[2]]
                )
     )
  )

exp.fit.405.nls[[2]] <- as.numeric(exp.fit.405.nls[[2]])
exp.fit.405.nls[[3]] <- as.numeric(exp.fit.405.nls[[3]])

# Now we join these newly found parameters into the alldata table --------------

alldata <- alldata %>%
  full_join(exp.fit.405.nls, by = "cell.unique")

################################################################################
# 488 Fitting procedure
################################################################################

# Logarithmic linearization and linear regression ------------------------------

exp.fit.488.lin <- alldata %>%
  filter(is.na(epi.488.int) == FALSE) %>%
  group_by(cell.unique) %>%
  do(epi.488.int.fit.lin.d = 
     exp(as.numeric(lm(log(epi.488.int) ~ tp.postinj, data = .)[[1]][[1]])),
     epi.488.int.fit.lin.k =
     as.numeric(lm(log(epi.488.int) ~ tp.postinj, data = .)[[1]][[2]])
  )

exp.fit.488.lin[[2]] <- as.numeric(exp.fit.488.lin[[2]])
exp.fit.488.lin[[3]] <- as.numeric(exp.fit.488.lin[[3]])

# Now we join these newly found parameters into the alldata table --------------

alldata <- alldata %>%
  full_join(exp.fit.488.lin, by = "cell.unique")

# They serves as starting points for NLS procedure -----------------------------

exp.fit.488.nls <- alldata %>%
  filter(is.na(epi.488.int) == F) %>%
  filter(is.na(epi.488.int.fit.lin.d) == F) %>%
  filter(is.na(epi.488.int.fit.lin.k) == F) %>%
  group_by(cell.unique) %>%
  do(epi.488.int.fit.nls.k = 
     as.numeric(try(summary(nls(formula = epi.488.int ~ (exp(k*tp.postinj)*C),
                                data = .,
                                start = list(k = unique(.$epi.488.int.fit.lin.k),
                                             C = unique(.$epi.488.int.fit.lin.d)
                                )
                            )
                    )$parameters[[1]]
                )
     ),
     epi.488.int.fit.nls.d =
     as.numeric(try(summary(nls(formula = epi.488.int ~ (exp(k*tp.postinj)*C),
                                data = .,
                                start = list(k = unique(.$epi.488.int.fit.lin.k),
                                             C = unique(.$epi.488.int.fit.lin.d)
                                )
                            )
                    )$parameters[[2]]
                )
     )
  )

exp.fit.488.nls[[2]] <- as.numeric(exp.fit.488.nls[[2]])
exp.fit.488.nls[[3]] <- as.numeric(exp.fit.488.nls[[3]])

# Now we join these newly found parameters into the alldata table --------------

alldata <- alldata %>%
  full_join(exp.fit.488.nls, by = "cell.unique")

################################################################################
# 561 Fitting procedure
################################################################################

# Logarithmic linearization and linear regression ------------------------------

exp.fit.561.lin <- alldata %>%
  filter(is.na(epi.561.int) == FALSE) %>%
  group_by(cell.unique) %>%
  do(epi.561.int.fit.lin.d = 
     exp(as.numeric(lm(log(epi.561.int) ~ tp.postinj, data = .)[[1]][[1]])),
     epi.561.int.fit.lin.k =
     as.numeric(lm(log(epi.561.int) ~ tp.postinj, data = .)[[1]][[2]])
  )

exp.fit.561.lin[[2]] <- as.numeric(exp.fit.561.lin[[2]])
exp.fit.561.lin[[3]] <- as.numeric(exp.fit.561.lin[[3]])

# Now we join these newly found parameters into the alldata table --------------

alldata <- alldata %>%
  full_join(exp.fit.561.lin, by = "cell.unique")

# They serves as starting points for NLS procedure -----------------------------

exp.fit.561.nls <- alldata %>%
  filter(is.na(epi.561.int) == F) %>%
  filter(is.na(epi.561.int.fit.lin.d) == F) %>%
  filter(is.na(epi.561.int.fit.lin.k) == F) %>%
  group_by(cell.unique) %>%
  do(epi.561.int.fit.nls.k = 
     as.numeric(try(summary(nls(formula = epi.561.int ~ (exp(k*tp.postinj)*C),
                                data = .,
                                start = list(k = unique(.$epi.561.int.fit.lin.k),
                                             C = unique(.$epi.561.int.fit.lin.d)
                                )
                            )
                    )$parameters[[1]]
                )
     ),
     epi.561.int.fit.nls.d =
     as.numeric(try(summary(nls(formula = epi.561.int ~ (exp(k*tp.postinj)*C),
                                data = .,
                                start = list(k = unique(.$epi.561.int.fit.lin.k),
                                             C = unique(.$epi.561.int.fit.lin.d)
                                )
                            )
                    )$parameters[[2]]
                )
     )
  )

exp.fit.561.nls[[2]] <- as.numeric(exp.fit.561.nls[[2]])
exp.fit.561.nls[[3]] <- as.numeric(exp.fit.561.nls[[3]])

# Now we join these newly found parameters into the alldata table --------------

alldata <- alldata %>%
  full_join(exp.fit.561.nls, by = "cell.unique")

################################################################################
# 640 Fitting procedure
################################################################################

# Logarithmic linearization and linear regression ------------------------------

exp.fit.640.lin <- alldata %>%
  filter(is.na(epi.640.int) == FALSE) %>%
  group_by(cell.unique) %>%
  do(epi.640.int.fit.lin.d = 
     exp(as.numeric(lm(log(epi.640.int) ~ tp.postinj, data = .)[[1]][[1]])),
     epi.640.int.fit.lin.k =
     as.numeric(lm(log(epi.640.int) ~ tp.postinj, data = .)[[1]][[2]])
  )

exp.fit.640.lin[[2]] <- as.numeric(exp.fit.640.lin[[2]])
exp.fit.640.lin[[3]] <- as.numeric(exp.fit.640.lin[[3]])

# Now we join these newly found parameters into the alldata table --------------

alldata <- alldata %>%
  full_join(exp.fit.640.lin, by = "cell.unique")

# They serves as starting points for NLS procedure -----------------------------

exp.fit.640.nls <- alldata %>%
  filter(is.na(epi.640.int) == F) %>%
  filter(is.na(epi.640.int.fit.lin.d) == F) %>%
  filter(is.na(epi.640.int.fit.lin.k) == F) %>%
  group_by(cell.unique) %>%
  do(epi.640.int.fit.nls.k = 
     as.numeric(try(summary(nls(formula = epi.640.int ~ (exp(k*tp.postinj)*C),
                                data = .,
                                start = list(k = unique(.$epi.640.int.fit.lin.k),
                                             C = unique(.$epi.640.int.fit.lin.d)
                                )
                            )
                    )$parameters[[1]]
                )
     ),
     epi.640.int.fit.nls.d =
     as.numeric(try(summary(nls(formula = epi.640.int ~ (exp(k*tp.postinj)*C),
                                data = .,
                                start = list(k = unique(.$epi.640.int.fit.lin.k),
                                             C = unique(.$epi.640.int.fit.lin.d)
                                )
                            )
                    )$parameters[[2]]
                )
     )
  )

exp.fit.640.nls[[2]] <- as.numeric(exp.fit.640.nls[[2]])
exp.fit.640.nls[[3]] <- as.numeric(exp.fit.640.nls[[3]])

# Now we join these newly found parameters into the alldata table --------------

alldata <- alldata %>%
  full_join(exp.fit.640.nls, by = "cell.unique")

# Change data type of wrongly imported variables -------------------------------

alldata <- alldata %>%
  mutate_at(vars(matches("*405|488|561|640*")),
            as.double
  )

################################################################################
# Create new column combining LIN and NLS, add LIN if NLS did not work
################################################################################

alldata <- alldata %>% 
  mutate(epi.405.int.nls.lin.k = if_else(is.na(epi.405.int.fit.nls.k),
                                         epi.405.int.fit.lin.k,
                                         epi.405.int.fit.nls.k
         )
  ) %>%
  mutate(epi.488.int.nls.lin.k = if_else(is.na(epi.488.int.fit.nls.k),
                                         epi.488.int.fit.lin.k,
                                         epi.488.int.fit.nls.k
         )
  ) %>% 
  mutate(epi.561.int.nls.lin.k = if_else(is.na(epi.561.int.fit.nls.k),
                                         epi.561.int.fit.lin.k,
                                         epi.561.int.fit.nls.k
         )
  ) %>% 
  mutate(epi.640.int.nls.lin.k = if_else(is.na(epi.640.int.fit.nls.k),
                                         epi.640.int.fit.lin.k,
                                         epi.640.int.fit.nls.k
         )
  ) %>%
  mutate(epi.405.int.nls.lin.d = if_else(is.na(epi.405.int.fit.nls.d),
                                         epi.405.int.fit.lin.d,
                                         epi.405.int.fit.nls.d
         )
  ) %>%
  mutate(epi.488.int.nls.lin.d = if_else(is.na(epi.488.int.fit.nls.d),
                                         epi.488.int.fit.lin.d,
                                         epi.488.int.fit.nls.d
         )
  ) %>% 
  mutate(epi.561.int.nls.lin.d = if_else(is.na(epi.561.int.fit.nls.d),
                                         epi.561.int.fit.lin.d,
                                         epi.561.int.fit.nls.d
         )
  ) %>% 
  mutate(epi.640.int.nls.lin.d = if_else(is.na(epi.640.int.fit.nls.d),
                                         epi.640.int.fit.lin.d,
                                         epi.640.int.fit.nls.d
         )
  )

################################################################################
# Single cell standard deviation (SD) between fit and raw data.
# The difference between each datapoint and its fit is calculated. This is later
# used to compare goodness of fit and check, e.g. whether there is a correlation
# with SD and injection volume.
################################################################################

# Calculate difference between fit and raw data for each channel ---------------

alldata <- alldata %>%
  mutate(epi.405.int.fit.dif =
         ifelse(is.na(epi.405.int.nls.lin.k),
         NA,
         epi.405.int -
         (epi.405.int.nls.lin.d * exp(tp.postinj * epi.405.int.nls.lin.k))
         )
  ) %>%
  mutate(epi.488.int.fit.dif =
         ifelse(is.na(epi.488.int.nls.lin.k),
         NA,
         epi.488.int -
         (epi.488.int.nls.lin.d * exp(tp.postinj * epi.488.int.nls.lin.k))
         )
  ) %>%
  mutate(epi.561.int.fit.dif =
         ifelse(is.na(epi.561.int.nls.lin.k),
         NA,
         epi.561.int -
         (epi.561.int.nls.lin.d * exp(tp.postinj * epi.561.int.nls.lin.k))
         )
  ) %>%
  mutate(epi.640.int.fit.dif =
         ifelse(is.na(epi.640.int.nls.lin.k),
         NA,
         epi.640.int -
         (epi.640.int.nls.lin.d * exp(tp.postinj * epi.640.int.nls.lin.k))
         )
  )

# Calculate standard deviation -------------------------------------------------
  
stdev <- alldata %>%
  select(cell.unique,
         matches("^epi.(405|488|561|640).int.fit.dif$")
  ) %>%
  distinct() %>%
  droplevels() %>%
  group_by(cell.unique) %>%
  summarise(sd.fit.405 = sd(epi.405.int.fit.dif, na.rm = T),
            sd.fit.488 = sd(epi.488.int.fit.dif, na.rm = T),
            sd.fit.561 = sd(epi.561.int.fit.dif, na.rm = T),
            sd.fit.640 = sd(epi.640.int.fit.dif, na.rm = T)
  )

# Now we join these newly found standard deviations into the alldata table -----

alldata <- alldata %>%
  full_join(stdev, by = "cell.unique")

################################################################################
# SINGLE-ANALYTE PER EXPERIMENT
# Outlier removal
# mean
# standard deviation (sd)
# standard error of the mean (sem)
# Cell number after outlier removal and total cell number
################################################################################

# Calculate IRQ for all coeffiecients analyzed ---------------------------------

allcoef.iqr <- alldata %>%
  ungroup() %>% 
  select(exp.analyte.unique,
         matches("^epi.(405|488|561|640).int.nls.lin.(d|k)$")
  ) %>%
  distinct() %>%
  droplevels() %>%
  group_by(exp.analyte.unique) %>%
  summarise(k.405.IQR = IQR(epi.405.int.nls.lin.k, na.rm = T),             # IQR
            d.405.IQR = IQR(epi.405.int.nls.lin.d, na.rm = T),
            k.488.IQR = IQR(epi.488.int.nls.lin.k, na.rm = T),
            d.488.IQR = IQR(epi.488.int.nls.lin.d, na.rm = T),
            k.561.IQR = IQR(epi.561.int.nls.lin.k, na.rm = T),
            d.561.IQR = IQR(epi.561.int.nls.lin.d, na.rm = T),
            k.640.IQR = IQR(epi.640.int.nls.lin.k, na.rm = T),
            d.640.IQR = IQR(epi.640.int.nls.lin.d, na.rm = T),

            k.405.q1 = quantile(epi.405.int.nls.lin.k, na.rm = T)[[2]],      #q1
            d.405.q1 = quantile(epi.405.int.nls.lin.d, na.rm = T)[[2]],
            k.488.q1 = quantile(epi.488.int.nls.lin.k, na.rm = T)[[2]],
            d.488.q1 = quantile(epi.488.int.nls.lin.d, na.rm = T)[[2]],
            k.561.q1 = quantile(epi.561.int.nls.lin.k, na.rm = T)[[2]],
            d.561.q1 = quantile(epi.561.int.nls.lin.d, na.rm = T)[[2]],
            k.640.q1 = quantile(epi.640.int.nls.lin.k, na.rm = T)[[2]],
            d.640.q1 = quantile(epi.640.int.nls.lin.d, na.rm = T)[[2]],
  
            k.405.q2 = quantile(epi.405.int.nls.lin.k, na.rm = T)[[4]],      #q2
            d.405.q2 = quantile(epi.405.int.nls.lin.d, na.rm = T)[[4]],
            k.488.q2 = quantile(epi.488.int.nls.lin.k, na.rm = T)[[4]],
            d.488.q2 = quantile(epi.488.int.nls.lin.d, na.rm = T)[[4]],
            k.561.q2 = quantile(epi.561.int.nls.lin.k, na.rm = T)[[4]],
            d.561.q2 = quantile(epi.561.int.nls.lin.d, na.rm = T)[[4]],
            k.640.q2 = quantile(epi.640.int.nls.lin.k, na.rm = T)[[4]],
            d.640.q2 = quantile(epi.640.int.nls.lin.d, na.rm = T)[[4]]
  )
 
alldata <- alldata %>%
  full_join(allcoef.iqr, by = "exp.analyte.unique")
 
# Actual removal of outliers ---------------------------------------------------

alldata <- alldata %>%
  group_by(cell.unique) %>% 
  mutate(epi.405.int.nls.lin.k.out = if_else(epi.405.int.nls.lin.k >
                                             (k.405.q1 - 1.5*k.405.IQR) &
                                             epi.405.int.nls.lin.k <
                                             (k.405.q2 + 1.5*k.405.IQR),
                                             epi.405.int.nls.lin.k,
                                             as.double(NA)
                                      )
  ) %>% 
  mutate(epi.405.int.nls.lin.d.out = if_else(epi.405.int.nls.lin.d >
                                             (d.405.q1 - 1.5*d.405.IQR) &
                                             epi.405.int.nls.lin.d <
                                             (d.405.q2 + 1.5*d.405.IQR),
                                             epi.405.int.nls.lin.d,
                                             as.double(NA)
                                      )
  ) %>%
  mutate(epi.488.int.nls.lin.k.out = if_else(epi.488.int.nls.lin.k >
                                             (k.488.q1 - 1.5*k.488.IQR) &
                                             epi.488.int.nls.lin.k <
                                             (k.488.q2 + 1.5*k.488.IQR),
                                             epi.488.int.nls.lin.k,
                                             as.double(NA)
                                      )
  ) %>% 
  mutate(epi.488.int.nls.lin.d.out = if_else(epi.488.int.nls.lin.d >
                                             (d.488.q1 - 1.5*d.488.IQR) &
                                             epi.488.int.nls.lin.d <
                                             (d.488.q2 + 1.5*d.488.IQR),
                                             epi.488.int.nls.lin.d,
                                             as.double(NA)
                                      )
  ) %>%
  mutate(epi.561.int.nls.lin.k.out = if_else(epi.561.int.nls.lin.k >
                                             (k.561.q1 - 1.5*k.561.IQR) &
                                             epi.561.int.nls.lin.k <
                                             (k.561.q2 + 1.5*k.561.IQR),
                                             epi.561.int.nls.lin.k,
                                             as.double(NA)
                                      )
  ) %>% 
  mutate(epi.561.int.nls.lin.d.out = if_else(epi.561.int.nls.lin.d >
                                             (d.561.q1 - 1.5*d.561.IQR) &
                                             epi.561.int.nls.lin.d <
                                             (d.561.q2 + 1.5*d.561.IQR),
                                             epi.561.int.nls.lin.d,
                                             as.double(NA)
                                      )
  ) %>%
  mutate(epi.640.int.nls.lin.k.out = if_else(epi.640.int.nls.lin.k >
                                             (k.640.q1 - 1.5*k.640.IQR) &
                                             epi.640.int.nls.lin.k <
                                             (k.640.q2 + 1.5*k.640.IQR),
                                             epi.640.int.nls.lin.k,
                                             as.double(NA)
                                      )
  ) %>% 
  mutate(epi.640.int.nls.lin.d.out = if_else(epi.640.int.nls.lin.d >
                                             (d.640.q1 - 1.5*d.640.IQR) &
                                             epi.640.int.nls.lin.d <
                                             (d.640.q2 + 1.5*d.640.IQR),
                                             epi.640.int.nls.lin.d,
                                             as.double(NA)
                                      )
  )

# Calculate mean, sd, sem, nr. per channel and total cell nr. ------------------

allcoef.meansd <- alldata %>%
  ungroup() %>%
  select(exp.analyte.unique,
         matches("^epi.(405|488|561|640).int.nls.lin.(d|k).out$")
  ) %>%
  distinct() %>% 
  droplevels() %>%
  arrange() %>% 
  group_by(exp.analyte.unique) %>% 
  summarise(epi.405.int.nls.lin.k.out.mean =                              # mean
            mean(epi.405.int.nls.lin.k.out, na.rm = T),
            epi.405.int.nls.lin.d.out.mean =
            mean(epi.405.int.nls.lin.d.out, na.rm = T),
            epi.488.int.nls.lin.k.out.mean =
            mean(epi.488.int.nls.lin.k.out, na.rm = T),
            epi.488.int.nls.lin.d.out.mean =
            mean(epi.488.int.nls.lin.d.out, na.rm = T),
            epi.561.int.nls.lin.k.out.mean =
            mean(epi.561.int.nls.lin.k.out, na.rm = T),
            epi.561.int.nls.lin.d.out.mean =
            mean(epi.561.int.nls.lin.d.out, na.rm = T),
            epi.640.int.nls.lin.k.out.mean =
            mean(epi.640.int.nls.lin.k.out, na.rm = T),
            epi.640.int.nls.lin.d.out.mean =
            mean(epi.640.int.nls.lin.d.out, na.rm = T),

            epi.405.int.nls.lin.k.out.sd =                                  # sd
            sd(epi.405.int.nls.lin.k.out, na.rm = T),
            epi.405.int.nls.lin.d.out.sd =
            sd(epi.405.int.nls.lin.d.out, na.rm = T),
            epi.488.int.nls.lin.k.out.sd =
            sd(epi.488.int.nls.lin.k.out, na.rm = T),
            epi.488.int.nls.lin.d.out.sd =
            sd(epi.488.int.nls.lin.d.out, na.rm = T),
            epi.561.int.nls.lin.k.out.sd =
            sd(epi.561.int.nls.lin.k.out, na.rm = T),
            epi.561.int.nls.lin.d.out.sd =
            sd(epi.561.int.nls.lin.d.out, na.rm = T),
            epi.640.int.nls.lin.k.out.sd =
            sd(epi.640.int.nls.lin.k.out, na.rm = T),
            epi.640.int.nls.lin.d.out.sd =
            sd(epi.640.int.nls.lin.d.out, na.rm = T),
    
            epi.405.int.nls.lin.k.out.sem =                                # sem
              sd(epi.405.int.nls.lin.k.out, na.rm = T) /
              sqrt(sum(!is.na(epi.405.int.nls.lin.k.out))),
            epi.405.int.nls.lin.d.out.sem =
              sd(epi.405.int.nls.lin.d.out, na.rm = T) /
              sqrt(sum(!is.na(epi.405.int.nls.lin.d.out))),
            epi.488.int.nls.lin.k.out.sem =
              sd(epi.488.int.nls.lin.k.out, na.rm = T) /
              sqrt(sum(!is.na(epi.488.int.nls.lin.k.out))),
            epi.488.int.nls.lin.d.out.sem =
              sd(epi.488.int.nls.lin.d.out, na.rm = T) /
              sqrt(sum(!is.na(epi.488.int.nls.lin.d.out))),
            epi.561.int.nls.lin.k.out.sem =
              sd(epi.561.int.nls.lin.k.out, na.rm = T) /
              sqrt(sum(!is.na(epi.561.int.nls.lin.k.out))),
            epi.561.int.nls.lin.d.out.sem =
              sd(epi.561.int.nls.lin.d.out, na.rm = T) /
              sqrt(sum(!is.na(epi.561.int.nls.lin.d.out))),
            epi.640.int.nls.lin.k.out.sem =
              sd(epi.640.int.nls.lin.k.out, na.rm = T) /
              sqrt(sum(!is.na(epi.640.int.nls.lin.k.out))),
            epi.640.int.nls.lin.d.out.sem =
              sd(epi.640.int.nls.lin.d.out, na.rm = T) /
              sqrt(sum(!is.na(epi.640.int.nls.lin.d.out))),

            epi.405.int.nls.lin.k.out.n =                 # cell nr. per channel
            sum(!is.na(epi.405.int.nls.lin.k.out)),
            epi.488.int.nls.lin.k.out.n =
            sum(!is.na(epi.488.int.nls.lin.k.out)),
            epi.561.int.nls.lin.k.out.n =
            sum(!is.na(epi.561.int.nls.lin.k.out)),
            epi.640.int.nls.lin.k.out.n =
            sum(!is.na(epi.640.int.nls.lin.k.out)),
    
            ntot = length(epi.640.int.nls.lin.k.out)            # total cell nr.

  )

alldata <- alldata %>%
  full_join(allcoef.meansd, by = "exp.analyte.unique")

################################################################################
# SINGLE ANALYTE, SO MULTIPLE EXPERIMENTS COMBINED
# Outlier removal
# mean
# standard deviation (sd)
# standard error of the mean (sem)
# Cell number after outlier removal and total cell number
################################################################################

# Calculate IRQ for all coeffiecients analyzed ---------------------------------

allcoef.iqr.analyte.trans <- alldata %>%
  ungroup() %>% 
  select(analyte.trans,
         matches("^epi.(405|488|561|640).int.nls.lin.(d|k)$")
  ) %>%
  distinct() %>%
  droplevels() %>%
  group_by(analyte.trans) %>%
  summarise(k.405.IQR.trans = IQR(epi.405.int.nls.lin.k, na.rm = T),       # IQR
            d.405.IQR.trans = IQR(epi.405.int.nls.lin.d, na.rm = T),
            k.488.IQR.trans = IQR(epi.488.int.nls.lin.k, na.rm = T),
            d.488.IQR.trans = IQR(epi.488.int.nls.lin.d, na.rm = T),
            k.561.IQR.trans = IQR(epi.561.int.nls.lin.k, na.rm = T),
            d.561.IQR.trans = IQR(epi.561.int.nls.lin.d, na.rm = T),
            k.640.IQR.trans = IQR(epi.640.int.nls.lin.k, na.rm = T),
            d.640.IQR.trans = IQR(epi.640.int.nls.lin.d, na.rm = T),

            k.405.q1.trans = quantile(epi.405.int.nls.lin.k, na.rm = T)[[2]], #q1
            d.405.q1.trans = quantile(epi.405.int.nls.lin.d, na.rm = T)[[2]],
            k.488.q1.trans = quantile(epi.488.int.nls.lin.k, na.rm = T)[[2]],
            d.488.q1.trans = quantile(epi.488.int.nls.lin.d, na.rm = T)[[2]],
            k.561.q1.trans = quantile(epi.561.int.nls.lin.k, na.rm = T)[[2]],
            d.561.q1.trans = quantile(epi.561.int.nls.lin.d, na.rm = T)[[2]],
            k.640.q1.trans = quantile(epi.640.int.nls.lin.k, na.rm = T)[[2]],
            d.640.q1.trans = quantile(epi.640.int.nls.lin.d, na.rm = T)[[2]],

            k.405.q2.trans = quantile(epi.405.int.nls.lin.k, na.rm = T)[[4]], #q2
            d.405.q2.trans = quantile(epi.405.int.nls.lin.d, na.rm = T)[[4]],
            k.488.q2.trans = quantile(epi.488.int.nls.lin.k, na.rm = T)[[4]],
            d.488.q2.trans = quantile(epi.488.int.nls.lin.d, na.rm = T)[[4]],
            k.561.q2.trans = quantile(epi.561.int.nls.lin.k, na.rm = T)[[4]],
            d.561.q2.trans = quantile(epi.561.int.nls.lin.d, na.rm = T)[[4]],
            k.640.q2.trans = quantile(epi.640.int.nls.lin.k, na.rm = T)[[4]],
            d.640.q2.trans = quantile(epi.640.int.nls.lin.d, na.rm = T)[[4]]
  )
 
alldata <- alldata %>%
  full_join(allcoef.iqr.analyte.trans, by = "analyte.trans")
 
# Actual removal of outliers ---------------------------------------------------

alldata <- alldata %>%
  group_by(cell.unique) %>% 
  mutate(epi.405.int.nls.lin.k.out.trans =                        # Outliers 405
         if_else(epi.405.int.nls.lin.k >
                 (k.405.q1.trans - 1.5*k.405.IQR.trans) &
                 epi.405.int.nls.lin.k <
                 (k.405.q2.trans + 1.5*k.405.IQR.trans),
                 
                 epi.405.int.nls.lin.k,
                 as.double(NA)
         )
  ) %>%
  mutate(epi.405.int.nls.lin.d.out.trans = 
         if_else(epi.405.int.nls.lin.d >
                 (d.405.q1.trans - 1.5*d.405.IQR.trans) &
                 epi.405.int.nls.lin.d <
                 (d.405.q2.trans + 1.5*d.405.IQR.trans),
                 
                 epi.405.int.nls.lin.d,
                 as.double(NA)
         )
  ) %>%
  mutate(epi.488.int.nls.lin.k.out.trans =                        # Outliers 488
         if_else(epi.488.int.nls.lin.k >
                 (k.488.q1.trans - 1.5*k.488.IQR.trans) &
                 epi.488.int.nls.lin.k <
                 (k.488.q2.trans + 1.5*k.488.IQR.trans),
                 
                 epi.488.int.nls.lin.k,
                 as.double(NA)
         )
  ) %>%
  mutate(epi.488.int.nls.lin.d.out.trans = 
         if_else(epi.488.int.nls.lin.d >
                 (d.488.q1.trans - 1.5*d.488.IQR.trans) &
                 epi.488.int.nls.lin.d <
                 (d.488.q2.trans + 1.5*d.488.IQR.trans),
                 
                 epi.488.int.nls.lin.d,
                 as.double(NA)
         )
  ) %>%
  mutate(epi.561.int.nls.lin.k.out.trans =                        # Outliers 561
         if_else(epi.561.int.nls.lin.k >
                 (k.561.q1.trans - 1.5*k.561.IQR.trans) &
                 epi.561.int.nls.lin.k <
                 (k.561.q2.trans + 1.5*k.561.IQR.trans),
                 
                 epi.561.int.nls.lin.k,
                 as.double(NA)
         )
  ) %>%
  mutate(epi.561.int.nls.lin.d.out.trans = 
         if_else(epi.561.int.nls.lin.d >
                 (d.561.q1.trans - 1.5*d.561.IQR.trans) &
                 epi.561.int.nls.lin.d <
                 (d.561.q2.trans + 1.5*d.561.IQR.trans),
                 
                 epi.561.int.nls.lin.d,
                 as.double(NA)
         )
  ) %>%
  mutate(epi.640.int.nls.lin.k.out.trans =                        # Outliers 640
         if_else(epi.640.int.nls.lin.k >
                 (k.640.q1.trans - 1.5*k.640.IQR.trans) &
                 epi.640.int.nls.lin.k <
                 (k.640.q2.trans + 1.5*k.640.IQR.trans),
                 
                 epi.640.int.nls.lin.k,
                 as.double(NA)
         )
  ) %>%
  mutate(epi.640.int.nls.lin.d.out.trans = 
         if_else(epi.640.int.nls.lin.d >
                 (d.640.q1.trans - 1.5*d.640.IQR.trans) &
                 epi.640.int.nls.lin.d <
                 (d.640.q2.trans + 1.5*d.640.IQR.trans),
                 
                 epi.640.int.nls.lin.d,
                 as.double(NA)
         )
  )

################################################################################
# Single-analyte: coefficients are now calculated as mean and stdev
################################################################################

allcoef.meansd.trans <- alldata %>%
  ungroup() %>%
  select(analyte.trans,
         matches("^epi.(405|488|561|640).int.nls.lin.(d|k).out.trans$")
  ) %>%
  distinct() %>% 
  droplevels() %>%
  arrange() %>% 
  group_by(analyte.trans) %>% 
  summarise(epi.405.int.nls.lin.k.out.mean.trans =                        # mean
            mean(epi.405.int.nls.lin.k.out.trans, na.rm = T),
            epi.405.int.nls.lin.d.out.mean.trans =
            mean(epi.405.int.nls.lin.d.out.trans, na.rm = T),
            epi.488.int.nls.lin.k.out.mean.trans =
            mean(epi.488.int.nls.lin.k.out.trans, na.rm = T),
            epi.488.int.nls.lin.d.out.mean.trans =
            mean(epi.488.int.nls.lin.d.out.trans, na.rm = T),
            epi.561.int.nls.lin.k.out.mean.trans =
            mean(epi.561.int.nls.lin.k.out.trans, na.rm = T),
            epi.561.int.nls.lin.d.out.mean.trans =
            mean(epi.561.int.nls.lin.d.out.trans, na.rm = T),
            epi.640.int.nls.lin.k.out.mean.trans =
            mean(epi.640.int.nls.lin.k.out.trans, na.rm = T),
            epi.640.int.nls.lin.d.out.mean.trans =
            mean(epi.640.int.nls.lin.d.out.trans, na.rm = T),

            epi.405.int.nls.lin.k.out.sd.trans =                            # sd
            sd(epi.405.int.nls.lin.k.out.trans, na.rm = T),
            epi.405.int.nls.lin.d.out.sd.trans =
            sd(epi.405.int.nls.lin.d.out.trans, na.rm = T),
            epi.488.int.nls.lin.k.out.sd.trans =
            sd(epi.488.int.nls.lin.k.out.trans, na.rm = T),
            epi.488.int.nls.lin.d.out.sd.trans =
            sd(epi.488.int.nls.lin.d.out.trans, na.rm = T),
            epi.561.int.nls.lin.k.out.sd.trans =
            sd(epi.561.int.nls.lin.k.out.trans, na.rm = T),
            epi.561.int.nls.lin.d.out.sd.trans =
            sd(epi.561.int.nls.lin.d.out.trans, na.rm = T),
            epi.640.int.nls.lin.k.out.sd.trans =
            sd(epi.640.int.nls.lin.k.out.trans, na.rm = T),
            epi.640.int.nls.lin.d.out.sd.trans =
            sd(epi.640.int.nls.lin.d.out.trans, na.rm = T),

            epi.405.int.nls.lin.k.out.sem.trans =                          # sem
              sd(epi.405.int.nls.lin.k.out.trans, na.rm = T) /
              sqrt(sum(!is.na(epi.405.int.nls.lin.k.out.trans))),
            epi.405.int.nls.lin.d.out.sem.trans =
              sd(epi.405.int.nls.lin.d.out.trans, na.rm = T) /
              sqrt(sum(!is.na(epi.405.int.nls.lin.d.out.trans))),
            epi.488.int.nls.lin.k.out.sem.trans =
              sd(epi.488.int.nls.lin.k.out.trans, na.rm = T) /
              sqrt(sum(!is.na(epi.488.int.nls.lin.k.out.trans))),
            epi.488.int.nls.lin.d.out.sem.trans =
              sd(epi.488.int.nls.lin.d.out.trans, na.rm = T) /
              sqrt(sum(!is.na(epi.488.int.nls.lin.d.out.trans))),
            epi.561.int.nls.lin.k.out.sem.trans =
              sd(epi.561.int.nls.lin.k.out.trans, na.rm = T) /
              sqrt(sum(!is.na(epi.561.int.nls.lin.k.out.trans))),
            epi.561.int.nls.lin.d.out.sem.trans =
              sd(epi.561.int.nls.lin.d.out.trans, na.rm = T) /
              sqrt(sum(!is.na(epi.561.int.nls.lin.d.out.trans))),
            epi.640.int.nls.lin.k.out.sem.trans =
              sd(epi.640.int.nls.lin.k.out.trans, na.rm = T) /
              sqrt(sum(!is.na(epi.640.int.nls.lin.k.out.trans))),
            epi.640.int.nls.lin.d.out.sem.trans =
              sd(epi.640.int.nls.lin.d.out.trans, na.rm = T) /
              sqrt(sum(!is.na(epi.640.int.nls.lin.d.out.trans))),

            epi.405.int.nls.lin.k.out.n.trans =           # Cell nr. per channel
            sum(!is.na(epi.405.int.nls.lin.k.out.trans)),
            epi.488.int.nls.lin.k.out.n.trans =
            sum(!is.na(epi.488.int.nls.lin.k.out.trans)),
            epi.561.int.nls.lin.k.out.n.trans =
            sum(!is.na(epi.561.int.nls.lin.k.out.trans)),
            epi.640.int.nls.lin.k.out.n.trans =
            sum(!is.na(epi.640.int.nls.lin.k.out.trans)),

            ntot.trans =                                        # Total cell nr.
            length(epi.640.int.nls.lin.k.out.trans)

  )

alldata <- alldata %>%
  full_join(allcoef.meansd.trans, by = "analyte.trans")

################################################################################
# DATA EXPORT
################################################################################

setwd("../02_01_single_cell_fit_oldschool")

write_csv(alldata,
	         paste(getwd(),
	         	     "/",
	         	     "single_cell_fit.csv",
	         	     sep = ""
	         	),
	          col_names = TRUE
)

tempdata <- alldata %>%
  ungroup() %>% 
  select(c("ntot.trans",
           "analyte.trans",
           "epi.405.int.nls.lin.k.out.mean.trans",
           "epi.488.int.nls.lin.k.out.mean.trans",
           "epi.561.int.nls.lin.k.out.mean.trans",
           "epi.640.int.nls.lin.k.out.mean.trans",
         )
  ) %>%
  distinct() %>%
  droplevels() %>% 
  mutate("halflife405_min" = log(2)/epi.405.int.nls.lin.k.out.mean.trans*-60,
         "halflife488_min" = log(2)/epi.488.int.nls.lin.k.out.mean.trans*-60,
         "halflife561_min" = log(2)/epi.561.int.nls.lin.k.out.mean.trans*-60,
         "halflife640_min" = log(2)/epi.640.int.nls.lin.k.out.mean.trans*-60
  )
  

write_csv(tempdata,
          paste(getwd(),
                "/",
                "degradation_rates_each_analyte.csv",
                sep = ""
          ),
          col_names = TRUE
)