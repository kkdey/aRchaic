

library(aRchaic)
library(Logolas)
library(ggplot2)
base_probs_list <- get(load("../../ancient-damage/utilities/base_probs_temp_2.rda"))

aRchaic_view(file = "../data/moderns/HG00096.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.q30.csv",
             output_dir = "../utilities/", title = "GBR1",
             logo.control = list(mut_ranges = c(1,1.5),
                                 break_ranges = c(0.5,0.5),
                                 base_probs_list = base_probs_list))

aRchaic_view(file = "../data/Pinhasi/IR1.q30.csv",
             output_dir = "../utilities/", title = "IR1",
             logo.control = list(mut_ranges = c(1,1.5),
                                 break_ranges = c(0.5,0.5),
                                 base_probs_list = base_probs_list))
