
###############  Test aRchaic_cluster() - remove strand info   ###################

labs <- c(rep("Skoglund",5), rep("moderns-50",50))

folders = c("../docs/data/Skoglund/", "../docs/data/moderns-50/")
K = 2
tol = 1
labs = labs
run_from = "gom"
output_dir = "../docs/output/skoglund_moderns_2/"
run_index = NULL
breaks = c(-1, seq(1,20,1))
flanking_bases = 1
gom_method = "independent"
topic_cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
               "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
               "brown4","darkorchid","magenta","yellow", "azure1","azure4")
structure.control = list()
logo.control = list()
topics.control = list()
positive_logo = FALSE
# save_plot = TRUE,
structure_width = 5
structure_height = 8


clus <- aRchaic_cluster(folders = c("../data/Skoglund/", "../data/moderns-50/"),
                        K = 2, 
                        tol = 1,
                        labs = labs,
                        run_from = "gom",
                        output_dir = "../output/skoglund_moderns_2/")