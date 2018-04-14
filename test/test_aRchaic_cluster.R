
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



theta_pool <- topic_clus$theta
max_pos = 20
flanking_bases=1
mutlogo.control = list()
breaklogo.control = list()
base_probs_list = NULL
clip = 0
mut_ranges  = c(0, 0)
break_ranges = c(0, 0)
logoport_x = 0.7
logoport_y= 0.5
logoport_width= 1.2
logoport_height= 1.1
breaklogoport_x = 0.5
breaklogoport_y = 0.4
breaklogoport_width=0.7
breaklogoport_height=1
lineport_x = 0.4
lineport_y=0.5
lineport_width=1
lineport_height=1
output_dir = NULL
filename = NULL
output_width = 1200
output_height = 700


labs <- c(rep("Skoglund",5), rep("moderns-50",50))
clus <- aRchaic_cluster(folders = c("../docs/data/Skoglund/", 
                                    "../docs/data/moderns-50/"),
                        K = 2, 
                        tol = 1,
                        labs = labs,
                        run_from = "gom",
                        output_dir = "../docs/output/skoglund_moderns_2/")




labs <- c(rep("Skoglund",5), rep("moderns-50",50))
base_probs_list <- get(load("../output/base_probs_temp_2.rda"))

clus <- aRchaic_cluster(folders = c("../docs/data/Skoglund/", "../docs/data/moderns-50/"),
                        K = 2, 
                        tol = 1,
                        labs = labs,
                        run_from = "gom",
                        logo.control = list(base_probs_list = base_probs_list,
                                            mut_ranges = c(1,1.5),
                                            break_ranges = c(0.5,0.5)),
                        output_dir = "../docs/output/skoglund_moderns_2_background/")


