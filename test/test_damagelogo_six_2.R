
###################   damagelogo 6 test    ####################


topic_clus <- get(load("../../ancient-damage/utilities/fig1_rep/Cluster3/model.rda"))
theta_pool <- topic_clus$theta

base_probs_list <- get(load("../../ancient-damage/utilities/base_probs_temp_2.rda"))

sig_names = NULL
max_pos = 20
flanking_bases=1
ranges = c(0.5,1)
mutlogo.control = list()
breaklogo.control = list()
base_probs_list = NULL
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
output_width = 1200
output_height = 700
