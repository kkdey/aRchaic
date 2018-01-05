
#################   damagelogo six   ##########################

topic_clus <- get(load("../../ancient-damage/utilities/moderns_Pinhasi/clus_2/model.rda"))
omega <- topic_clus$omega
theta <- topic_clus$theta

base_probs_list <- get(load("../../ancient-damage/utilities/base_probs_temp_2.rda"))
# base_probs_mat <- matrix(NA, 10, 3)
# rownames(base_probs_mat) <- c("A", "C", "G", "T", "C->A", "C->G", "C->T", "T->A", "T->C", "T->G")
# for(l in 1:length(base_probs_list)){
#   base_probs_mat[match(names(base_probs_list[[l]]), rownames(base_probs_mat)),l] <- base_probs_list[[l]]
# }
#



damageLogo_six(theta,
               base_probs_list = base_probs_list,
               ranges = c(1,3),
               logoport_x = 0.7,
               logoport_y= 0.5,
               logoport_width= 1.2,
               logoport_height= 1.2,
               breaklogoport_x = 0.6,
               breaklogoport_y = 0.45,
               breaklogoport_width=0.4,
               breaklogoport_height=1,
               lineport_x = 0.8,
               lineport_y=0.5,
               lineport_width=1,
               lineport_height=1.4)


theta_pool <- theta
sig_names = NULL
max_pos = 20
flanking_bases=1
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






theta_pool <- theta
tab <- prop_patterns_list[[1]]
colnames(tab) <- c("left \n flank", "mismatch", "right \n flank")

library(Logolas)
cols = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category ==
                                       'qual',]
col_vector = unlist(mapply(RColorBrewer::brewer.pal, cols$maxcolors, rownames(cols)))
total_chars = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O",
                "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "zero", "one", "two",
                "three", "four", "five", "six", "seven", "eight", "nine", "dot", "comma",
                "dash", "colon", "semicolon", "leftarrow", "rightarrow")

set.seed(20)
color_profile_1 <- list("type" = "per_symbol",
                      "col" = sample(col_vector, length(total_chars), replace=FALSE))

rownames(tab)[match(c("C->A", "C->G", "C->T",
        "T->A", "T->C", "T->G"), rownames(tab))] <- c("C>A", "C>G", "C>T",
                                                      "T>A", "T>C", "T>G")

breaks_theta_vec = breaks_theta[,l, drop=FALSE]
colnames(breaks_theta_vec) <- c("strand break")


color_profile = list("type" = "per_row",
                     "col" = RColorBrewer::brewer.pal(4,name ="Spectral"))

prob_mutation[l,]

pos_data <- data.frame(position = as.numeric(names(prob_mutation[l,])),
                       val = as.numeric(prob_mutation[l,]))

get_viewport_logo(1, 3)

seekViewport(paste0("plotlogo", 1))
vp1 <- viewport(x=0.7, y=0.5, width=1.2, height=1.1)
pushViewport(vp1)
nlogomaker(tab,
           logoheight = 'log',
           color_profile = color_profile_1,
           frame_width = 1,
           xlab = "",
           main_fontsize = 14,
           xlab_fontsize=13,
           pop_name = 'Mismatch and \n flanking base composition',
           control = list(epsilon=0.25,gap_ylab=3.5, gap_xlab = 4,
                          round_off = 1, posbins = 3, negbins = 2),
           newpage = FALSE
)
upViewport(0)


seekViewport(paste0("plotlogo", 2))
vp2 <- viewport(x=0.5, y=0.4, width=0.7, height=1)
pushViewport(vp2)
nlogomaker(breaks_theta_vec,logoheight = 'ic_log',color_profile = color_profile,
           frame_width = 1, main_fontsize = 14,
           control = list(gap_ylab=3.5, epsilon = 0.01, round_off = 0, symm = TRUE),
           xaxis = FALSE, yaxis = FALSE, col_line_split = "white",
           newpage = FALSE, xlab = "",
           pop_name = "strand break \n composition")
upViewport(0)

seekViewport(paste0("plotlogo", 3))

vp3 = viewport(x = 0.4, y = 0.5, width=1, height=1)
p <- ggplot(data=pos_data, aes(x=position,y=val)) + geom_point() +
  ggtitle("Probability of mismatch along the read") +
  labs(x="position from end of read",y="probability of mismatch") +
  scale_x_continuous(limits = c(0, 20))
print(p, vp = vp3)


set.seed(20)
cols = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category ==
                                       'qual',]
col_vector = unlist(mapply(RColorBrewer::brewer.pal, cols$maxcolors, rownames(cols)))
total_chars = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O",
                "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "zero", "one", "two",
                "three", "four", "five", "six", "seven", "eight", "nine", "dot", "comma",
                "dash", "colon", "semicolon", "leftarrow", "rightarrow")
color_profile <- list("type" = "per_symbol",
                      "col" = sample(col_vector, length(total_chars), replace=FALSE))

grid.newpage()
nlogomaker(pwm1,
           logoheight = 'log',
           color_profile = color_profile,
           frame_width = 1,
           xlab = "Position",
           bg = base_probs_mat,
           pop_name = '',
           control = list(epsilon=0.25,gap_ylab=3.5,
                          round_off = 1),
           xaxis = FALSE
)

