
#################   damagelogo six   ##########################

topic_clus <- get(load("../../ancient-damage/utilities/Lazaridis_moderns/clus_2/model.rda"))
omega <- topic_clus$omega
theta <- topic_clus$theta

sig_names = NULL
ic.scale=TRUE
max_pos = 20
flanking_bases=1
yscale_change = TRUE
xaxis=TRUE
yaxis=TRUE
xlab = " "
xaxis_fontsize=10
xlab_fontsize=20
title_aligner = 15
y_fontsize=20
title_fontsize = 20
mut_width=2
start=0.0001
renyi_alpha = 1
inflation_factor = c(2,1,2)
base_probs_list = NULL
pop_names=paste0("Cluster ",1:dim(theta)[2])
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



