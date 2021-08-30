#!/usr/bin/env Rscript

library(ggplot2)
library(ggpubr)



lower_ratio = 1.8
upper_ratio = 2.3


p1_file <- '/Users/lding/Documents/X-SCID/data/generate_features_and_outputs_v2/ML_features/P1_all_data_combine_features.txt'
p1 <- read.table(p1_file, sep='\t', header=TRUE, row.names = 1)
p1$log_output <- log2(p1$output)

pl1 <- ggplot(p1, aes(x=log_output)) + geom_density()
# Add mean line
pl1 <- pl1 + geom_vline(aes(xintercept=mean(log_output)),
                        color="blue", linetype="dashed", size=1)
# Add mean-sd, mean+sd line
pl1 <- pl1 + geom_vline(aes(xintercept=mean(log_output)-lower_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl1 <- pl1 + geom_vline(aes(xintercept=mean(log_output)+upper_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl1



p2_file <- '/Users/lding/Documents/X-SCID/data/generate_features_and_outputs_v2/ML_features/P2_all_data_combine_features.txt'
p2 <- read.table(p2_file, sep='\t', header=TRUE, row.names = 1)
p2$log_output <- log2(p2$output)

pl2 <- ggplot(p2, aes(x=log_output)) + geom_density()
# Add mean line
pl2 <- pl2 + geom_vline(aes(xintercept=mean(log_output)),
                        color="blue", linetype="dashed", size=1)
# Add mean-sd, mean+sd line
pl2 <- pl2 + geom_vline(aes(xintercept=mean(log_output)-lower_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl2 <- pl2 + geom_vline(aes(xintercept=mean(log_output)+upper_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl2



p3_file <- '/Users/lding/Documents/X-SCID/data/generate_features_and_outputs_v2/ML_features/P3_all_data_combine_features.txt'
p3 <- read.table(p3_file, sep='\t', header=TRUE, row.names = 1)
p3$log_output <- log2(p3$output)

pl3 <- ggplot(p3, aes(x=log_output)) + geom_density()
# Add mean line
pl3 <- pl3 + geom_vline(aes(xintercept=mean(log_output)),
                        color="blue", linetype="dashed", size=1)
# Add mean-sd, mean+sd line
pl3 <- pl3 + geom_vline(aes(xintercept=mean(log_output)-lower_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl3 <- pl3 + geom_vline(aes(xintercept=mean(log_output)+upper_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl3



p4_file <- '/Users/lding/Documents/X-SCID/data/generate_features_and_outputs_v2/ML_features/P4_all_data_combine_features.txt'
p4 <- read.table(p4_file, sep='\t', header=TRUE, row.names = 1)
p4$log_output <- log2(p4$output)

pl4 <- ggplot(p4, aes(x=log_output)) + geom_density()
# Add mean line
pl4 <- pl4 + geom_vline(aes(xintercept=mean(log_output)),
                        color="blue", linetype="dashed", size=1)
# Add mean-sd, mean+sd line
pl4 <- pl4 + geom_vline(aes(xintercept=mean(log_output)-lower_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl4 <- pl4 + geom_vline(aes(xintercept=mean(log_output)+upper_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl4



p5_file <- '/Users/lding/Documents/X-SCID/data/generate_features_and_outputs_v2/ML_features/P5_all_data_combine_features.txt'
p5 <- read.table(p5_file, sep='\t', header=TRUE, row.names = 1)
p5$log_output <- log2(p5$output)

pl5 <- ggplot(p5, aes(x=log_output)) + geom_density()
# Add mean line
pl5 <- pl5 + geom_vline(aes(xintercept=mean(log_output)),
                        color="blue", linetype="dashed", size=1)
# Add mean-sd, mean+sd line
pl5 <- pl5 + geom_vline(aes(xintercept=mean(log_output)-lower_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl5 <- pl5 + geom_vline(aes(xintercept=mean(log_output)+upper_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl5



p6_file <- '/Users/lding/Documents/X-SCID/data/generate_features_and_outputs_v2/ML_features/P6_all_data_combine_features.txt'
p6 <- read.table(p6_file, sep='\t', header=TRUE, row.names = 1)
p6$log_output <- log2(p6$output)

pl6 <- ggplot(p6, aes(x=log_output)) + geom_density()
# Add mean line
pl6 <- pl6 + geom_vline(aes(xintercept=mean(log_output)),
              color="blue", linetype="dashed", size=1)
# Add mean-sd, mean+sd line
pl6 <- pl6 + geom_vline(aes(xintercept=mean(log_output)-lower_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl6 <- pl6 + geom_vline(aes(xintercept=mean(log_output)+upper_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl6



p7_file <- '/Users/lding/Documents/X-SCID/data/generate_features_and_outputs_v2/ML_features/P7_all_data_combine_features.txt'
p7 <- read.table(p7_file, sep='\t', header=TRUE, row.names = 1)
p7$log_output <- log2(p7$output)

pl7 <- ggplot(p7, aes(x=log_output)) + geom_density()
# Add mean line
pl7 <- pl7 + geom_vline(aes(xintercept=mean(log_output)),
                        color="blue", linetype="dashed", size=1)
# Add mean-sd, mean+sd line
pl7 <- pl7 + geom_vline(aes(xintercept=mean(log_output)-lower_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl7 <- pl7 + geom_vline(aes(xintercept=mean(log_output)+upper_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl7



p8_file <- '/Users/lding/Documents/X-SCID/data/generate_features_and_outputs_v2/ML_features/P8_all_data_combine_features.txt'
p8 <- read.table(p8_file, sep='\t', header=TRUE, row.names = 1)
p8$log_output <- log2(p8$output)

pl8 <- ggplot(p8, aes(x=log_output)) + geom_density()
# Add mean line
pl8 <- pl8 + geom_vline(aes(xintercept=mean(log_output)),
                        color="blue", linetype="dashed", size=1)
# Add mean-sd, mean+sd line
pl8 <- pl8 + geom_vline(aes(xintercept=mean(log_output)-lower_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl8 <- pl8 + geom_vline(aes(xintercept=mean(log_output)+upper_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl8



p9_file <- '/Users/lding/Documents/X-SCID/data/generate_features_and_outputs_v2/ML_features/P9_all_data_combine_features.txt'
p9 <- read.table(p9_file, sep='\t', header=TRUE, row.names = 1)
p9$log_output <- log2(p9$output)

pl9 <- ggplot(p9, aes(x=log_output)) + geom_density()
# Add mean line
pl9 <- pl9 + geom_vline(aes(xintercept=mean(log_output)),
                        color="blue", linetype="dashed", size=1)
# Add mean-sd, mean+sd line
pl9 <- pl9 + geom_vline(aes(xintercept=mean(log_output)-lower_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl9 <- pl9 + geom_vline(aes(xintercept=mean(log_output)+upper_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl9



p10_file <- '/Users/lding/Documents/X-SCID/data/generate_features_and_outputs_v2/ML_features/P10_all_data_combine_features.txt'
p10 <- read.table(p10_file, sep='\t', header=TRUE, row.names = 1)
p10$log_output <- log2(p10$output)

pl10 <- ggplot(p10, aes(x=log_output)) + geom_density()
# Add mean line
pl10 <- pl10 + geom_vline(aes(xintercept=mean(log_output)),
                        color="blue", linetype="dashed", size=1)
# Add mean-sd, mean+sd line
pl10 <- pl10 + geom_vline(aes(xintercept=mean(log_output)-lower_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl10 <- pl10 + geom_vline(aes(xintercept=mean(log_output)+upper_ratio*sd(log_output)),
                        color="orange", linetype="dashed", size=1)
pl10



ggarrange(pl1, pl2, pl3, pl4, pl5, pl6, pl7, pl8, pl9, pl10,
          labels = c("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"),
          ncol = 4, nrow = 3)

ggsave("/Users/lding/Desktop/log2_output_density.pdf", width = 30, height = 16, dpi=300, useDingbats = FALSE)
