library(tidyverse)
library("ggplot2")
library("ggdendro")
library(proxy)
library(stringr)
library("RColorBrewer")
library("gplots")
library("ggheatmap")
library(reshape2)
library(ggnewscale)
getwd()

# clade info ----
clade_info <- read.csv('E:/adcadmic/Prof_guo/all_ec2/clade_info_2.csv', 
                         header = FALSE,
                         # row.names = 1
)
colnames(clade_info) <- c("sample", "clade")

color_1 <- c()#for clade order
color_2 <- c()
color_3 <- c()

for (x in rownames(clade_info)) {
  if (clade_info[x, "clade"] == "B2") {
    color_1 <- c(color_1, 'B2')
  } else if (clade_info[x, "clade"] == "N9") {
    color_1 <- c(color_1, 'N9')
  } else {
    color_1 <- c(color_1, 'Non-B2')
  }
}

for (x in rownames(clade_info)) {
  if (clade_info[x, "clade"] == "B2") {
    color_2 <- c(color_2, 1)
    
  } else if (clade_info[x, "clade"] == "D") {
    color_2 <- c(color_2, 2)
  } else if (clade_info[x, "clade"] == "E") {
    color_2 <- c(color_2, 3)
  } else if (clade_info[x, "clade"] == "A") {
    color_2 <- c(color_2, 4)
  } else if (clade_info[x, "clade"] == "B1") {
    color_2 <- c(color_2, 5)
  } else {
    color_2 <- c(color_2, 6)
  }
  # else if (clade_info[x, "clade"] == "N9") {
  #   color_2 <- c(color_2, 3)
  # } else {
  #   color_2 <- c(color_2, 2)
  # }
}

for (x in rownames(clade_info)) {
  if (clade_info[x, "clade"] == "B2") {
    color_3 <- c(color_3, 'B2')
  } else if (clade_info[x, "clade"] == "N9") {
    color_3 <- c(color_3, 'N9')
  } else {
    color_3 <- c(color_3, 'Non-B2')
  }
}


clade_info$color_1 <- color_1
clade_info$color_2 <- color_2
clade_info$color_3 <- color_3

# clade order by clade
clade_order <- clade_info$sample[order(clade_info$color_2)]
# clade_order

# get tables ----
# first_round_VF <- read.table('E:/adcadmic/Prof_guo/all_ec/VFDB/cov_filtered_vf_637.csv',
#                              header = FALSE,
#                              sep = ','
# )
# 
# colnames(first_round_VF) <- first_round_VF[1, ]
# first_round_VF <- first_round_VF[-c(1), ]
# 
# rownames(first_round_VF) <- first_round_VF[, 1]
# first_round_VF <- first_round_VF[, -c(1)]
# 
# second_round_VF <- read.table('E:/adcadmic/Prof_guo/all_ec/VFDB/to_predict_vf.csv',
#                               header = FALSE,
#                               sep = ','
# )
# 
# colnames(second_round_VF) <- second_round_VF[1, ]
# second_round_VF <- second_round_VF[-c(1), ]
# 
# rownames(second_round_VF) <- second_round_VF[, 1]
# second_round_VF <- second_round_VF[, -c(1)]
# 
# 
# needed_vfs <- c(rownames(first_round_VF), rownames(second_round_VF)) 
# needed_vfs <- unique(needed_vfs)
# 
# combined_tbl <- matrix(nrow = length(rownames(clade_info_2)),
#                        ncol = length(needed_vfs)) %>% as.data.frame()
# 
# rownames(combined_tbl) <- rownames(clade_info_2)
# colnames(combined_tbl) <- needed_vfs
# 
# combined_tbl[is.na(combined_tbl)] = 0.0 # 0 filling
# 
# for (x in rownames(combined_tbl)) {
#   if (grepl('GCA', x)) {
#     for (y in rownames(second_round_VF)) {
#       combined_tbl[x, y] <- second_round_VF[y, x]
#     }
#   } else {
#     for (y in rownames(first_round_VF)) {
#       combined_tbl[x, y] <- first_round_VF[y, x]
#     }
#   }
# }
# 
# to_keep_ <- c()
# 
# for (x in colnames(combined_tbl)) {
#   if (combined_tbl[, x] %>% as.numeric() %>% sum() > 0) {
#     to_keep_ <- c(to_keep_, x)
#   }
# }
# 
# combined_tbl <- combined_tbl[, to_keep_]
# 
# write.csv(combined_tbl,
#           'E:/adcadmic/Prof_guo/all_ec2/vf/1127-vf.csv',
#           row.names = TRUE)

# read the get table ----
heatmap_df_VF <- read.csv("E:/adcadmic/Prof_guo/all_ec2/vf/1127-vf.csv",
                          #sep = ",",
                          encoding = "UTF-8",
                          header = FALSE,
                          #row.names = FALSE
                          )

heatmap_df_VF[1, "V1"] <- "strain"
colnames(heatmap_df_VF) <- heatmap_df_VF[1, ]
heatmap_df_VF <- heatmap_df_VF[-c(1),]

# VFinfo ----
customized_vfinfo <- read.csv("E:/adcadmic/Prof_guo/5strains-20210824/92_EC/VFDB/bowtie2ed/processed_vfinfo.txt",
                              encoding = "UTF-8",
                              header = 1,
                              row.names = 1
)


vfinfo_needed <- data.frame()

# colnames(heatmap_df_VF)[-1][1]

for (x in row.names(customized_vfinfo)) {
  if (x %in% colnames(heatmap_df_VF)) {
    vfinfo_needed <- rbind(vfinfo_needed, customized_vfinfo[x,])
  }
}

vf_order <- rownames(vfinfo_needed)[order(vfinfo_needed$cate)]

# df melting ----

melted_df <- melt(heatmap_df_VF, id = "strain")

melted_df <- transform(melted_df, value = as.numeric(value))

melted_df$variable <- as.character(melted_df$variable)

melted_df_ <- melted_df

# plotting ----

VF_cate_color_hasp_map <- data.frame(vf_cate = vfinfo_needed$cate %>% table() %>% names() %>% as.vector(),
                                     color = c('#FF7474', '#FFFB25','#25FDFF', 
                                               '#BA25FF', '#FF2535', '#FF8825',
                                               '#B0F100', '#06DC66', '#C3C3C3', 
                                               '#2C3E50')
)
rownames(VF_cate_color_hasp_map) <- VF_cate_color_hasp_map$vf_cate

VF_cate_color <- c()

for (x in rownames(vfinfo_needed)) {
  VF_cate_color <- c(VF_cate_color,
                     VF_cate_color_hasp_map[
                       vfinfo_needed[x,"cate"],
                       "color"
                     ]
  )
}

hm <- ggplot(melted_df_) +
  geom_tile(aes(x = strain,
                y = variable,
                fill = value)) +
  scale_x_discrete(limits = clade_order) +
  scale_y_discrete(limits = vf_order) +
  #geom_raster(aes(fill = value), interpolate=TRUE) +
  #scale_fill_manual(breaks=c("[0,75)", "[75,80)", "[80,85)",
  #                           "[85,90)", "[90,95)", "[95,98)",
  #                           "[98,100)", "[100,Inf]"),
  #                    #c("\[0, 75)", "\[75, 80)", "\[80,85)", 
  #                    #       "\[85, 90)", "\[90, 95)", "\[95, 98)", 
  #                    #       "\[98, 100)"),
  #                  values = c("white", "lightblue", "blue",
  #                             "darkblue", "lightgreen", "green",
  #                             "darkgreen", "red")) +
  scale_fill_gradient2(low="blue", 
                       mid="gold", 
                       high="red", 
                       midpoint=87.5, #(100-75)/2 + 75
                       limits = c(75, 100),
                       na.value="white",#"transparent",
                       name = "Value"
  ) +
  theme_bw() +
  theme(axis.text = element_text(size = 4),
        # axis.text.x = element_text(angle = 90),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),#remove y axis labels
        axis.text.x = element_blank()
  ) +
  labs(x = "strain (n = 1127)", y = "VFs (n = 465)") +
  new_scale_fill() +
  geom_tile(data = clade_info,# stain clade color
            aes(x = sample, 
                y = -2,
                height = 5,
                fill = clade
            ),
            # fill = color_2
  ) +
  #scale_colour_gradient2(low="red",
  #                       mid="cyan",
  #                       high="purple",
  #                       midpoint=2,
  #                       limits = c(1, 3),
  #                       na.value="white"#"transparent"
  #                     ) +
  scale_fill_manual(name = "Clade",
                    values = c("purple", "salmon", "red", 
                               'green', 'blue', 'cyan')
  ) +
  #theme(axis.ticks.x = element_blank()) +
  new_scale_fill() +
  geom_tile(data = vfinfo_needed,# stain clade color
            aes(x = -20, 
                y = rownames(vfinfo_needed),
                width = 40,
                fill = cate
            ),
            # fill = VF_cate_color
  ) +
  scale_fill_manual(name = "Virulence Factor Category",
                    values = c('#FF7474', '#FFFB25','#25FDFF', 
                               '#BA25FF', '#FF2535', '#FF8825',
                               '#B0F100', '#06DC66', '#C3C3C3',
                               '#2C3E50')
  )

# vfinfo_needed$cate %>% table() %>% names()

#hm
ggsave("hm.pdf", 
       hm,
       height = 16, 
       width = 9)

