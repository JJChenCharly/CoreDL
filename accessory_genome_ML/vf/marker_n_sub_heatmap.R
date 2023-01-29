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
# getwd()

# clade info ----
clade_info <- read.csv("E:/adcadmic/Prof_guo/all_ec2/clade_info_2.csv",
                       header = FALSE)
colnames(clade_info) <- c("sample", "clade")

color_2 <- c()#for clade order
color_3 <- c()

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
  } else {
    color_3 <- c(color_3, 'Non-B2')
  }
}

clade_info$color_2 <- color_2
clade_info$color_3 <- color_3


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

for (x in row.names(customized_vfinfo)) {
  if (x %in% colnames(heatmap_df_VF)) {
    vfinfo_needed <- rbind(vfinfo_needed, customized_vfinfo[x,])
  }
}

vf_order <- rownames(vfinfo_needed)[order(vfinfo_needed$cate)]

# marker table ----
marker_tab <- read.csv("Vf_marker.csv",
                       #sep = ",",
                       encoding = "UTF-8",
                       header = FALSE,
                       # row.names = FALSE
)
marker_tab[1, 1] <- 'strain'

colnames(marker_tab) <- marker_tab[1, ]
marker_tab <- marker_tab[-c(1), ]

# get needed columns ----
needed_cols <- c('strain')

for (x in colnames(marker_tab)[2: 465]) {
  if (sum(marker_tab[, x] %>% as.numeric() != 0)) {
    needed_cols <- c(needed_cols, x)
  }
}

heatmap_df_VF <- heatmap_df_VF[, needed_cols]
marker_tab <- marker_tab[, needed_cols]

# clade order ----
clade_order <- clade_info$sample[order(clade_info$color_2)]

# VFinfo ----
customized_vfinfo <- read.csv("E:/adcadmic/Prof_guo/5strains-20210824/92_EC/VFDB/bowtie2ed/processed_vfinfo.txt",
                              encoding = "UTF-8",
                              header = 1,
                              row.names = 1
)


vfinfo_needed <- data.frame()

for (x in row.names(customized_vfinfo)) {
  if (x %in% colnames(heatmap_df_VF)[-1]) {
    vfinfo_needed <- rbind(vfinfo_needed, customized_vfinfo[x,])
  }
}

# vfinfo_needed$pro_Name %>% unique() %>% length()
vf_order <- rownames(vfinfo_needed)[order(vfinfo_needed$cate)]

# df melting ----

melted_df <- melt(heatmap_df_VF, id = "strain")

melted_df <- transform(melted_df, value = as.numeric(value))

melted_df$variable <- as.character(melted_df$variable)

melted_df_ <- melted_df

# ggplt2 ----
vfinfo_needed$cate %>% table() %>% names() %>% length()
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

melted_df_$value %>% unique() %>% sort()

hm <- ggplot(melted_df_) +
  geom_tile(aes(x = strain,
                y = variable,
                fill = value)) +
  scale_x_discrete(limits = clade_order) +
  scale_y_discrete(limits = vf_order) +
  scale_fill_gradient2(low="blue", 
                       mid="gold", 
                       high="red", 
                       midpoint=90, #(100-75)/2 + 75
                       limits = c(80, 100),
                       na.value="white",#"transparent",
                       name = "Value"
  ) +
  theme_bw() +
  theme(axis.text = element_text(size = 4),
        # axis.text.x = element_text(angle = 90),
        # axis.text.y = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),#remove y axis labels
        axis.text.x = element_blank()
  ) +
  labs(x = "strain (n = 1127)", y = "VFs (n = 253)") +
  new_scale_fill() +
  geom_tile(data = clade_info,# stain clade color
            aes(x = sample, 
                y = -1,
                height = 3,
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
            aes(x = -10, 
                # y = rownames(vfinfo_needed),
                y = rownames(vfinfo_needed),
                width = 20,
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
ggsave("E:/adcadmic/Prof_guo/all_ec2/vf/sub_heatmap.pdf", 
       hm,
       height = 9, 
       width = 16)

# marker heatmap ----
# marker tab normalization or nudge ----
marker_tab_wei <- marker_tab[, 2: length(colnames(marker_tab))]
marker_tab_wei <- sapply( marker_tab_wei, as.numeric ) %>% as.data.frame()

# marker_tab_wei[1, ] %>% unname() %>% unlist() %>% unique() %>% sort()
marker_tab_wei[1, ] %>% unname() %>% unlist() %>% unique() %>% abs() %>% max()

# for (x in rownames(marker_tab_wei)) {
#   range_ <-  marker_tab_wei[x, ] %>% unname() %>% unlist() %>% unique() %>% abs() %>% max()
#   
#   if (range_ %>% min() == 0) {# if all weights > 0
#     for (y in colnames(marker_tab_wei)) {
#       marker_tab_wei[x, y] <- marker_tab_wei[x, y]/max(range_)
#     }
#   } else if (range_ %>% max() == 0) {# if all weights < 0
#     for (y in colnames(marker_tab_wei)) {
#       marker_tab_wei[x, y] <- (marker_tab_wei[x, y]/min(range_))*-1
#     }
#   } else {
#     for (y in colnames(marker_tab_wei)) {
#       if (marker_tab_wei[x, y] <= 0) {
#         marker_tab_wei[x, y] <- (marker_tab_wei[x, y]/min(range_))*-1
#       } else {
#         marker_tab_wei[x, y] <- marker_tab_wei[x, y]/max(range_)
#       }
#       
#     }
#   }
#   
# }

for (x in rownames(marker_tab_wei)) {
  limit_ <-  marker_tab_wei[x, ] %>% unname() %>% unlist() %>% unique() %>% abs() %>% max()
  
  for (y in colnames(marker_tab_wei)) {
    marker_tab_wei[x, y] <- marker_tab_wei[x, y]/limit_
  }

}

marker_tab_wei[1, ] %>% unname() %>% unlist() %>% unique() %>% sort()

marker_tab_wei$strain <- marker_tab$strain
# melt and plot ----
marker_melt <- melt(marker_tab_wei, id = "strain")

marker_melt$variable <- as.character(marker_melt$variable)

marker_melt$value <- as.numeric(marker_melt$value)

marker_melt$value %>% unique() %>% sort()

hm_marker <- ggplot(marker_melt) +
  geom_tile(aes(x = strain,
                y = variable,
                fill = value)) +
  scale_x_discrete(limits = marker_tab$strain) +
  scale_y_discrete(limits = vf_order) +
  scale_fill_gradient2(low="blue", 
                       mid="grey", 
                       high="red", 
                       midpoint=0, #(100-75)/2 + 75
                       limits = c(-1.001, 1.001),
                       na.value="white",#"transparent",
                       name = "weights"
  ) +
  theme_bw() +
  theme(axis.text = element_text(size = 4),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),#remove y axis labels
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  ) +
  labs(x = "clade", y = "VFs (n = 253)") +
  new_scale_fill() +
  geom_tile(data = clade_info,# stain clade color
            aes(x = clade, 
                y = -1.1,
                height = 3,
                fill = clade
            ),
            # fill = color_2
  ) +
  scale_fill_manual(name = "Clade",
                    values = c("purple", "salmon", "red", 
                               'green', 'blue', 'cyan')
  ) +
  #theme(axis.ticks.x = element_blank()) +
  new_scale_fill() +
  geom_tile(data = vfinfo_needed,# stain clade color
            aes(x = 0.4, 
                y = rownames(vfinfo_needed),
                width = 0.2,
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
ggsave("E:/adcadmic/Prof_guo/all_ec2/vf/marker_heatmap.pdf", 
       hm_marker,
       height = 16, 
       width = 9)
