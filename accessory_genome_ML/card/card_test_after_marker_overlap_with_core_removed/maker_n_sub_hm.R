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
clade_info <- read.csv("E:/adcadmic/Prof_guo/all_ec2/clade_info_2.csv",
                       header = FALSE)
colnames(clade_info) <- c("sample", "clade")

color_1 <- c()
color_2 <- c()
color_3 <- c()

for (x in rownames(clade_info)) {
  if (clade_info[x, "clade"] == "B2") {
    color_1 <- c(color_1, 1)
  } else if (clade_info[x, "clade"] == "N9") {
    color_1 <- c(color_1, 3)
  } else {
    color_1 <- c(color_1, 2)
  }
}

for (x in rownames(clade_info)) {
  if (clade_info[x, "clade"] == "B2") {
    color_2 <- c(color_2, "red")
  } else if (clade_info[x, "clade"] == "N9") {
    color_2 <- c(color_2, "cyan")
  } else {
    color_2 <- c(color_2, "purple")
  }
}

for (x in rownames(clade_info)) {
  if (clade_info[x, "clade"] == "B2") {
    color_3 <- c(color_3, 1)
    
  } else if (clade_info[x, "clade"] == "D") {
    color_3 <- c(color_3, 2)
  } else if (clade_info[x, "clade"] == "E") {
    color_3 <- c(color_3, 3)
  } else if (clade_info[x, "clade"] == "A") {
    color_3 <- c(color_3, 4)
  } else if (clade_info[x, "clade"] == "B1") {
    color_3 <- c(color_3, 5)
  } else {
    color_3 <- c(color_3, 6)
  }
  # else if (clade_info[x, "clade"] == "N9") {
  #   color_3 <- c(color_3, 3)
  # } else {
  #   color_3 <- c(color_3, 2)
  # }
}


clade_info$color_1 <- color_1
clade_info$color_2 <- color_2
clade_info$color_3 <- color_3

# clade order by clade
clade_order <- clade_info$sample[order(clade_info$color_3)]

# read the get table ----
heatmap_df_CARD <- read.csv("res_mech_filtered.csv",
                            #sep = ",",
                            encoding = "UTF-8",
                            header = FALSE,
                            #row.names = FALSE
)

colnames(heatmap_df_CARD) <- heatmap_df_CARD[1, ]
heatmap_df_CARD <- heatmap_df_CARD[-c(1),]

gene_multi_cat <- table(heatmap_df_CARD$gene)

# list those belongs to more than 1 cate
multi_cate <- c()

for(x in names(gene_multi_cat)) {
  if (gene_multi_cat[x] > 1) {
    multi_cate <- c(multi_cate, x)
    # print(x)
    # print(gene_multi_cat[x])
  }
}

# give unique marker

gene_unique <- heatmap_df_CARD$gene

# on all
for (to_process in multi_cate) {
  to_p <-  which(gene_unique %in% to_process)
  gene_unique[to_p[1]] <- substr(gene_unique[to_p[1]], 
                                 1,
                                 nchar(gene_unique[to_p[1]]) - 1
  )
  step <- 1
  for (posit in to_p[2: length(to_p)]) {
    #print(posit)
    #print(step)
    #print(to_p[posit])
    #print(paste0(gene_unique[to_p[1]],
    #             strrep("*", step)))
    gene_unique[posit] <- paste0(gene_unique[to_p[1]],
                                 strrep("*", step))
    step <- step + 1
  }
}

# verification of replication removal
length(gene_unique)#268
length(levels(factor(gene_unique)))#268

heatmap_df_CARD$gene <- gene_unique

# CARD info ----
CARD_info <- heatmap_df_CARD[, 1: 2]
rownames(CARD_info) <- CARD_info$gene# has double star

# heat map redo ----
heatmap_df_CARD$gene # has double star
heatmap_df_CARD_redo <- heatmap_df_CARD[, 2: 1129]
rownames(heatmap_df_CARD_redo) <- heatmap_df_CARD_redo[, 1]

heatmap_df_CARD_redo <- heatmap_df_CARD_redo %>% t() %>% as.data.frame()
heatmap_df_CARD_redo <- heatmap_df_CARD_redo[-c(1), ]


# marker table ----
marker_tab <- read.csv("E:/adcadmic/Prof_guo/all_ec2/card/card_test_after_marker_overlap_with_core_removed/res_marker.csv",
                       #sep = ",",
                       encoding = "UTF-8",
                       header = FALSE,
                       # row.names = FALSE
)

marker_tab[1, 1] <- 'strain'

colnames(marker_tab) <- marker_tab[1, ]
marker_tab <- marker_tab[-c(1), ]

marker_tab %>% colnames() # NO double star
# colnames(marker_tab) %>% table() %>% unname() %>% unlist()

# get needed columns ----
needed_cols <- c()
for (x in colnames(marker_tab)[2: 256]) {
  if (sum(marker_tab[, x] %>% as.numeric() != 0)) {
    needed_cols <- c(needed_cols, x)
  }
}

heatmap_df_CARD_redo <- heatmap_df_CARD_redo[, needed_cols]

needed_cols <- c('strain', needed_cols)
marker_tab <- marker_tab[, needed_cols]


# 'marA*' %in% needed_cols

# 'marA' %in% colnames(heatmap_df_CARD) # T
# 'marA' %in% needed_cols # F
# 'marA**' %in% colnames(heatmap_df_CARD) # F
# "Escherichia coli soxS with mutation conferring antibiotic resistance**"  %in% CARD_info$gene # T

# grepl('\\*', 'marA*')
# grepl('\\**', 'marA**')
# gsub('\\*', '', 'marA*')
# gsub('\\*', '\\*\\*', 'marA*')

for (x in needed_cols) {
  if (grepl('\\*', x)) {
    heatmap_df_CARD_redo[gsub('\\*', '', x)] <- heatmap_df_CARD_redo[x]
    marker_tab[gsub('\\*', '', x)] <- marker_tab[x]
  }
  
  if (grepl('\\*', x) & (gsub('\\*', '\\*\\*', x) %in% CARD_info$gene)) {
    heatmap_df_CARD_redo[gsub('\\*', '\\*\\*', x)] <- heatmap_df_CARD_redo[x]
    marker_tab[gsub('\\*', '\\*\\*', x)] <- marker_tab[x]
  }
}

# CARD order ----
CARD_info <- CARD_info[colnames(heatmap_df_CARD_redo), ]

CARD_order <- CARD_info$gene[order(CARD_info$resistance_mechanism)]

# df melting ----

heatmap_df_CARD_redo$strain <- rownames(heatmap_df_CARD_redo)

melted_df <- melt(heatmap_df_CARD_redo, id = "strain")

melted_df$variable <- as.character(melted_df$variable)

melted_df_ <- melted_df

# ggplot ----
hm <- ggplot(melted_df_) +
  geom_tile(aes(x = strain,
                y = variable,
                fill = value)) +
  scale_fill_manual(name = "RGI Score",
                    breaks = c('0', '1', '2') %>% rev(),
                    values = c('red', 'gold', 'lightblue')) +
  scale_x_discrete(limits = clade_order) +
  scale_y_discrete(limits = CARD_order) +
  # scale_fill_gradient2(low="lightblue",
  #                      mid="gold",
  #                      high="red",
  #                      midpoint=1, #(100-75)/2 + 75
  #                      limits = c(0, 2),
  #                      na.value="white"#"transparent"
  #                      ) +
  theme_bw() +
  theme(axis.text = element_text(size = 4),
        # axis.text.x = element_text(angle = 90,
        #                            size = 3),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
  ) +
  labs(x = "strain (n = 1127)", y = "AROs (n = 184)") +
  new_scale_fill() +
  geom_tile(data = clade_info,# stain clade color
            aes(x = sample,
                y = -1,
                height = 3,
                fill = clade
            )
  ) +
  scale_fill_manual(name = "Clade",
                    values = c("purple", "salmon", "red", 
                               'green', 'blue', 'cyan')
  ) +
  new_scale_fill() +
  geom_tile(data = CARD_info,
            aes(x = -25,
                y = gene,
                width = 50,
                fill = resistance_mechanism
            )
  ) +
  scale_fill_manual(name = "Resistance Mechanism",
                    values = brewer.pal(n = 6,
                                        name = 'Accent')
                    # c('#FF7474', '#FFFB25','#25FDFF',
                    #          '#BA25FF', '#FF2535', '#FF8825')
  )

#hm
ggsave("E:/adcadmic/Prof_guo/all_ec2/card/card_test_after_marker_overlap_with_core_removed/sub_heatmap.pdf", 
       hm,
       height = 16, 
       width = 9)


# marker heatmap ----
# marker tab normalization or nudge ----
marker_tab_wei <- marker_tab[, 2: length(colnames(marker_tab))]
marker_tab_wei <- sapply( marker_tab_wei, as.numeric ) %>% as.data.frame()

# marker_tab_wei[2, "TEM-206"]
# marker_tab_wei[2, ] %>% range()

# marker_tab_wei[1, ] %>% unname() %>% unlist() %>% unique() %>% sort()



for (x in rownames(marker_tab_wei)) {
  range_ <-  marker_tab_wei[x, ] %>% unname() %>% unlist() %>% unique() %>% sort()
  
  if (range_ %>% min() == 0) {# if all weights > 0
    for (y in colnames(marker_tab_wei)) {
      marker_tab_wei[x, y] <- marker_tab_wei[x, y]/max(range_)
    }
  } else if (range_ %>% max() == 0) {# if all weights < 0
    for (y in colnames(marker_tab_wei)) {
      marker_tab_wei[x, y] <- (marker_tab_wei[x, y]/min(range_))*-1
    }
  } else {
    for (y in colnames(marker_tab_wei)) {
      if (marker_tab_wei[x, y] <= 0) {
        marker_tab_wei[x, y] <- (marker_tab_wei[x, y]/min(range_))*-1
      } else {
        marker_tab_wei[x, y] <- marker_tab_wei[x, y]/max(range_)
      }
      
    }
  }
  
}

marker_tab_wei[1, ] %>% unname() %>% unlist() %>% unique() %>% sort()

marker_tab_wei$strain <- marker_tab$strain

# melt and plot ----
marker_melt <- melt(marker_tab_wei, id = "strain")

str(marker_melt)

marker_melt$variable <- as.character(marker_melt$variable)

marker_melt$value <- as.numeric(marker_melt$value)

marker_melt$value %>% unique() %>% sort()

hm_marker <- ggplot(marker_melt) +
  geom_tile(aes(x = strain,
                y = variable,
                fill = value)) +
  scale_x_discrete(limits = clade_info$clade %>% unique()) +
  scale_y_discrete(limits = CARD_order) +
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
        # axis.text.y = element_blank(),#remove y axis labels
        axis.text.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
  ) +
  labs(x = "strain (n = 1127)", y = "AROs (n = 184)") +
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
  geom_tile(data = CARD_info,# stain clade color
            aes(x = 0.4, 
                y = gene,
                width = 0.2,
                fill = resistance_mechanism
            ),
            # fill = VF_cate_color
  ) +
  scale_fill_manual(name = "Resistance Mechanism",
                    values = brewer.pal(n = 6,
                                        name = 'Accent')
  )
ggsave("E:/adcadmic/Prof_guo/all_ec2/card/card_test_after_marker_overlap_with_core_removed/marker_heatmap.pdf", 
       hm_marker,
       height = 16, 
       width = 9)
