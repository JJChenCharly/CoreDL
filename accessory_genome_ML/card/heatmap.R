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
library(RColorBrewer)

# CARDDB df ----
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

# CARD order ----
CARD_info <- heatmap_df_CARD[, 1: 2]
rownames(CARD_info) <- CARD_info$gene

CARD_order <- CARD_info$gene[order(CARD_info$resistance_mechanism)]

# df melting ----

a <- as.vector(colnames(heatmap_df_CARD)[2: length(colnames(heatmap_df_CARD)
)
]
)

heatmap_df_CARD_gg <- t(rbind(a,
                              heatmap_df_CARD[, 2: length(colnames(heatmap_df_CARD
                              )
                              )
                              ]
)
) %>% data.frame()

heatmap_df_CARD_gg[1, 1] <- "strain"
colnames(heatmap_df_CARD_gg) <- heatmap_df_CARD_gg[1, ] 
heatmap_df_CARD_gg <- heatmap_df_CARD_gg[-c(1),]

melted_df <- melt(heatmap_df_CARD_gg, id = "strain")
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
  labs(x = "strain (n = 1127)", y = "AROs (n = 268)") +
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
ggsave("hm.pdf", 
       hm,
       height = 16, 
       width = 9)
