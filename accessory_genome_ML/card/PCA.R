library(devtools)
library(ggfortify)
library("dplyr")
library(vegan)
library(combinat)
library(ape)
library(ggpca)

# Read csv ----
heatmap_df_CARD <- read.csv("res_mech_filtered.csv",
                            #sep = ",",
                            encoding = "UTF-8",
                            header = FALSE,
                            #row.names = FALSE
                            )

# once
to_keep <- c()
to_keep_row <- c()
for (x in rownames(heatmap_df_CARD)) {
  if (!(heatmap_df_CARD[x, "V2"] %in% to_keep)) {
    to_keep_row <- c(to_keep_row, x)
    to_keep <- c(to_keep, heatmap_df_CARD[x, "V2"])
  }
}

heatmap_df_CARD <- heatmap_df_CARD[to_keep_row, ]

heatmap_df_CARD[1, "V2"] <- "strain"
colnames(heatmap_df_CARD) <- heatmap_df_CARD[1, ]
heatmap_df_CARD <- heatmap_df_CARD[-c(1),-c(1)]

rownames(heatmap_df_CARD) <- heatmap_df_CARD[, 1]
heatmap_df_CARD <- heatmap_df_CARD[, 2: length(colnames(heatmap_df_CARD))]

# clade info ----
clade_info <- read.csv("E:/adcadmic/Prof_guo/all_ec2/clade_info_2.csv",
                       header = FALSE)
colnames(clade_info) <- c("sample", "clade")

color_1 <- c()
color_2 <- c()

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
  } else {
    color_2 <- c(color_2, "purple")
  }
}



clade_info$color_1 <- color_1
clade_info$color_2 <- color_2

row.names(clade_info) <- clade_info[, 1]
clade_info <- clade_info[, 2: 4]
str(clade_info)


# transform to numeric
heatmap_df_CARD_to_num <- mutate_all(heatmap_df_CARD,
                                     function(x) as.numeric(x)
) %>% t() %>% data.frame()

colnames(heatmap_df_CARD_to_num) <- rownames(heatmap_df_CARD)


# clade info for heatmap ----
color_3 <- c()

for (x in rownames(heatmap_df_CARD_to_num)) {
  if (clade_info[x, "clade"] == "B2") {
    color_3 <- c(color_3, "red")
  } else if (clade_info[x, "clade"] == "N9") {
    color_3 <- c(color_3, "cyan")
  } else {
    color_3 <- c(color_3, "purple")
  }
}

heatmap_df_CARD_to_num$color_3 <- color_3
heatmap_df_CARD_to_num$color_2 <- color_2

heatmap_df_CARD_to_num$color_1 <- rep('1', 1127)

for (x in rownames(heatmap_df_CARD_to_num)) {
  heatmap_df_CARD_to_num[x, "color_1"] <- clade_info[x, "clade"]
}

# dist and permanova ----
heatmap_df_CARD_to_num[261]
heatmap_df_CARD_to_num[260]

bray_dist <- vegdist( heatmap_df_CARD_to_num[1: 260],
                      method = "bray"#"bray"
)
# as.matrix(bray_dist)

for_ggplot <- cmdscale(bray_dist, 
                       eig = TRUE, 
                       x.ret = TRUE)
for_ggplot_points <- for_ggplot$points
for_ggplot_df <- data.frame(Strain = rownames(for_ggplot_points),
                            X = for_ggplot_points[, 1],
                            Y = for_ggplot_points[, 2])

for_ggplot_df$group <- rep('1', 1127)
for (x in rownames(for_ggplot_df)) {
  clade_ <- clade_info[x, "clade"]
  # ifelse(clade_ == "B2",
  #        clade_ <- 'B2',
  #        clade_ <- 'Non-B2')
  for_ggplot_df[x, "group"] <- clade_
}

PcoA_gg <- ggplot(data = for_ggplot_df,
                  mapping = aes(x = X,
                                y = Y,
                                color = group
                  )
) + 
  # geom_mark_ellipse(aes(fill = group,
  #                       # label = group
  #                       ),
  #                   expand = unit(0.5,"mm"),
  #                   label.buffer = unit(-5, 'mm'),
  #                   ) +
  scale_fill_manual(name = 'Clade',
                    breaks = c("B2", "D", "E", "A", "B1", "N9"),
                    values = c("red", "green", "blue", 
                               'purple', 'salmon', 'cyan')
  ) +
  geom_point(size = 2) +
  scale_color_manual(name = 'Clade',
                     breaks = c("B2", "D", "E", "A", "B1", "N9"),
                     values = c("red", "green", "blue", 
                                'purple', 'salmon', 'cyan')
  ) +
  stat_ellipse(geom = "polygon",
               mapping = aes(fill = group),
               alpha = 0.2,
               level = 0.99
  ) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),) +
  labs(x = "PC1", y = "PC2")


ggsave("PCoA_bray_ready.pdf", 
       PcoA_gg,
       height = 9, 
       width = 16)
# sigtest ----
go <- combn(unique(heatmap_df_CARD_to_num$color_1),m = 2)


for (x in seq(1, 15)) {
  pair <- go[, x] %>% as.vector()
  
  keep_ <- c()
  
  for (y in rownames(heatmap_df_CARD_to_num)) {
    if (heatmap_df_CARD_to_num[y, 'color_1'] %in% pair) {
      keep_ <- c(keep_, y)
    }
  }
  tmp <- heatmap_df_CARD_to_num[keep_, ]
  
  sigtest <- adonis2(vegdist( tmp[1: 260],
                              method = "bray"#"bray"
  ) ~ tmp$color_1,  
  # data = env_for_perm,
  permutations = 999,
  # method="bray",
  by = "terms"
  )
  print(pair)
  print(sigtest$`Pr(>F)`)
}

# for (x in colnames(heatmap_df_CARD)) {
#   if (heatmap_df_CARD[, x] %>% as.numeric() %>% sum() == 0) {
#     print(x)
#   }
# }# no zero column
