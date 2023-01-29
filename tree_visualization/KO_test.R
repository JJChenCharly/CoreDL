library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(phytools)
library(r2r)
library(gmodels)

# getwd()

# functions ----
# get every under set child first ----
get_nodes_under <- function(node_No, the_ggplot) {
  temp_df <- the_ggplot$data
  
  temp_df$c1 <- rep('a', length(rownames(temp_df)
  )
  )
  temp_df$c2 <- rep('a', length(rownames(temp_df)
  )
  )
  
  for (x in rownames(temp_df)) {
    parent_ <- temp_df[x, "parent"] %>% as.integer()
    
    if (temp_df[parent_, 'c1'] == 'a') {
      temp_df[parent_, 'c1'] <- x
    } else if (temp_df[parent_, 'c2'] == 'a') {
      temp_df[parent_, 'c2'] <- x
    } else {
      print(x)
    }
  }
  
  node_for_next_loop <- c(node_No)
  node_coll <- c(node_No)
  
  while (length(node_for_next_loop) > 0) {
    node_for_next_loop_togo <- node_for_next_loop
    node_for_next_loop <- c()
    
    for (nodes in node_for_next_loop_togo) {
      
      
      if (temp_df[nodes, 'c1'] == 'a') {
        next
      } else {
        node_for_next_loop <- c(node_for_next_loop, as.integer(temp_df[nodes, 'c1']
        )
        )
        node_coll <- c(node_coll,  as.integer(temp_df[nodes, 'c1']
        )               
        )
      }
      
      if (temp_df[nodes, 'c2'] == 'a') {
        next
      } else {
        node_for_next_loop <- c(node_for_next_loop, as.integer(temp_df[nodes, 'c2']
        )
        )
        node_coll <- c(node_coll,  as.integer(temp_df[nodes, 'c2']
        )
        )
      }
      
    }
  }
  
  return(node_coll)
}

# get tip under ----
get_tip_under <- function(char_lst_of_members, total_tip_count) {
  coll <- c()
  
  for (x in char_lst_of_members) {
    if (as.numeric(x) <= total_tip_count) {
      coll <- c(coll, x)
    }
  }
  
  return(coll)
}

# bootstrap value reading ----
# bootstrap_ <- read.raxml(file = "rooted-op/RAxML_bipartitionsBranchLabels.1135")
bootstrap_ <- read.raxml(file = "KOed_top1/raxml/RAxML_bipartitionsBranchLabels.KOed_1127")

# a quick look ----
tree_1127 <- bootstrap_

# clade color 
tip_color <- tree_1127@phylo$tip.label

for (x in seq(1, length(tip_color))) {
  if (grepl("YN-", tip_color[x])) {
    tip_color[x] <- "red"
  } else if (grepl("GCA", tip_color[x])) {
    tip_color[x] <- "black"
  } else {
    tip_color[x] <- "#0ABAB5"
  }
}

tree_plt <- ggtree(tree_1127,
                   branch.length='none'
                   # layout = "circular",
                   # open.angle = 90
) +
  #geom_label(aes(label = node),
  #           size = 2) +
  geom_text(#aes(label=ifelse(bootstrap < 95, bootstrap, NA)),
    # aes(label=ifelse(bootstrap < 95, bootstrap, NA)),
    aes(label=node),
    size = 1.5,
    hjust = 0
  ) +
  geom_tiplab(size = 2,
              # hjust = -5,
              align = FALSE,# False
              color = tip_color
  )  
#geom_rootedge(rootedge = 0.005)
#geom_label(aes(label = label),
#           size = 2,
#           hjust = -0.5)

ggsave("KOed_top1/temp_tree.pdf",
       tree_plt,
       height = 80,
       width = 40,
       limitsize = FALSE
)

# reroot ----
tree_reroot <- root(tree_1127, 
                    #outgroup = c(247),
                    node = 1332,
                    #resolve.root = TRUE,#!!!!!!!!!!
                    #interactive = TRUE
                    #edgelabel = TRUE
)

tip_color <- tree_1127@phylo$tip.label

for (x in seq(1, length(tip_color))) {
  if (grepl("YN-", tip_color[x])) {
    tip_color[x] <- "red"
  } else if (grepl("GCA", tip_color[x])) {
    tip_color[x] <- "black"
  } else {
    tip_color[x] <- "#0ABAB5"
  }
}

# #1F18C0 klein blue
# #0ABAB5 tiffanny blue
# #86FFFF lighter tiffanny blue

tree_plt_reroot <- ggtree(tree_reroot,
                          # colour = rep("cyan", 272)
                          # layout = "circular",
                          # open.angle = 90,
                          # layout="ellipse", 
                          branch.length="none"
) +
  #geom_label(aes(label = node),
  #           size = 2) +
  geom_text(
    # aes(label=ifelse(bootstrap < 95, bootstrap, NA)),
    aes(label = node),
    size = 2,
    hjust = -0.2,
    #fontface = "bold"
    #vjust = -1
  ) +
  geom_tiplab(size = 3,
              #hjust = -3,
              align = TRUE,
              color = tip_color
  ) #+
# geom_rootedge(rootedge = 0.002)


ggsave("KOed_top1/temp_tree.pdf",
       tree_plt_reroot,
       height = 100,
       width = 20,
       limitsize = FALSE
)

# make new node ----
tree_plt_reroot$data

tree_plt_reroot$data[2252, ]# last row
tree_plt_reroot$data[2253, ]# empty row

tree_plt_reroot$data[2253, "parent"] <- 1128# reduced tree has only 1127 tipes
tree_plt_reroot$data[2253, "node"] <- 2253
tree_plt_reroot$data[2253, "branch.length"] <- 0.002825
tree_plt_reroot$data[2253, "isTip"] <- FALSE
tree_plt_reroot$data[2253, "x"] <- 0 # temporary
tree_plt_reroot$data[2253, "y"] <- tree_plt_reroot$data[c(1668, 1667), "y"] %>% sum()/2
tree_plt_reroot$data[2253, "branch"] <- 0
tree_plt_reroot$data[2253, "angle"] <- 0

tree_plt_reroot$data[1668, "parent"] <- 2253
tree_plt_reroot$data[1667, "parent"] <- 2253

ggsave("KOed_top1/temp_tree.pdf",
       tree_plt_reroot,
       height = 100,
       width = 15,
       limitsize = FALSE
)

nodes_under_2253 <- get_nodes_under(node_No = 2253, 
                                    the_ggplot = tree_plt_reroot)

nodes_under_1478 <- get_nodes_under(node_No = 1478, 
                                    the_ggplot = tree_plt_reroot)# NON-B2 node

# clade position ----
# tree_plt_reroot$data[1128, ] # x = 0, y = 262.
tree_plt_reroot$data[1478, ]
correction <- tree_plt_reroot$data[1478, "x"] %>% as.double()

for (x in as.character(nodes_under_2253)) {
  tree_plt_reroot$data[x, "x"] <- tree_plt_reroot$data[x, "x"] + (correction/2)
}

for (x in as.character(nodes_under_1478)) {
  tree_plt_reroot$data[x, "x"] <- tree_plt_reroot$data[x, "x"] - (correction/2)
}

ggsave("KOed_top1/temp_tree.pdf",
       tree_plt_reroot 
       %>% flip(1478, 2253)
       %>% flip(1479, 1477)
       %>% flip(1475, 2041)
       # %>% flip(1296, 2052)
       # %>% flip(1307, 1305)
       ,
       height = 80,
       width = 40,
       limitsize = FALSE
)

# tree_plt_reroot$data[304, ]$label
# tree_plt_reroot$data[300, ]$label
# tree_plt_reroot$data[117, ]$label


# clade info save ----
flipped_tree <- tree_plt_reroot %>% flip(1478, 2253) %>% flip(1479, 1477) %>% flip(1475, 2041)

df <- data.frame(node_no = c(2253, 1479, 2041, 1254, 1475, 1856))
rownames(df) <- c("B2", "D", "E", "A", "B1", "N9")

clade_info <- data.frame(strain = c('UTI89'),
                         clade = c('B2'))

# 'UTI89' %in% clade_info[, 'strain']

for (x in rownames(df)) {
  # print(typeof(x))
  ndoes <- get_nodes_under(df[x, "node_no"], 
                           flipped_tree)
  
  tips <- get_tip_under(ndoes, 1127) %>% as.vector()
  
  for (y in tips) {
    if (flipped_tree$data$label[y] %in% clade_info[, 'strain']) {
      next
    } else {
      clade_info <- rbind(clade_info, data.frame(strain = c(flipped_tree$data$label[y]),
                                                 clade = c(x)
      )
      )
    }
  }
}

write.table(clade_info,
            "KOed_top1/clade_info_koed.csv",
            sep = ",",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# comparison and comfirmation ----
clade_info_2 <- read.csv('E:/adcadmic/Prof_guo/all_ec2/clade_info_2.csv', 
                         header = FALSE,
                         row.names = 1
)

clade_info_KO <- read.csv('E:/adcadmic/Prof_guo/all_ec2/KOed_top1/clade_info_koed.csv', 
                          header = FALSE,
                          row.names = 1
)

for (x in rownames(clade_info_2)) {
  if (clade_info_2[x, "V2"] != clade_info_KO[x, "V2"]) {
    print(x)
  }
}