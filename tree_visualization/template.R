library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(phytools)
library(r2r)
# library(gmodels)

# getwd()

# functions ----
# search its nodes ----
# get every under ----
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
      # print(x
      next
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

# get_nodes_under <- function(node_No, the_ggplot, tip_count) {
#   tip_node_under <- c(node_No)
#   
#   in_togo <- c(node_No)
#   count_ <- 0
#   in_togo2 <- c()
#   
#   while (length(in_togo) > 0) {
#     
#     if (count_ == 0) {
#       in_togo <- c(node_No)
#     } else {
#       in_togo <- c()
#       in_togo <- in_togo2
#       in_togo2 <- c()
#     }
#     
#     for (x in in_togo) {
#       for (p in row.names(the_ggplot$data)) {
#         if (the_ggplot$data[p, "parent"] == x) {
#           tip_node_under <- c(tip_node_under, the_ggplot$data[p, "node"])
#           
#           if (the_ggplot$data[p, "node"] > tip_count) {
#             in_togo2 <- c(in_togo2, the_ggplot$data[p, "node"]) 
#           }
#         }
#       }
#     }
#     count_ <- count_ + 1
#   }
#   
#   return(tip_node_under)
# }

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


# return extend ----
get_extend <- function(the_plot,node_id){
  temp_df <- the_plot$data
  
  tip_count <- temp_df$isTip %>% sum()
  
  gloubal_max <- temp_df[1: tip_count, 'x'] %>% max()
  
  tip_list <- get_nodes_under(the_ggplot = the_plot,
                              node_No = node_id) %>% get_tip_under(tip_count)
  
  local_max <- temp_df[tip_list, 'x'] %>% max()
  
  return (gloubal_max - local_max)
}

# tree_plt_reroot$data$isTip %>% sum()
# test_m <- tree_plt_reroot$data[1: 544, 'x'] %>% max()
# 
# test <- get_nodes_under(651, tree_plt_reroot) %>% get_tip_under(544)
# 
# test_v <- tree_plt_reroot$data[test, 'x'] %>% max()
# 
# test_m - test_v

# bootstrap value reading ----
# bootstrap_ <- read.raxml(file = "rooted-op/RAxML_bipartitionsBranchLabels.1135")
bootstrap_ <- read.raxml(file = "raxml/RAxML_bipartitionsBranchLabels.core_genome_S_aureus_600")
# bootstrap_ @phylo %>% str()

# str(bootstrap_)
# bootstrap_@data$bootstrap
# bootstrap_@data$node

# a quick look ----
#tree_599 <- read.tree("91/RAxML_bestTree.91")# no bootstrap info
tree_599 <- bootstrap_


# class(tree_599)
# tree_599[1]
# tree_599

# clade color 

tip_color <- tree_599@phylo$tip.label

for (x in seq(1, length(tip_color))) {
  if (grepl("YN-", tip_color[x])) {
    tip_color[x] <- "red"
  } else if (grepl("GCA", tip_color[x])) {
    tip_color[x] <- "black"
  } else {
    tip_color[x] <- "#0ABAB5"
  }
}

tree_plt <- ggtree(tree_599,
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

# tree_plt

ggsave("temp_tree.pdf",
       tree_plt,
       height = 80,
       width = 35,
       limitsize = FALSE
)

# node 1689 parent of B2

# drop some tips not sure which clade to go, draw the quick tree first ----
to_drop <- c()

to_drop <- c(to_drop, tree_plt$data[1013, "label"]$label)

nodes_under_1087 <- get_nodes_under(node_No = 1013, 
                                    the_ggplot = tree_plt
                                    # tip_count = 1135
)
tips_under_1087 <- get_tip_under(char_lst_of_members = nodes_under_1087, 600)#找TIP需要小于总菌数，这个数代表菌数


to_drop <- c(to_drop, tree_plt$data[tips_under_1087, "label"]$label)

to_drop %>% length()
to_drop
# reduce bootstrap
bootstrap_ <- drop.tip(bootstrap_, to_drop)

# second look ----
#tree_599 <- read.tree("91/RAxML_bestTree.91")# no bootstrap info
tree_590 <- bootstrap_

# class(tree_590)
# tree_590[1]
# tree_590

# clade color 
tip_color <- tree_590@phylo$tip.label

for (x in seq(1, length(tip_color))) {
  if (grepl("YN-", tip_color[x])) {
    tip_color[x] <- "red"
  } else if (grepl("GCA", tip_color[x])) {
    tip_color[x] <- "black"
  } else {
    tip_color[x] <- "#0ABAB5"
  }
}

tree_plt <- ggtree(tree_590,
                   # layout = "circular",
                   # open.angle = 90
) +
  #geom_label(aes(label = node),
  #           size = 2) +
  geom_text(
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

# tree_plt

ggsave("temp_tree.pdf",
       tree_plt,
       height = 80,
       width = 35,
       limitsize = FALSE
)


# reroot ----
tree_reroot <- root(tree_590, 
                    #outgroup = c(247),
                    node = 606,
                    #resolve.root = TRUE,#!!!!!!!!!!
                    #interactive = TRUE
                    #edgelabel = TRUE
)

tip_color <- tree_590@phylo$tip.label

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
                          # branch.length="none"
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


ggsave("temp_tree.pdf",
       tree_plt_reroot,
       height = 100,
       width = 20,
       limitsize = FALSE
)



# make new node ----
tree_plt_reroot$data

tree_plt_reroot$data[1194, ]# last row
tree_plt_reroot$data[1195, ]# empty row

# new node between 1488 1863 both non B2 of tree_plt_reroot

tree_plt_reroot$data[1195, "parent"] <- 599# reduced tree has only 1127 tipes
tree_plt_reroot$data[1195, "node"] <- 1195
tree_plt_reroot$data[1195, "branch.length"] <- 0.002825
tree_plt_reroot$data[1195, "isTip"] <- FALSE
tree_plt_reroot$data[1195, "x"] <- 0 # temporary
tree_plt_reroot$data[1195, "y"] <- tree_plt_reroot$data[c(606,607), "y"] %>% sum()/2
tree_plt_reroot$data[1195, "branch"] <- 0
tree_plt_reroot$data[1195, "angle"] <- 0

# change parent of 1488 1863
tree_plt_reroot$data[606, "parent"] <- 1195
tree_plt_reroot$data[607, "parent"] <- 1195

ggsave("temp_tree.pdf",
       tree_plt_reroot,
       height = 100,
       width = 15,
       limitsize = FALSE
)

nodes_under_1087 <- get_nodes_under(node_No = 1195, 
                                    the_ggplot = tree_plt_reroot
                                    # tip_count = 1127
)
# nodes_under_1087 %>% length()# 319
# "368" %in% nodes_under1087

nodes_under_1674 <- get_nodes_under(node_No = 611, 
                                    the_ggplot = tree_plt_reroot
                                    # 1127
)# B2 node
# nodes_under_1674[[1]][1] %>% typeof()

# tree_plt_reroot$data[c(1127, 1128), ]

tree_plt_reroot$data$isTip %>% sum()

a <- get_tip_under(nodes_under_1087, 598)
a %>% length()

b <- get_tip_under(nodes_under_1674, 598)
b %>% length()

# clade position ----
# tree_plt_reroot$data[807, ] # x = 0, y = 262.
# tree_plt_reroot$data[817, ]
correction <- tree_plt_reroot$data[611, "x"] %>% as.double()

for (x in as.character(nodes_under_1087)) {
  tree_plt_reroot$data[x, "x"] <- tree_plt_reroot$data[x, "x"] + (correction/2)
}

for (x in as.character(nodes_under_1674)) {
  tree_plt_reroot$data[x, "x"] <- tree_plt_reroot$data[x, "x"] - (correction/2)
}

ggsave("temp_tree.pdf",
       tree_plt_reroot 
       # %>% flip(1674, 1087)
       # %>% flip(1863, 1488)
       # %>% flip(1297, 1299)
       # %>% flip(1296, 2052)
       # %>% flip(1307, 1305)
       ,
       height = 80,
       width = 40,
       limitsize = FALSE
)

# tree_plt_reroot$data[1087, ]


# draw ----
aes_df <- data.frame(id=c(1195,612,739,810),
                     annote=c("A", "B", "C", "D"),
                     annote_no = c("", "", "", ""),
                     colr=c("18-1552", "14-3205", "13-5412", "16-4535"),
                     offset_1=c(0.01, 0.01, 0.01, 0.01),
                     offset_2=rep(0.0135, 4),#0.0135
                     offset_text=rep(-0.003, 4)
)

# set1 <- c(1, 2, 3, 4, 5)
# set2 <- c(2, 4)
# result <- setdiff(set1, set2)

extend_c <- c()
for (x in aes_df$id){
  extend_c <- c(extend_c, get_extend(tree_plt_reroot, x))
  
}
extend_c
# 0.02781408 > 0.02806981 

flipped_tree <- tree_plt_reroot + #%>% flip(1674, 1087) %>% flip(1863, 1488) %>% flip(1297, 1299) %>% flip(1296, 2052)
  geom_cladelab(
    #data = aes_df,
    #mapping = aes(
    #node = id,
    #label = annote_no,
    #colour = colr,
    #offset = offset,
    #offset.text = offset_text,
    #hjust = 10
    #),
    align = TRUE,
    node = aes_df$id,
    # label = c("B2", "D", "E", "A", "B1", "N9"),
    label = rep(0, 4),
    offset = rep(0.016, 4),#0.0113
    barcolor = c("18-1552", "14-3205", "13-5412", "#ffff00"),
    barsize = 80,
    extend = rep(0.5, 4),
    show.legend = FALSE,
  ) +
  geom_cladelab(
    #data = aes_df,
    #mapping = aes(
    #  node = id,
    #  label = annote,
    #  offset = offset_2,
    #  offset.text = offset_text,
    #  hjust = rep(0.5, 6),
    #  ),
    node = aes_df$id,
    label = c("A", "B", "C", "D"),
    offset = rep(0, 4),
    offset.text =rep(0.0145, 4),
    textcolour = "black",
    fontsize = 50, # size of clade bar txt
    align = TRUE,
    barcolor = "white",
    hjust = rep(-0.1, 4),
    barsize = 0,
    angle = 90,
    show.legend = FALSE
  ) +
  #geom_balance(node = 138,
  #             fill = c("red", "purple"),
  #             alpha = .1,
  #             color = "white"
  #extend = .001
  #             )
  geom_hilight(node= aes_df$id, 
               # fill=c("red", "purple", "cyan"), 
               fill= c("18-1552", "14-3205", "13-5412", "#ffff00"),
               alpha=.1, 
               # extend = c(0.0205, 0, 0, 0, 0, 0.0203),
               extend = extend_c,
               # extend = rep(0, 4),
               align = "left"
  ) + 
  geom_treescale(x = 0.01,
                 y = 600)
#geom_text(aes(label=node),
#          size = 1.5,
#          hjust = - 1) 

ggsave("annotated.pdf",
       flipped_tree,
       height = 80,
       width = 40,
       limitsize = FALSE
)

# set1 <- c(1, 2, 3, 4, 5)
# set2 <- c(2, 4)
# result <- setdiff(set1, set2)

tip_count <- tree_plt_reroot$data$isTip %>% sum()
set1 <- c(get_nodes_under(739,tree_plt_reroot)) %>% get_tip_under(tip_count)
set2 <- c(get_nodes_under(810,tree_plt_reroot)) %>% get_tip_under(tip_count)
class3 <- setdiff(set1, set2)
class3 %>% max()
# member count ----
for (x in aes_df$id) {
  nodes <- get_nodes_under(x, flipped_tree)
  tips <- get_tip_under(nodes, 598)
  
  print(paste0(as.character(x), " has >>>>>>"
  )
  )
  print(as.character(length(tips)))
}

# [1] "1674 has >>>>>>" B2
# [1] "190"
# [1] "1863 has >>>>>>" D
# [1] "190"
# [1] "1299 has >>>>>>" E
# [1] "190"
# [1] "2052 has >>>>>>" A
# [1] "183"
# [1] "1296 has >>>>>>" B1
# [1] "188"
# [1] "1489 has >>>>>>" N9
# [1] "186"
# 190*3 + 183 + 188 + 186 == 1127 # TRUE
# clade info save ----
df <- data.frame(node_no = aes_df$id)
rownames(df) <- c("A", "B", "C", "D")

clade_info <- list()
  
  # data.frame(strain = c(),
  #            clade = c())
for (x in rownames(df)) {
  # print(typeof(x))
  ndoes <- get_nodes_under(df[x, "node_no"], 
                           flipped_tree)
  tips <- get_tip_under(ndoes, 598) %>% as.vector()
  
  for (y in tips) {
    clade_info[flipped_tree$data$label[y]] <- x
  }
}

# names(clade_info) %>% length()
# clade_info["GCA_001019395.2_ASM101939v2_genomic.fna"]

to_save <- data.frame(value = clade_info %>% unlist() %>% as.vector())
row.names(to_save) <- names(clade_info)


write.table(to_save,
            "S_aureus_598.csv",
            sep = ",",
            quote = FALSE,
            row.names = TRUE,
            col.names = FALSE)









# compare with the prediction result from the first round and get confusion matirx ----

res_pre_based_on_1 <- read.csv('E:/adcadmic/Prof_guo/all_ec/VFDB/res_of_predict.csv',
                               # row.names = 'X'
)
needed_id <- flipped_tree$data[1: 1127, ]$label
needed_id2 <- c()

for (x in seq(1, length(needed_id))) {# remove Non GCA
  if (grepl('GCA', needed_id[x])) {
    needed_id2 <- c(needed_id2, needed_id[x])
  } else {
    next
  }
}

colSums(is.na(res_pre_based_on_1)) # no NA

to_drop_id <- c()

for (x in seq(1, length(rownames(res_pre_based_on_1)))) {
  if (!(res_pre_based_on_1[x, "X"] %in% needed_id2)) {
    to_drop_id <- c(to_drop_id, x)
  }
}

to_drop_id %>% length()

res_pre_based_on_1 <- res_pre_based_on_1[-to_drop_id, ]
rownames(res_pre_based_on_1) <- res_pre_based_on_1$X
# res_pre_based_on_1 <- res_pre_based_on_1[, -c(1)]

# real tree clade info
clade_info_2 <- read.csv('clade_info_2.csv', 
                         header = FALSE,
                         row.names = 1
)

cross_table <- res_pre_based_on_1

for (x in rownames(res_pre_based_on_1)) {
  cross_table[x, "X"] <- clade_info_2[x, "V2"]
}


CrossTable(cross_table$X, cross_table$clade)
