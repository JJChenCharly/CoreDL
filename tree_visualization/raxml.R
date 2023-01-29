library(ape)
library(ggplot2)
library(ggtree)
library(treeio)
library(phytools)
library(r2r)
library(gmodels)

# getwd()

# functions ----
# search its nodes ----
# get every under ----
get_nodes_under <- function(node_No, the_ggplot, tip_count) {
  tip_node_under <- c(node_No)
  
  in_togo <- c(node_No)
  count_ <- 0
  in_togo2 <- c()
  
  while (length(in_togo) > 0) {
    
    if (count_ == 0) {
      in_togo <- c(node_No)
    } else {
      in_togo <- c()
      in_togo <- in_togo2
      in_togo2 <- c()
    }
    
    for (x in in_togo) {
      for (p in row.names(the_ggplot$data)) {
        if (the_ggplot$data[p, "parent"] == x) {
          tip_node_under <- c(tip_node_under, the_ggplot$data[p, "node"])
          
          if (the_ggplot$data[p, "node"] > tip_count) {
            in_togo2 <- c(in_togo2, the_ggplot$data[p, "node"]) 
          }
        }
      }
    }
    count_ <- count_ + 1
  }
  
  return(tip_node_under)
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
bootstrap_ <- read.raxml(file = "op_server/RAxML_bipartitionsBranchLabels.1135")
# bootstrap_ @phylo %>% str()

# str(bootstrap_)
# bootstrap_@data$bootstrap
# bootstrap_@data$node

# a quick look ----
#tree_1135 <- read.tree("91/RAxML_bestTree.91")# no bootstrap info
tree_1135 <- bootstrap_

# class(tree_1135)
# tree_1135[1]
# tree_1135

# clade color 
tip_color <- tree_1135@phylo$tip.label

for (x in seq(1, length(tip_color))) {
  if (grepl("YN-", tip_color[x])) {
    tip_color[x] <- "red"
  } else if (grepl("GCA", tip_color[x])) {
    tip_color[x] <- "black"
  } else {
    tip_color[x] <- "#0ABAB5"
  }
}

tree_plt <- ggtree(tree_1135,
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

to_drop <- c(to_drop, tree_plt$data[357, "label"]$label)

nodes_under_1497 <- get_nodes_under(node_No = 1497, 
                                    the_ggplot = tree_plt, 
                                    tip_count = 1135)
tips_under_1497 <- get_tip_under(char_lst_of_members = nodes_under_1497, 1135)


to_drop <- c(to_drop, tree_plt$data[tips_under_1497, "label"]$label)

to_drop %>% length()
to_drop
# reduce bootstrap
bootstrap_ <- drop.tip(bootstrap_, to_drop)

# second look ----
#tree_1135 <- read.tree("91/RAxML_bestTree.91")# no bootstrap info
tree_1127 <- bootstrap_

# class(tree_1127)
# tree_1127[1]
# tree_1127

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
                   # layout = "circular",
                   # open.angle = 90
) +
  #geom_label(aes(label = node),
  #           size = 2) +
  geom_text(#aes(label=ifelse(bootstrap < 95, bootstrap, NA)),
    aes(label=ifelse(bootstrap < 95, bootstrap, NA)),
    # aes(label=node),
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
tree_reroot <- root(tree_1127, 
                    #outgroup = c(247),
                    node = 1673,
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
                          # branch.length="none"
) +
  #geom_label(aes(label = node),
  #           size = 2) +
  geom_text(
    aes(label=ifelse(bootstrap < 95, bootstrap, NA)),
    # aes(label = node),
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

tree_plt_reroot$data[2252, ]# last row
tree_plt_reroot$data[2253, ]# empty row

# new node between 1488 1863 both non B2 of tree_plt_reroot

tree_plt_reroot$data[2253, "parent"] <- 1128# reduced tree has only 1127 tipes
tree_plt_reroot$data[2253, "node"] <- 2253
tree_plt_reroot$data[2253, "branch.length"] <- 0.002825
tree_plt_reroot$data[2253, "isTip"] <- FALSE
tree_plt_reroot$data[2253, "x"] <- 0 # temporary
tree_plt_reroot$data[2253, "y"] <- tree_plt_reroot$data[c(1488, 1863), "y"] %>% sum()/2
tree_plt_reroot$data[2253, "branch"] <- 0
tree_plt_reroot$data[2253, "angle"] <- 0

# change parent of 1488 1863
tree_plt_reroot$data[1488, "parent"] <- 2253
tree_plt_reroot$data[1863, "parent"] <- 2253

ggsave("temp_tree.pdf",
       tree_plt_reroot,
       height = 100,
       width = 15,
       limitsize = FALSE
)

nodes_under_2253 <- get_nodes_under(node_No = 2253, 
                                    the_ggplot = tree_plt_reroot,
                                    tip_count = 1127)
# nodes_under_2253 %>% length()# 319
# "368" %in% nodes_under2253

nodes_under_1674 <- get_nodes_under(node_No = 1674, 
                                    the_ggplot = tree_plt_reroot,
                                    1127)# B2 node
# nodes_under_1674[[1]][1] %>% typeof()

tree_plt_reroot$data[c(1127, 1128), ]

a <- get_tip_under(nodes_under_2253, 1127)
a %>% length()

b <- get_tip_under(nodes_under_1674, 1127)
b %>% length()

# clade position ----
tree_plt_reroot$data[1128, ] # x = 0, y = 262.
tree_plt_reroot$data[1674, ]
correction <- tree_plt_reroot$data[1674, "x"] %>% as.double()

for (x in as.character(nodes_under_2253)) {
  tree_plt_reroot$data[x, "x"] <- tree_plt_reroot$data[x, "x"] + (correction/2)
}

for (x in as.character(nodes_under_1674)) {
  tree_plt_reroot$data[x, "x"] <- tree_plt_reroot$data[x, "x"] - (correction/2)
}

ggsave("temp_tree.pdf",
       tree_plt_reroot 
       %>% flip(1674, 2253)
       %>% flip(1863, 1488)
       %>% flip(1297, 1299)
       %>% flip(1296, 2052)
       # %>% flip(1307, 1305)
       ,
       height = 80,
       width = 40,
       limitsize = FALSE
)

tree_plt_reroot$data[2253, ]

# draw ----
aes_df <- data.frame(id=c(1674, 1863, 1299, 2052, 1296, 1489),
                     annote=c("B2", "D", "E", "A", "B1", "N9"),
                     annote_no = c("", "", "", "", "", ""),
                     colr=c("red", "purple", "purple", "purple", "purple", "cyan"),
                     offset_1=c(0.01, 0.01, 0.01, 0.01, 0.01, 0.01),
                     offset_2=rep(0.0135, 6),#0.0135
                     offset_text=rep(-0.003, 6)
)

flipped_tree <- tree_plt_reroot %>% flip(1674, 2253) %>% flip(1863, 1488) %>% flip(1297, 1299) %>% flip(1296, 2052) +
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
    node = c(1674, 1863, 1299, 2052, 1296, 1489),
    # label = c("B2", "D", "E", "A", "B1", "N9"),
    label = rep("", 6),
    offset = rep(0.016, 6),#0.0113
    barcolor = c("red", "purple", "purple", "purple", "purple", "cyan"),
    barsize = 80,
    extend = rep(0.4, 6),
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
    node = c(1674, 1863, 1299, 2052, 1296, 1489),
    label = c("B2", "D", "E", "A", "B1", "N9"),
    offset = rep(0, 6),
    offset.text =rep(0.0145, 6),
    textcolour = "black",
    fontsize = 50, # size of clade bar txt
    align = TRUE,
    barcolor = "white",
    hjust = rep(-0.01, 6),
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
  geom_hilight(node=c(1674, 1863, 1299, 2052, 1296, 1489), 
               # fill=c("red", "purple", "cyan"), 
               fill=c("red", "#CA6F1E", "#1ABC9C", "#8E44AD", "#F4D03F", "cyan"),
               alpha=.1, 
               # extend = c(0.0205, 0, 0, 0, 0, 0.0203),
               extend = c(0.024, 0.0265, 0.0068, 0.0052, 0.0042, 0),
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

# member count ----
for (x in c(1674, 1863, 1299, 2052, 1296, 1489)) {
  nodes <- get_nodes_under(x, flipped_tree, 1127)
  tips <- get_tip_under(nodes, 1127)
  
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
190*3 + 183 + 188 + 186 == 1127 # TRUE
# clade info save ----
df <- data.frame(node_no = c(1674, 1863, 1299, 2052, 1296, 1489))
rownames(df) <- c("B2", "D", "E", "A", "B1", "N9")

clade_info <- data.frame(strain = c(),
                         clade = c())
for (x in rownames(df)) {
  # print(typeof(x))
  ndoes <- get_nodes_under(df[x, "node_no"], 
                           flipped_tree,
                           1127)
  tips <- get_tip_under(ndoes, 1127) %>% as.vector()
  
  for (y in tips) {
    clade_info <- rbind(clade_info, data.frame(strain = c(flipped_tree$data$label[y]),
                                               clade = c(x)
    )
    )
  }
}

write.table(clade_info,
            "clade_info_2.csv",
            sep = ",",
            quote = FALSE,
            row.names = FALSE,
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
