library(ape)
library(ggplot2)
library(ggtree)
library(dplyr)

# Reading trees from Newick format files
#ran_1 <- read.tree("/home/piyumal/test/iqtree_tree_3_partition.tree")
# ran_2 <- read.tree("/home/piyumal/test/mcmctree_3_partition.tree")


ran_1 <- read.tree("/home/piyumal/test/iqtree_xenarthra_test_xenarthra/tree_all_parts_after_br_fixed_HKY_G5_2/Xenarthra_5parts_iqtree.nex.rooted.mcmctree.tree")
ran_2 <- read.tree("/home/piyumal/test/iqtree_xenarthra_test_xenarthra/Xenarthra_old.tree")



ran_1 <- root(ran_1,  edgelabel = TRUE)
ran_2 <- root(ran_2,  edgelabel = TRUE)




d1 <- fortify(ran_1, branch.length='none')
d2 <- fortify(ran_2, branch.length='none')

d2$x <- max(d2$x) - d2$x + max(d1$x) + 4

dd = bind_rows(d1,d2) %>% 
  filter(!is.na(label)) 

# Note: you need to use dd, d1, d2 tables with filters to draw trees.

ggtree(ran_spe,branch.length='none') + 
  geom_tree() + 
  geom_tree(data = d2)+
  geom_tiplab(data=d1, size = 2) +
  #geom_tiplab(data=d2, size = 1) +
  geom_nodelab(size = 2, nudge_x = -0.5, nudge_y = 0.5)+
  geom_nodelab(data=d2, size = 2, nudge_x = 0.5, nudge_y = 0.5 )+
  geom_line(aes(x, y, group=label), data=dd[c(1:33,66:98),], alpha=0.3)
#geom_line(aes(x, y, group=label), data=dd[c(1:90,179:268),], alpha=0.3)