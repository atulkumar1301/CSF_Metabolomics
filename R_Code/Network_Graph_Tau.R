library(data.table)
library(igraph)

####Increased Metabolites
Nodes <- fread ("TauPET_Nodes_Positive.txt")
Edges <- fread ("TauPET_Edges_Positive.txt")
net.i <- graph_from_data_frame(d=Edges, vertices=Nodes, directed=T)
#l <- layout_with_fr(net)
#l <- layout_with_lgl(net)
E(net.i)$width <- E(net.i)$Weight
colrs.i <- c ("red", "gold")
V(net.i)$color <- colrs.i[V(net.i)$Type]
l.i <- layout_with_kk(net.i)
par(family = "Times")
plot(net.i, vertex.label=V(net.i)$Nodes, layout = l.i, vertex.label.family = "Times", vertex.label.dist	= -pi/2, vertex.label.color="black", vertex.label.font	= 4,
     main = expression ("Increased Metabolites associated with TauPET with and without adjustment"))

####Decreased Metabolites
Nodes <- fread ("TauPET_Nodes_Negative.txt")
Edges <- fread ("TauPET_Edges_Negative.txt")
net.d <- graph_from_data_frame(d=Edges, vertices=Nodes, directed=T)
#l <- layout_with_fr(net)
#l <- layout_with_lgl(net)
E(net.d)$width <- E(net.d)$Weight
colrs.d <- c ("red", "lightsteelblue2")
V(net.d)$color <- colrs.d[V(net.d)$Type]
l.d <- layout_with_kk(net.d)
par(family = "Times")
plot(net.d, vertex.label=V(net.d)$Nodes, layout = l.d, vertex.label.family = "Times", vertex.label.dist	= -pi/2, vertex.label.color="black", vertex.label.font	= 4,
     main = expression ("Decreased Metabolites associated with TauPET with and without adjustment"))
