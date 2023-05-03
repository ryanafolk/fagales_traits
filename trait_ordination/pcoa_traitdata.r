library(ape)
library(phytools)

# Distance matrix should have same column and row number as taxon number in tree
distancematrix <- read.csv(file="distancematrix.csv", row.names = NULL, header = TRUE)

library(dplyr)
distancematrix %>% select(-contains(".1")) -> distancematrix
distancematrix %>% distinct(X, .keep_all = TRUE) -> distancematrix
setdiff(names(distancematrix), unique(names(distancematrix)))
row.names(distancematrix) <- distancematrix$X
distancematrix$X <- NULL

##############
# Regular PCoA
##############

# Replace NA with 0
distancematrix[is.na(distancematrix)] = 0

# Text-labeled plot
plot(cmdscale(distancematrix), xlab="Coordinate 1", ylab="Coordinate 2", main="MDS", pch = 19, cex = 0.1)
text(cmdscale(distancematrix), labels = row.names(distancematrix), cex = .1, pos = 3, offset = 0.1)

# Colored plot -- define colors by regex (make sure there are no nested matches in species names or change patterns accordingly
names <- row.names(distancematrix)
names <- gsub(" .*", "", names)
names <- gsub("Betulaceae", "blue", names)
names <- gsub("Betuloideae", "blue", names)
names <- gsub("Coryloideae", "blue", names)
names <- gsub("Alnus", "blue", names)
names <- gsub("Betula", "blue", names)
names <- gsub("Carpinus", "blue", names)
names <- gsub("Corylus", "blue", names)
names <- gsub("Ostrya", "blue", names)
names <- gsub("Ostryopsis", "blue", names)

names <- gsub("Casuarinaceae", "red", names)
names <- gsub("Casuarina", "red", names)
names <- gsub("Allocasuarina", "red", names)
names <- gsub("Ceuthostoma", "red", names)
names <- gsub("Gymnostoma", "red", names)

names <- gsub("Fagaceae", "forestgreen", names)
names <- gsub("Castanea", "forestgreen", names)
names <- gsub("Chrysolepis", "forestgreen", names)
names <- gsub("Fagus", "forestgreen", names)
names <- gsub("Lithocarpus", "forestgreen", names)
names <- gsub("Quercus", "forestgreen", names)
names <- gsub("Castanopsis", "forestgreen", names)
names <- gsub("Cyclobalanopsis", "forestgreen", names)
names <- gsub("Formanodendron", "forestgreen", names)


names <- gsub("Juglandaceae", "plum", names)
names <- gsub("Carya", "plum", names)
names <- gsub("Juglans", "plum", names)
names <- gsub("Annamocarya", "plum", names)
names <- gsub("Cyclocarya", "plum", names)
names <- gsub("Engelhardia", "plum", names)
names <- gsub("Platycarya", "plum", names)
names <- gsub("Pterocarya", "plum", names)

names <- gsub("Myricaceae", "darkgoldenrod1", names)
names <- gsub("Comptonia", "darkgoldenrod1", names)
names <- gsub("Myrica", "darkgoldenrod1", names)
names <- gsub("Canacomyrica", "darkgoldenrod1", names)

plot(cmdscale(distancematrix), xlab="Coordinate 1", ylab="Coordinate 2", main="MDS", col = names, pch = 19)
legend("bottomleft", title="Families", c("Betulaceae", "Casuarinaceae", "Fagaceae", "Juglandaceae", "Myricaceae"), fill = c("blue", "red", "forestgreen", "plum", "darkgoldenrod1"), horiz = TRUE, cex = 0.4)

mds <- cmdscale(distancematrix)
dataset <- read.csv(file="out_droppedmissing_coded_onevalue.csv", row.names = NULL, header = TRUE)
dataset %>% distinct(taxon, .keep_all = TRUE) -> dataset
row.names(dataset) <- dataset$taxon
dataset$taxon <- NULL

MDS1_loadings <- apply(dataset, 2, function(x) summary(lm(x ~ mds[,1]))$r.squared)
MDS2_loadings <- apply(dataset, 2, function(x) summary(lm(x ~ mds[,2]))$r.squared)

noquote(format(sort(MDS1_loadings, decreasing = TRUE), scientific = FALSE))
noquote(format(sort(MDS2_loadings, decreasing = TRUE), scientific = FALSE))


# Read tree
tree = read.newick("fagales.cut.rooted.tre")
# Remove long branches that will affect the plot
tree <- drop.tip(tree, c("Fagaceae_Quercus_tinkhamii", "Fagaceae_Quercus_leiophylla", "Fagaceae_Quercus_gambleana", "Fagaceae_Quercus_crispifolia", "Fagaceae_Quercus_neopalmeri", "Fagaceae_Lithocarpus_imperialis", "Fagaceae_Lithocarpus_sericobalanos", "Fagaceae_Lithocarpus_ruminatus", "Fagaceae_Lithocarpus_revolutus", "Fagaceae_Lithocarpus_hatusimae", "Fagaceae_Lithocarpus_echinifer", "Fagaceae_Lithocarpus_turbinatus", "Fagaceae_Lithocarpus_kalkmanii", "Fagaceae_Lithocarpus_beccarianus", "Fagaceae_Lithocarpus_truncatus", "Fagaceae_Colombobalanus_excelsa", "Fagaceae_Formanodendron_doichangensis"))

# Remove families from names and format
tree$tip.label <- gsub(".*aceae_", "", tree$tip.label)
tree$tip.label <- gsub("\\.", "_", tree$tip.label)
tree$tip.label <- gsub("__", "_", tree$tip.label)
tree$tip.label <- gsub("_$", "", tree$tip.label, perl=TRUE)


plotdata <- mds[,1] 
names(plotdata) <- gsub("\\.", "_", names(plotdata))
names(plotdata) <- gsub(" ", "_", names(plotdata))
names(plotdata) <- gsub("__", "_", names(plotdata))
names(plotdata) <- gsub("_$", "", names(plotdata), perl=TRUE)

tree.reduced <- treedata(tree, plotdata)$phy
temp <- treedata(tree, plotdata)$data
plotdata.reduced <- temp[,1]
names(plotdata.reduced) <- row.names(temp)

contMap(tree.reduced, plotdata.reduced, fsize = 0.1, lwd = 0.6, outline = FALSE)


## reduced taxa plot
#
#subset <- as.matrix(dist_subset(distance, c(grep(".*maulensis.*", colnames(distance), value = TRUE), grep(".*angustifolia.*", colnames(distance), value = TRUE), grep(".*scarlatina.*", colnames(distance), value = TRUE), grep(".*lutea.*", colnames(distance), value = TRUE), grep(".*angustifolia.*", colnames(distance), value = TRUE), grep(".*chilensis.*", colnames(distance), value = TRUE), grep(".*arzae.*", colnames(distance), value = TRUE), grep(".*australis.*", colnames(distance), value = TRUE), grep(".*fulgens.*", colnames(distance), value = TRUE), grep(".*davidii.*", colnames(distance), value = TRUE))))
#
#names <- row.names(subset)
#names <- gsub(".*amoena.*", "turquoise", names)
#names <- gsub(".*ornata.*", "olivedrab1", names)
#names <- gsub(".*lutea.*", "darkgoldenrod1", names)
#names <- gsub(".*scarlatina.*", "red", names)
#names <- gsub(".*davidii.*", "hotpink", names)
#names <- gsub(".*chilensis.*", "darkorchid3", names)
#names <- gsub(".*quilapilun.*", "steelblue1", names)
#names <- gsub(".*australis.*", "plum", names)
#names <- gsub(".*fulgens.*", "salmon", names)
#names <- gsub(".*maulensis.*", "forestgreen", names)
#names <- gsub(".*angustifolia.*", "darkorange1", names)
#names <- gsub(".*arzae.*", "royalblue1", names)
#
## Color plot
#plot(cmdscale(subset), xlab="Coordinate 1", ylab="Coordinate 2", main="MDS", col = names, pch = 19)
#
## Text-labeled plot
#plot(cmdscale(distance), xlab="Coordinate 1", ylab="Coordinate 2", main="MDS")
#text(cmdscale(distance), labels = row.names(distance), cex=.5, pos=3)


##############
# Phylogenetic PCoA
##############

# Make matrix names more regular
colnames(distancematrix) <- gsub("\\.", "_", colnames(distancematrix))
colnames(distancematrix) <- gsub("__", "_", colnames(distancematrix))
row.names(distancematrix) <- gsub("\\.", "_", row.names(distancematrix))
row.names(distancematrix) <- gsub(" ", "_", row.names(distancematrix))
row.names(distancematrix) <- gsub("__", "_", row.names(distancematrix))
row.names(distancematrix) <- gsub("_$", "", row.names(distancematrix), perl=TRUE)


# Match tree and data sampling
shared_list <- intersect(tree$tip.label, colnames(distancematrix))
tree.reduced <- drop.tip(tree,tree$tip.label[-match(shared_list, tree$tip.label)])
distancematrix.reduced <- distancematrix[shared_list,shared_list]
# Check that they match
setdiff(tree.reduced$tip.label,row.names(distancematrix.reduced))
setdiff(row.names(distancematrix.reduced), tree.reduced$tip.label)

# Run MDS
trait_mds = phyl.pca(tree.reduced, distancematrix.reduced)
trait_mds = data.frame(trait_mds$S)

PC1 = trait_mds$PC1
names(PC1) <- row.names(trait_mds)
PC1 = data.frame(PC1)
PC2 = trait_mds$PC2
names(PC2) <- row.names(trait_mds)
PC2 = data.frame(PC2)
PC3 = trait_mds$PC3
names(PC3) <- row.names(trait_mds)
PC3 = data.frame(PC3)
plot(trait_mds$PC1, trait_mds$PC2)
ggplot(trait_mds, aes(Re(PC1),Re(PC2))) + geom_point() + geom_text(aes(label=row.names(trait_mds)), size = 1, hjust = 0, nudge_x = 0.05)
ggplot(trait_mds, aes(Re(PC1),Re(PC2))) + geom_point() + geom_text(aes(label=row.names(trait_mds)), size = 1, hjust = 0, position=position_jitter(width=.1,height=.1))

plotdata <- Re(PC1$PC1) # Sometimes phyl.pca returns complex numbers, so we can just look at the real portion (in practice relative scaling doesn't differ)
names(plotdata) <- row.names(trait_mds)

contMap(tree.reduced, plotdata, fsize = 0.1, lwd = 0.6, outline = FALSE)


