# First, read in your packages and your core alignment file

library(ape)
library(phangorn)

#read alignment
alignment <- read.dna("core_gene_alignment.aln", format = "fasta")
#your filepath here

#make matrix for NJ tree and ML tree
dist_matrix <- dist.dna(alignment, model = "JC69")  # Jukes-Cantor distance model

# First let's try an unrooted tree, which allows you to see how your sequences 
# are associated relative to each other, rather than to a common ancestor

# Read IQ-TREE .treefile (unrooted ML tree)
tree_unrooted <- read.tree("example.treefile")
plot(tree_unrooted, main = "Unrooted Tree", type = "unrooted")

# Now let's make a neighbor joining tree. Neighbor joining trees allow you to 
# examine distances between sequences in terms of how many base pairs are changed 
# between sequences, and builds a tree grouping the shortest possible distances

nj_tree <- nj(dist_matrix)
plot(nj_tree, main = "Neighbor Joining Tree")

# Now let's try a maximum likelihood tree. Maximum likelihood trees use different 
# modeling techniques, here the GTR, the General Time Reversible model which will 
# use the overall frequncies of the bases and the frequency that they vary 
# between sequences (e.g., all As, all Ts,  all switches from A-T or T-A, ) 

# Here we have the gamma parameter enabled which allows us to make heterogenous 
# rates of evolution between different areas (some regions are conserved, some 
# vary more) using a Gamma distribution.

phydat_alignment <- phyDat(alignment, type = "DNA")

# Start with NJ tree
start_tree <- nj(dist_matrix)

# Fit a model to examine rates of substitution (e.g. GTR + Gamma)
fit <- pml(start_tree, data = phydat_alignment)
fit <- optim.pml(fit, model = "GTR", optGamma = TRUE, control = pml.control(trace = 0))
fit$gamma  # here we can look at our different rates

# Plot ML tree
plot(fit$tree, main = "Maximum Likelihood Tree")

# bootstrap consensus tree: randomly resamples your alignment and builds trees 
# from each sample, then compares the trees and record the frequency of each 
# branch across all trees. then it builds consensus tree that will contain all 
# of the branches that occur in at least a given percentage of sample trees.

boot_tree <- bootstrap.phyDat(phydat_alignment, FUN = function(x) optim.pml(pml(nj_tree, x), model = "GTR")$tree, bs = 100)
cons_tree <- consensus(boot_tree, p = 0.5) 
# Here is where you decide what percentage of trees need to contain the branch 
# for it to appear in the consensus tree
?bootstrap.phyDat
plot(cons_tree, main = "Custom Bootstrap Consensus")

write.tree(nj_tree, file = "nj_tree.nwk")
write.tree(fit$tree, file = "ml_tree.nwk")
write.tree(boot_tree, file = "boot_trees.nwk")

# Now we can take these over to Microreact!