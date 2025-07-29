# read in packages
library(ape)
library(phangorn)

# read in alignment
alignment <- read.dna("core_gene_alignment.aln", format = "fasta")

# make matrix for NJ tree and ML tree
dist_matrix <- dist.dna(alignment, model = "JC69") # Jukes-Cantor distance model

# Start with a neighbor joining tree: allows us to examine distances between
# sequences in terms of how many base pairs are changed between sequences, and
# builds a tree grouping the shortest possible distances

nj_tree <- nj(dist_matrix)
plot(nj_tree, main = "Neighbor Joining Tree")

# Next, make a maximum likelihood tree: they use different modeling techniques
# (here the GTR, or General Time Reversible model) which uses the overall 
# frequencies of the bases and the frequency that they vary between sequences
# (e.g., all As, all Ts,  all switches from A-T or T-A, etc.) 

# Enabled the gamma parameter, which allows us to make heterogeneous rates of
# evolution between different areas (some regions are conserved, some 
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

# Next, we make a bootstrap consensus tree, which randomly resamples our 
# alignment and builds trees from each sample, then compares the trees and 
# records the frequency of each branch across all trees. then it builds a 
# consensus tree that will contain all of the branches that occur in at least 
# a given percentage of sample trees.
# the p refers to the percentages of trees neded to contain the branch for it 
# to appear in the consensus tree: here, we set it to 0.5

boot_tree <- bootstrap.phyDat(phydat_alignment, FUN = function(x) optim.pml(pml(nj_tree, x), model = "GTR")$tree, bs = 100)
cons_tree <- consensus(boot_tree, p = 0.5) 
plot(cons_tree, main = "Custom Bootstrap Consensus")

# write trees into .nwk files so they can be visualized on Microreact.org
write.tree(nj_tree, file = "nj_tree.nwk")
write.tree(fit$tree, file = "ml_tree.nwk")
write.tree(boot_tree, file = "boot_trees.nwk")

# can also make an unrooted tree -> allows us to see how our sequences 
# are associated relative to each other, rather than to a common ancestor

# to visualize the unrooted ML tree, for example:
tree_unrooted <- read.tree("ml_tree.nwk")
plot(tree_unrooted, main = "Unrooted Tree", type = "unrooted")