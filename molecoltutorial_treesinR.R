##############################################################
# code for basic analysis and estimation of phylogenies in R
# thanks to Arun Sethuraman for catching a few errors
##############################################################

#setwd()
install.packages("ape")
install.packages("phangorn")
library(ape)
library(phangorn)

#read in sequence data, convert to phyDat
mammals <- read.dna("mammals.dna", format="interleaved")
mammals_phyDat <- phyDat(mammals, type = "DNA", levels = NULL)
mammals10 <- subset(mammals_phyDat,1:10)
mammals10_phyDat <- phyDat(mammals10, type = "DNA", levels = NULL)

#model testing
mt <- modelTest(mammals10)
print(mt)

#estimate a distance matrix using a Jules-Cantor Model
dna_dist <- dist.ml(mammals10, model="JC69")

#quick and dirty UPGMA and NJ trees
mammals_UPGMA <- upgma(dna_dist)
mammals_NJ  <- NJ(dna_dist)
plot(mammals_UPGMA, main="UPGMA")
plot(mammals_NJ, main = "Neighbor Joining")

#parsimony searches 
mammals_optim <- optim.parsimony(mammals_NJ, mammals)
mammals_pratchet <- pratchet(mammals10) #returning error
plot(mammals_optim)
plot(mammals_pratchet)

#ml estimation w/ distance matrix
fit <- pml(mammals_NJ, mammals10)
print(fit)
fitJC <- optim.pml(fit, model = "JC", rearrangement = "stochastic")
logLik(fitJC)
bs <- bootstrap.pml(fitJC, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")

#subsetted alignment bs example
mammals10_dm <- dist.ml(mammals10)
mammals10_NJ  <- NJ(mammals10_dm)
fit2 <- pml(mammals10_NJ, data=mammals10)
print(fit2)
fitJC2 <- optim.pml(fit2, model = "JC", rearrangement = "stochastic")
logLik(fitJC2)
bs_subset <- bootstrap.pml(fitJC2, bs=100, optNni=TRUE, multicore=TRUE, control = pml.control(trace=0))
bs2 <- plotBS(midpoint(fitJC2$tree), bs, p = 50, type="p")

#exporting trees
write.tree(bs2, file="bootstrap_example.tre")

