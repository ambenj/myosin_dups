## Prepare QTL plot
# Run with R 4.5.1

# load in library
library(qtl)


myh = read.cross(format = "csv", dir = ".", "MYH3QTL_numeric.csv", na.strings="-", genotypes=c("A", "H", "B"), estimate.map=FALSE)

# will give summary of all the data, make sure it looks good before moving forward
summary(myh)

# Estimate the sex-averaged recombination fraction between all pairs of genetic markers
est = est.rf(myh)

# setting step=0 means that you are using markers only.  If I set step=2 then it will estimate genotype probability every 2cM
genp = calc.genoprob(myh, step=0)

# genomic scan using single QTL model
out.norm=scanone(genp, pheno.col=1, model='normal', upper=TRUE)

sorted = out.norm[order(as.numeric(as.character(out.norm$chr))),]

plot(sorted)
# draw signficance line
abline(h=4.32, lty=5, col='red1')

#runs a permutation test to test for significance
operm <- scanone(genp, method="hk", n.perm=10000)

#calculates significance thresholds and p-values automatically
#MYH3 is the only region that comes up
summary(out.norm, perms=operm, alpha=0.01, pvalues=TRUE)

#finds LOD significance threshold
# pvalue<0.05 is LOD 3.55
summary(operm, alpha=c(0.01, 0.05, 0.2))
