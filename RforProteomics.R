source("https://bioconductor.org/biocLite.R")
biocLite("RforProteomics", dependencies = TRUE)

library("RforProteomics")
library("RColorBrewer")
library("ggplot2")
library("reshape2")
library("mzR")
library("msdata")
library("rpx")
library("MALDIquant")
library("MALDIquantForeign")
library("IPPD")
library("BRAIN")
library(Rdisop)

#Experiment info
px1 <- PXDataset("PXD000001")
px1
pxfiles(px1)

#download mzTab data
mztab <- pxget(px1, "PXD000001_mztab.txt")
mztab

#Load mzTab peptide data
qnt <- readMzTabData(mztab, what = "PEP", version = "0.9")
sampleNames(qnt) <- reporterNames(TMT6)
head(exprs(qnt))

#removes missing values
qnt <- filterNA(qnt)
processingData(qnt)

#combine into protein using the accession feature meta data sum the peptide intensities
protqnt <- combineFeatures(qnt, groupBy = fData(qnt)$accession, fun = sum)

cls <- brewer.pal(5, "Set1")
matplot(t(tail(exprs(protqnt), n=5)), type = "b",
        lty=1, col=cls,
        ylab = "Protein intensity (summed peptides)",
        xlab = "TMT reporters")
legend("topright", tail(featureNames(protqnt), n=5), lty=1, bty="n", cex=0.8, col = cls)

qntS <- normalise(qnt, "sum")
qntV <- normalise(qntS, "vsn")
qntV2 <- normalise(qnt, "vsn")

acc <- c("P00489", "P00924", "P02769", "P62894", "ECA")

idx <- sapply(acc, grep, fData(qnt)$accession)
idx2 <- sapply(idx, head, 10)
small <- qntS[unlist(idx2),]

idx3 <- sapply(idx, head, 10)
medium <- qntV[unlist(idx3),]

m <- exprs(medium)
colnames(m) <- c("126", "127", "128", "129", "130", "131")
rownames(m) <- fData(medium)$accession
rownames(m)[grep("CYC", rownames(m))] <- "CYT"
rownames(m)[grep("ENO", rownames(m))] <- "ENO"
rownames(m)[grep("ALB", rownames(m))] <- "BSA"
rownames(m)[grep("PYGM", rownames(m))] <- "PHO"
rownames(m)[grep("ECA", rownames(m))] <- "Background"

cls <- c(brewer.pal(length(unique(rownames(m)))-1, "Set1"),"grey")
names(cls) <- unique(rownames(m))
wbcol <- colorRampPalette(c("white", "darkblue"))(256)         

heatmap(m, col=wbcol, RowSideColors = cls[rownames(m)])

#MALDIquant package
#Loading data
datapath <- file.path(system.file("Examples", package="readBrukerFlexData"),
                      "2010_05_19_Gibb_C8_A1")
dir(datapath)
sA1<- importBrukerFlex(datapath, verbose=F)
s <- sA1[[1]]
summary(mass(s))
summary(intensity(s))
head(as.matrix(s))
plot(s)

#preprocessing
#sqrt transform (variance stabilization)
s2<- transformIntensity(s, method="sqrt")
s2
#smooothing - 5 point moving average
s3<- smoothIntensity(s2, method="MovingAverage", halfWindowSize=2)
s3
#baseline subtraction
s4 <- removeBaseline(s3, method="SNIP")
s4

#Peak picking
p <- detectPeaks(s4)
length(p)
peak.data <- as.matrix(p)

par(mfrow=c(2,3))
xl <- range(mass(s))
# use same xlim on all plots for better comparison
plot(s, sub="", main="1: raw", xlim=xl)
plot(s2, sub="", main="2: variance stabilisation", xlim=xl)
plot(s3, sub="", main="3: smoothing", xlim=xl)
plot(s4, sub="", main="4: base line correction", xlim=xl)
plot(s4, sub="", main="5: peak detection", xlim=xl)
points(p)
top20 <- intensity(p) %in% sort(intensity(p), decreasing=TRUE)[1:20]
labelPeaks(p, index=top20, underline=TRUE)
plot(p, sub="", main="6: peak plot", xlim=xl)
labelPeaks(p, index=top20, underline=TRUE)

#Working with peptide sequences
atoms <- getAtomsFromSeq("SIVPSGASTGVHEALEMR")
unlist(atoms)

pepmol <- getMolecule(paste0(names(atoms),unlist(atoms),collapse = ""))
pepmol

library(OrgMassSpecR)
data(itraqdata)

simplottest <- itraqdata[featureNames(itraqdata) %in% paste0("X",46:47)]
par(mfrow=c(1,1))
sim <- SpectrumSimilarity(as(simplottest[[1]], "data.frame"),
                          as(simplottest[[2]], "data.frame"),
                          top.label = "itraqdata[['X46']]",
                          bottom.label = "itraqdata[['X47']]",
                          b=25)
title(main = paste("Spectrum similarity", round(sim, 3)))

MonoisotopicMass(formula = list(C=2,O=1, H=6))

molecule <- getMolecule("C2H5OH")
molecule$exactmass
