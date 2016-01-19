#load packages
library(ggplot2)

#colorBlind color palettes

# The palette with grey:
cbPalette <- c("#999999","#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#load data
compData <- read.csv("../data/compData2.csv")

#get rid of NAs

compData <- compData[!is.na(compData$phyllotaxy1),]
compData <- na.omit(compData)

#check data
summary(compData)

head(compData)
summary(compData)


#Subsetting Data for genotypic specific questions

compdataE2 <- subset(compData, genotype == "e2", select = c("genotype", "phyllotaxy1"))
compdataWT <- subset(compData, genotype == "wt", select = c("genotype", "phyllotaxy1"))
compdataHET <- subset(compData, genotype == "het", select = c("genotype", "phyllotaxy1"))

#e2
summary(compdataE2)
length(compdataE2$phyllotaxy1)

#wt
summary(compdataWT)
length(compdataWT$phyllotaxy1)

#het
summary(compdataHET)
length(compdataHET$phyllotaxy1)

#How many have spiral phenotype
dim(subset(compData, genotype == "wt" & phyllotaxy1 == "spiral"))
dim(subset(compData, genotype == "e2" & phyllotaxy1 == "spiral"))
dim(subset(compData, genotype == "het" & phyllotaxy1 == "spiral"))

# Making the Graph
## Organize Data
dataPlot1 <- subset(compData, select = c("genotype", "phyllotaxy1"))

### Making New Column that is identical to phyllotaxy1
dataPlot1$phyllotaxy <- dataPlot1$phyllotaxy1

# plot
ggplot(dataPlot1, aes(x = genotype, fill = phyllotaxy)) + 
  geom_bar() +
  scale_fill_manual(values=cbPalette) +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=30), 
        axis.title.y = element_text(face="bold", size=30), 
        axis.text.x  = element_text(size=16),
        axis.text.y  = element_text(size=16),
        legend.title = element_text(size=30),
        legend.text = element_text(size=20),
        strip.text.x = element_text(size = 20),
        legend.position="none") 
        
