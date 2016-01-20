#required libraries
library(ggplot2)
library(reshape2)
library(plyr)

## 
cbPalette <- c("#56B4E9","#999999", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

customTheme <- list(theme(axis.title.x = element_text(face="bold", size=30), 
                          axis.title.y = element_text(face="bold", size=30), 
                          axis.text.x  = element_text(size=16),
                          axis.text.y  = element_text(size=16),
                          legend.title = element_text(size=30),
                          legend.text = element_text(size=20),
                          strip.text.x = element_text(size = 20),
                          strip.text.y = element_text(size = 20),
                          legend.position="none") )

#Read in datasets
spiral <- read.csv("../data/phyllotaxyScoring_11March2015.csv")
leaf <- read.csv("../data/complexityFull_15March2015.csv")

#Merge
full <- merge(spiral, leaf, by = "id")

head(full)
dim(full)

#There are four types of plants
levels(full$line) 
is.na(full$line)

#We need to get rid of the plants that full$id == Na.
fullsub1 <- subset(full, !is.na(full$line))
levels(fullsub1$line)

names(fullsub1)
head(fullsub1)

#PHYLLOTAXY
#Summarize
leaf.graph <-ddply(fullsub1, .(line),summarise,
                 nspiral= length(phenotype[phenotype=="spiral"]), 
                 nNoSpiral = length(phenotype[phenotype=="noSpiral"]), 
                 total = nspiral + nNoSpiral,
                 percentSpiral= nspiral/(nspiral+nNoSpiral)) 
#plot
head(leaf.graph)
leaf.graph$line <- gsub("PIN1::GFP", "e2/e2, +PIN1::GFP", leaf.graph$line)

ggplot(leaf.graph, aes(reorder(line, percentSpiral), percentSpiral)) +
  geom_bar(stat ="identity") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=30), 
        axis.title.y = element_text(face="bold", size=30), 
        axis.text.x  = element_text(size=16),
        axis.text.y  = element_text(size=16),
        legend.title = element_text(size=30),
        legend.text = element_text(size=20),
        strip.text.x = element_text(size = 20),
        legend.position="none") +
  ylab("percent spiral")  +
  xlab("genotype") 

#Just simplify.
head(leaf.graph)

leaf.graph.sub <- leaf.graph[c(-2),]

ggplot(leaf.graph.sub, aes(reorder(line, percentSpiral), percentSpiral)) +
  geom_bar(stat ="identity") +
  theme_bw() +
  theme(axis.title.x = element_text(face="bold", size=30), 
        axis.title.y = element_text(face="bold", size=30), 
        axis.text.x  = element_text(size=16),
        axis.text.y  = element_text(size=16),
        legend.title = element_text(size=30),
        legend.text = element_text(size=20),
        strip.text.x = element_text(size = 20),
        legend.position="none") +
  ylab("percent spiral")  +
  xlab("genotype") 

## Chi- squared
 
## Or Phenotype / Genotype (expected) not taking e2 or 3130 into account

#28 no spiral, 89 yes spiral
actual = c(28,89)    

#Genotype: 32 no PIN1GFP, 84 yes
expected = c(37,84)     

results = rbind(actual, expected)
head(results)
chisq.test(results)


#LEAF COMPLEXITY
#This is just useless it doesn't really help the cause because complexity needs to be sampled at many numbers
#So you really don't even see the difference between WT and e2.  Also, it appears as though the PIN1::GFP 
#background may have an effect on complexity.

head(fullsub1)

#summarize
leaf.graph <- ddply(fullsub1, .(line), summarise,
                    N  = length(total),
                    Leaflet = mean(total),
                    sd = sd(total),
                    se = sd / sqrt(N) )



ggplot(leaf.graph, aes(line, Leaflet)) +
   geom_bar(position = "dodge", stat = "identity") + 
  geom_errorbar(aes(ymin=Leaflet - se, ymax=Leaflet + se), 
                width=.2, 
                colour="black", 
                position = position_dodge(.9)) +
  theme_bw()


###Subsetted without parentals

sub <- subset(leaf.graph, line == "no PIN1::GFP" | line == "PIN1::GFP")

ggplot(sub, aes(line, Leaflet)) +
  geom_bar(position = "dodge", stat = "identity") + 
  geom_errorbar(aes(ymin=Leaflet - se, ymax=Leaflet + se), 
                width=.2, 
                colour="black", 
                position = position_dodge(.9)) +
                theme_bw() +
                customTheme


#The Primary Leaflet numbers might be messing this up
names(fullsub1)

fullsub1$int.secTotal <- fullsub1$intercalary + fullsub1$secondary
fullsub1 <- subset(fullsub1, !is.na(fullsub1$int.secTotal))

leaf.graph2 <- ddply(fullsub1, .(line), summarise,
                    N  = length(int.secTotal),
                    secintLeaflet = mean(int.secTotal),
                    sd = sd(int.secTotal),
                    se = sd / sqrt(N) )


ggplot(leaf.graph2, aes(line, secintLeaflet)) +
  geom_bar(position = "dodge", stat = "identity") + 
  geom_errorbar(aes(ymin=secintLeaflet - se, ymax=secintLeaflet + se), 
                width=.2, 
                colour="black", 
                position = position_dodge(.9)) +
  theme_bw()

#RACHIS
#remove NA
fullsub3 <- subset(fullsub1, !is.na(fullsub1$rachisLength))

#summarize
rachis.graph <- ddply(fullsub3, .(line), summarise,
                     N  = length(rachisLength),
                     rachis = mean(rachisLength),
                     sd = sd(rachisLength),
                     se = sd / sqrt(N) )

ggplot(rachis.graph, aes(line, rachis)) +
  geom_bar(position = "dodge", stat = "identity") + 
  geom_errorbar(aes(ymin=rachis - se, ymax=rachis + se), 
                width=.2, 
                colour="black", 
                position = position_dodge(.9)) +
  theme_bw()

#This is also useless. 

#TERMINAL LEAF LENGTH

#remove NA
fullsub4 <- subset(fullsub1, !is.na(fullsub1$terminalLeafLength))

#summarize
TL.graph <- ddply(fullsub4, .(line), summarise,
                      N  = length(terminalLeafLength),
                      TL = mean(terminalLeafLength),
                      sd = sd(terminalLeafLength),
                      se = sd / sqrt(N) )

ggplot(TL.graph, aes(line, TL)) +
  geom_bar(position = "dodge", stat = "identity") + 
  geom_errorbar(aes(ymin=TL - se, ymax=TL + se), 
                width=.2, 
                colour="black", 
                position = position_dodge(.9)) +
  theme_bw()

#LEAVES PRESENT

fullsub5 <- subset(fullsub1, !is.na(fullsub1$leavesPresent))

#summarize
LP.graph <- ddply(fullsub5, .(line), summarise,
                  N  = length(leavesPresent),
                  LP = mean(leavesPresent),
                  sd = sd(leavesPresent),
                  se = sd / sqrt(N) )

ggplot(LP.graph, aes(line, LP)) +
  geom_bar(position = "dodge", stat = "identity") + 
  geom_errorbar(aes(ymin=LP - se, ymax=LP + se), 
                width=.2, 
                colour="black", 
                position = position_dodge(.9)) +
  theme_bw()

