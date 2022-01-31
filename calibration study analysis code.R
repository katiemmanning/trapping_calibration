##trap data analysis

#bring in data sets from github

pitfall <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Insect%20ID%202020_pitfall_order.csv",na.strings = NULL)

ramp <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Insect%20ID%202020_yellowramp_order.csv",na.strings = NULL)

jar <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Insect%20ID%202020_jarramp_order.csv",na.strings = NULL)

sticky <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Insect%20ID%202020_stickycard_order.csv",na.strings = NULL)

taxa <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Taxa_updated.csv")

#add trap type as a column on each data file
pitfall$Trap="pitfall"
ramp$Trap="ramp"
jar$Trap="jar"
sticky$Trap="sticky"

#combine data tables 
library (plyr)
pitfallramp <- rbind.fill (pitfall, ramp)
pitfallrampjar <-rbind.fill (pitfallramp, jar)
insects <- rbind.fill (pitfallrampjar, sticky)

#############
#NMDS of insect community between trap types
library (vegan)

#Create matrix of environmental variables
env.matrix<-insects[c(1:3,22)]
#create matrix of community variables
com.matrix<-insects[c(4:21)]

#ordination by NMDS
NMDS<-metaMDS(com.matrix, distance="bray", k=2, autotransform=FALSE, trymax=100)
stressplot(NMDS)

#NMDS visualization 

#what taxa to display using "taxa"
#most.abund<-as.vector(t(taxa[1,]))#greater than 100 caught
#most.abund<-most.abund[-1]
#bioind<-as.vector(t(taxa[2,]))
#bioind<-bioind[-1]
flying<-as.vector(t(taxa[1,]))
flying<-flying[-1]
crawling<-as.vector(t(taxa[2,]))
crawling<-crawling[-1]
include<-as.vector(t(taxa[3,]))
include<-include[-1]

#plot NMDS
plot(NMDS, disp='sites', type="n")
#add ellipsoids with ordiellipse
ordiellipse(NMDS, env.matrix$Trap, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "pitfall")
ordiellipse(NMDS, env.matrix$Trap, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "jar") 
ordiellipse(NMDS, env.matrix$Trap, draw="polygon", col="#F0E442",kind="sd", conf=0.95, label=FALSE, show.groups = "ramp") 
ordiellipse(NMDS, env.matrix$Trap, draw="polygon", col="#CC79A7",kind="sd", conf=0.95, label=FALSE, show.groups = "sticky")
#display ground trap data as solid shapes - pitfall=circle, ramp trap=square, jar=triangle, flying trap as triangle outline
points(NMDS, display="sites", select=which(env.matrix$Trap=="pitfall"),pch=19, col="#E69F00")
points(NMDS, display="sites", select=which(env.matrix$Trap=="jar"), pch=17, col="#009E73")
points(NMDS, display="sites", select=which(env.matrix$Trap=="ramp"), pch=15, col="#F0E442")
points(NMDS, display="sites", select=which(env.matrix$Trap=="sticky"), pch=25, col="#CC79A7")
#add legend
legend(1.46,1.45, title=NULL, pch=c(19,17,15,25), col=c("#E69F00","#009E73","#F0E442","#CC79A7"), cex=.7, legend=c("Pitfall", "Jar ramp", "Yellow ramp", "Yellow sticky card"))
#add insect taxa as text
ordilabel(NMDS, display="species", select =which (include==TRUE & crawling == TRUE), cex=0.6, col="black", fill="white")
ordilabel(NMDS, display="species", select =which (include==TRUE & flying == TRUE), cex=0.6, col="white", fill="black")

#bootstrapping and testing for differences between the groups (traps)
fit<-adonis(com.matrix ~ Trap, data = env.matrix, permutations = 999, method="bray")
fit

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix)
anova(betadisper(distances_data, env.matrix$Trap))
#P-value = 0.004 -- cannot assume homogeneity of multivariate dispersion


################
#calculate Abundance
insects.abun <- rowSums(insects[,4:21])
insects$abundance <- insects.abun

#calculate Richness
insects.rowsums <- rowSums(insects[,4:21]>0)
insects$richness <- insects.rowsums

#calculate Shannon diversity
diversity <-diversity(insects[,4:21])
insects$diversity <-diversity

#calculate Evenness
evenness <-diversity/log(specnumber(insects[,4:21]))
insects$evenness <- evenness

#######
#Mixed effects models
library(lme4)
library(lmerTest) #to obtain p values
library (emmeans) #for pairwise comparisons
library (multcompView) #to view letters

#richness
##AIC 646
richness.model<-lmer(richness ~ Trap + Date + (1 | Site), data=insects)
summary(richness.model)
anova(richness.model)
AIC(richness.model)
#pairwise comparison 
rich.emm<-emmeans(richness.model,pairwise~Trap)
rich.emm
#results: sig difference btw all
rich.cld<-multcomp::cld(rich.emm, alpha = 0.05, Letters = LETTERS)
rich.cld

#abundance
##AIC 1825
abundance.model<-lmer(abundance ~ Trap + Date + (1 | Site), data=insects)
summary(abundance.model)
anova(abundance.model)
AIC(abundance.model)
#pairwise comparison 
abun.emm<-emmeans(abundance.model,pairwise~Trap)
abun.emm
#results: no sig diff in abundance btw jar and pitfall (0.732); sig btw rest
abun.cld<-as.data.frame(multcomp::cld(abun.emm, alpha = 0.05, Letters = LETTERS))
abun.cld

#diversity
##AIC 82
diversity.model<-lmer(diversity ~ Trap + (1 | Site), data=insects)
summary(diversity.model)
anova(diversity.model)
AIC(diversity.model)
#pairwise comparison 
div.emm<-emmeans(diversity.model,pairwise~Trap)
div.emm
#results: no sig diff jar-pitfall (0.2652), jar-sticky (0.9288), pitfall-sticky (0.61); sig between rest
div.cld<-multcomp::cld(div.emm, alpha = 0.05, Letters = LETTERS)
div.cld

#evenness
##AIC -286
evenness.model<-lmer(evenness ~ Trap +(1 | Date), data=insects)
summary(evenness.model)
anova(evenness.model)
AIC(evenness.model)
#pairwise comparison 
even.emm<-emmeans(evenness.model,pairwise~Trap)
even.emm
#results: no sig diff between jar-ramp (0.1151), ramp-sticky (0.2065); sig btw rest
even.cld<-multcomp::cld(even.emm, alpha = 0.05, Letters = LETTERS)
even.cld

###########
library(ggplot2)
#abundance plot
abundance.plot<-ggplot(insects, aes(x =Trap, y = abundance, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Abundance (log10)")+
  scale_y_continuous(trans="log10")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=abun.cld, aes(y = 600, label = .group))
abundance.plot

#richness plot
richness.plot<-ggplot(insects, aes(x =Trap, y = richness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Richness")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=rich.cld, aes(y = 25, label = .group))
richness.plot

#diversity plot
diversity.plot<-ggplot(insects, aes(x =Trap, y = diversity, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Shannon Diversity")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=div.cld, aes(y = 2.5, label = .group))
diversity.plot

#evenness plot
evenness.plot<-ggplot(insects, aes(x =Trap, y = evenness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Evenness")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=even.cld, aes(y = 1.2, label = .group))
evenness.plot

#Mush plots together
library(ggpubr) 
figure3 <- ggarrange(richness.plot, abundance.plot, diversity.plot, evenness.plot,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "bottom")
pdf("figure3.pdf", height=6, width=8) #height and width in inches
figure3
dev.off()

figure3
#####
#functional group abundance by trap type

#input data
flying<-read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/flying.csv")
crawling<-read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/crawling.csv")

#calculating abundance for flying
flying.abun <- rowSums(flying[,2:32])
flying$abundance <- flying.abun

#calculating abundance for crawling
crawling.abun <- rowSums(crawling[,2:10])
crawling$abundance <- crawling.abun

#calculating richness for flying
flying.rich <- rowSums(flying[,2:32]>0)
flying$richness <- flying.rich

#calculating richness for crawling
crawling.rich <- rowSums(crawling[,2:10]>0)
crawling$richness <- crawling.rich

#abundance model for flying arthropods
#AIC = 1845
abundance.model_flying<-lm(abundance ~ Trap, data=flying)
summary(abundance.model_flying)
anova(abundance.model_flying)
AIC(abundance.model_flying)
#pairwise comparison
abun_f.emm<-emmeans(abundance.model_flying,pairwise~Trap)
abun_f.emm
#results: sig diff btw everything except jar-pitfall (0.9644) 
abun_f.cld<-multcomp::cld(abun_f.emm, alpha = 0.05, Letters = LETTERS)
abun_f.cld


#abundance model for crawling arthropods
#AIC = 1439
abundance.model_crawling<-lm(abundance ~ Trap, data=crawling)
summary(abundance.model_crawling)
anova(abundance.model_crawling)
AIC(abundance.model_crawling)
#pairwise comparison
abun_c.emm<-emmeans(abundance.model_crawling,pairwise~Trap)
abun_c.emm
#results: sig diff btw everything except jar-pitfall (0.8163) and pitfall-sticky (0.1255)
abun_c.cld<-multcomp::cld(abun_c.emm, alpha = 0.05, Letters = LETTERS)
abun_c.cld

#richness model for flying arthropods
#AIC = 668
richness.model_flying<-lm(richness ~ Trap, data=flying)
summary(richness.model_flying)
anova(richness.model_flying)
AIC(richness.model_flying)
#pairwise comparison
rich_f.emm<-emmeans(richness.model_flying,pairwise~Trap)
rich_f.emm
#results: sig diff btw everything except ramp-sticky (0.3286) and jar-pitfall(0.06)
rich_f.cld<-multcomp::cld(rich_f.emm, alpha = 0.05, Letters = LETTERS)
rich_f.cld

#richness model for crawling arthropods
#AIC = 452
richness.model_crawling<-glm(richness ~ Trap, data=crawling)
summary(richness.model_crawling)
anova(richness.model_crawling)
AIC(richness.model_crawling)
#pairwise comparison
rich_c.emm<-emmeans(richness.model_crawling,pairwise~Trap)
rich_c.emm
#results: sig diff btw everything except jar-pitfall (0.9814)
rich_c.cld<-multcomp::cld(rich_c.emm, alpha = 0.05, Letters = LETTERS)
rich_c.cld


##plot flying abundance
abundance.plot_flying<-ggplot(flying, aes(x =Trap, y = abundance, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(title="Flying", x="", y="Abundance (log10)")+
  scale_y_continuous(trans="log10")+
  theme (plot.title = element_text(hjust=0.5))+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=abun_f.cld, aes(y = 600, label = .group))
abundance.plot_flying

##plot crawling abundance
abundance.plot_crawling<-ggplot(crawling, aes(x =Trap, y = abundance, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(title="Ground-crawling", x="", y="")+
  scale_y_continuous(trans="log10")+
  theme (plot.title = element_text(hjust=0.5))+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=abun_c.cld, aes(y = 600, label = .group))
abundance.plot_crawling

##plot flying richness
richness.plot_flying<-ggplot(flying, aes(x =Trap, y = richness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(title="", x="", y="Richness")+
  theme (plot.title = element_text(hjust=0.5))+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=rich_f.cld, aes(y = 15, label = .group))
richness.plot_flying

##plot crawling richness
richness.plot_crawling<-ggplot(crawling, aes(x =Trap, y = richness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(title="", x="", y="")+
  theme (plot.title = element_text(hjust=0.5))+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=rich_c.cld, aes(y = 15, label = .group))
richness.plot_crawling

#mush together
figure4 <- ggarrange(abundance.plot_flying, abundance.plot_crawling,richness.plot_flying,richness.plot_crawling,
                     labels = c("A", "B", "C", "D"),
                     ncol = 2, nrow = 2,
                     common.legend = TRUE, legend = "bottom")
pdf("figure4.pdf", height=6, width=8) #height and width in inches
figure4
dev.off()

#####################
#species accumulation
library (BiodiversityR)
library(ggplot2)

#individual curves for each trap type
pitfall.com.matrix<-pitfall[c(4:52)]
pitfall_curve<-accumresult(pitfall.com.matrix, method = "exact", permutations = 1000)

jar.com.matrix<-jar[c(4:52)]
jar_curve<-accumresult(jar.com.matrix, method = "exact", permutations = 1000)

ramp.com.matrix<-ramp[c(4:52)]
ramp_curve<-accumresult(ramp.com.matrix, method = "exact", permutations = 1000)

sticky.com.matrix<-sticky[c(4:52)]
sticky_curve<-accumresult(sticky.com.matrix, method = "exact", permutations = 1000)

#first-order jackknife estimates are based on the number of singletons
#second-order jackknife estimates are based on the number of singletons and doubletons

#calculates species richness for each sample
specnumber(com.matrix) #ranges from 1 to 23

#calculates species richness by treatment (trap)
specnumber(com.matrix, groups = insects$Trap) #jar=32; pitfall=27; ramp=43; sticky=36

#total richness and jackknife
rich <- diversityresult(com.matrix, y=NULL, index = "richness")
rich # 49
j1 <- diversityresult(com.matrix, y=NULL, index = "jack1")
j1 # 55.957576
#88%
j2 <- diversityresult(com.matrix, y=NULL, index = "jack2")
j2 # 57.963452
#85%

#jar jackknife; richness = 32
j1.j <- diversityresult(jar.com.matrix, y=NULL, index = "jack1")
j1.j # 39.809524
#80%
j2.j <- diversityresult(jar.com.matrix, y=NULL, index = "jack2")
j2.j # 41.853659
#76%

#pitfall jackknife; richness = 27
j1.p <- diversityresult(pitfall.com.matrix, y=NULL, index = "jack1")
j1.p # 32.846154
#82%
j2.p <- diversityresult(pitfall.com.matrix, y=NULL, index = "jack2")
j2.p # 32.995951
#82%

#ramp jackknife; richness = 43
j1.r <- diversityresult(ramp.com.matrix, y=NULL, index = "jack1")
j1.r # 52.761905
#81%
j2.r <- diversityresult(ramp.com.matrix, y=NULL, index = "jack2")
j2.r # 60.42741
#71%

#sticky jackknife; richness = 36
j1.s <- diversityresult(sticky.com.matrix, y=NULL, index = "jack1")
j1.s # 42.833333
#84%
j2.s <- diversityresult(sticky.com.matrix, y=NULL, index = "jack2")
j2.s # 46.712544
#77%

#BiodiversityR::accumcomp
Accum.1 <- accumcomp(com.matrix, y=env.matrix, factor='Trap', 
                     method='random', conditioned=FALSE, plotit=FALSE)
Accum.1

#BiodiversityR::accumcomp.long
accum.long1 <- accumcomp.long(Accum.1, ci=NA, label.freq=5)
head(accum.long1)

#plot
#empty canvas
BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

figure5 <- ggplot(data=accum.long1, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_color_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  scale_shape_manual(values=c(19,17,15,25))+
  geom_line(aes(colour=Grouping), size=0.1) +
  geom_ribbon(aes(colour=Grouping, fill=after_scale(alpha(colour, 0.3))), 
              show.legend=FALSE, linetype = 0) + 
  geom_point(data=subset(accum.long1, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=3) +
  BioR.theme +
  labs(x = "Number of samples", y = "Richness", colour = "Trap", shape = "Trap")
figure5

pdf("Figure 5.pdf", height=6, width=8) #height and width in inches
figure5
dev.off()

