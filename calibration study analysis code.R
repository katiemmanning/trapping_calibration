##trap data analysis

#bring in order data sets from github

pitfall_order <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Insect%20ID%202020_pitfall_order.csv",na.strings = NULL)

ramp_order <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Insect%20ID%202020_yellowramp_order.csv",na.strings = NULL)

jar_order <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Insect%20ID%202020_jarramp_order.csv",na.strings = NULL)

sticky_order <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Insect%20ID%202020_stickycard_order.csv",na.strings = NULL)

taxa_order <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Order%20taxa.csv")

#add trap type as a column on each data file
pitfall_order$Trap="pitfall"
ramp_order$Trap="ramp"
jar_order$Trap="jar"
sticky_order$Trap="sticky"

#combine order data tables 
library (plyr)
pitfallramp_order <- rbind.fill (pitfall_order, ramp_order)
pitfallrampjar_order <-rbind.fill (pitfallramp_order, jar_order)
insects_order <- rbind.fill (pitfallrampjar_order, sticky_order)

#############
#NMDS of insect community by order between trap types
library (vegan)

#Create matrix of environmental variables
env.matrix_order<-insects_order[c(1:3,16)]
#create matrix of community variables
com.matrix_order<-insects_order[c(4:15)]

#ordination by NMDS
NMDS_order<-metaMDS(com.matrix_order, distance="bray", k=2, autotransform=FALSE, trymax=100)
stressplot(NMDS_order)
#stress=0.14

#order NMDS visualization 

#what taxa to display using "taxa"
flying<-as.vector(t(taxa_order[1,]))
flying<-flying[-1]
crawling<-as.vector(t(taxa_order[2,]))
crawling<-crawling[-1]
include<-as.vector(t(taxa_order[3,]))
include<-include[-1]

#plot order NMDS
#10x12
plot(NMDS_order, disp='sites', type="n")
#add ellipsoids with ordiellipse
ordiellipse(NMDS_order, env.matrix_order$Trap, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "pitfall")
ordiellipse(NMDS_order, env.matrix_order$Trap, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "jar") 
ordiellipse(NMDS_order, env.matrix_order$Trap, draw="polygon", col="#F0E442",kind="sd", conf=0.95, label=FALSE, show.groups = "ramp") 
ordiellipse(NMDS_order, env.matrix_order$Trap, draw="polygon", col="#CC79A7",kind="sd", conf=0.95, label=FALSE, show.groups = "sticky")
#display ground trap data as solid shapes - pitfall=circle, ramp trap=square, jar=triangle, flying trap as triangle outline
points(NMDS_order, display="sites", select=which(env.matrix_order$Trap=="pitfall"),pch=19, col="#E69F00")
points(NMDS_order, display="sites", select=which(env.matrix_order$Trap=="jar"), pch=17, col="#009E73")
points(NMDS_order, display="sites", select=which(env.matrix_order$Trap=="ramp"), pch=15, col="#F0E442")
points(NMDS_order, display="sites", select=which(env.matrix_order$Trap=="sticky"), pch=25, col="#CC79A7")
#add legend
legend(1.20,1.35, title=NULL, pch=c(19,17,15,25), col=c("#E69F00","#009E73","#F0E442","#CC79A7"), cex=.7, legend=c("Pitfall", "Jar ramp", "Yellow ramp", "Yellow sticky card"))
#add insect taxa as text
ordilabel(NMDS_order, display="species", select =which (include==TRUE & crawling == TRUE), cex=0.6, col="black", fill="white")
ordilabel(NMDS_order, display="species", select =which (include==TRUE & flying == TRUE), cex=0.6, col="white", fill="black")

#bootstrapping and testing for differences between the groups (traps)
fit<-adonis(com.matrix_order ~ Trap, data = env.matrix_order, permutations = 999, method="bray")
fit
#P-value = 0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix_order)
anova(betadisper(distances_data, env.matrix_order$Trap))
#P-value = 0.005918 -- cannot assume homogeneity of multivariate dispersion


################
#calculate order Abundance
insects.abun_order <- rowSums(insects_order[,4:15])
insects_order$abundance <- insects.abun_order

#calculate order Richness
insects.rowsums_order <- rowSums(insects_order[,4:15]>0)
insects_order$richness <- insects.rowsums_order

#calculate order Shannon diversity
diversity_order <-diversity(insects_order[,4:15])
insects_order$diversity <-diversity_order

#calculate order Evenness
evenness_order <-diversity_order/log(specnumber(insects_order[,4:15]))
insects_order$evenness <- evenness_order

#######
#Mixed effects models
library(lme4)
library(lmerTest) #to obtain p values
library (emmeans) #for pairwise comparisons
library (multcompView) #to view letters

#order richness
##AIC 567
richness.model_order<-lmer(richness ~ Trap + Date + (1 | Site), data=insects_order)
summary(richness.model_order)
anova(richness.model_order)
AIC(richness.model_order)
#pairwise comparison 
rich.emm_order<-emmeans(richness.model_order,pairwise~Trap)
rich.emm_order
#results: sig difference btw all
rich.cld_order<-multcomp::cld(rich.emm_order, alpha = 0.05, Letters = LETTERS)
rich.cld_order

#order abundance
##AIC 1794
abundance.model_order<-lmer(abundance ~ Trap + Date + (1 | Site), data=insects_order)
summary(abundance.model_order)
anova(abundance.model_order)
AIC(abundance.model_order)
#pairwise comparison 
abun.emm_order<-emmeans(abundance.model_order,pairwise~Trap)
abun.emm_order
#results: no sig diff in abundance btw jar and pitfall (0.8055); sig btw rest
abun.cld_order<-multcomp::cld(abun.emm_order, alpha = 0.05, Letters = LETTERS)
abun.cld_order

#order diversity
##AIC 128 (103 w/o Date)
diversity.model_order<-lmer(diversity ~ Trap + Date + (1 | Site), data=insects_order)
summary(diversity.model_order)
anova(diversity.model_order)
AIC(diversity.model_order)
#pairwise comparison 
div.emm_order<-emmeans(diversity.model_order,pairwise~Trap)
div.emm_order
#results: no sig diff jar-pitfall (0.4304), jar-sticky (0.8797), pitfall-sticky (0.1163); sig between rest
div.cld_order<-multcomp::cld(div.emm_order, alpha = 0.05, Letters = LETTERS)
div.cld_order

#order evenness
##AIC -186 (-206 w/o Date)
evenness.model_order<-lmer(evenness ~ Trap + Date + (1 | Site), data=insects_order)
summary(evenness.model_order)
anova(evenness.model_order)
AIC(evenness.model_order)
#pairwise comparison 
even.emm_order<-emmeans(evenness.model_order,pairwise~Trap)
even.emm_order
#results: no sig diff between jar-pitfall (0.1139), jar-ramp (0.8769),jar-sticky (0.0743), ramp-sticky (0.3328); sig btw rest
even.cld_order<-multcomp::cld(even.emm_order, alpha = 0.05, Letters = LETTERS)
even.cld_order

###########
library(ggplot2)
#order abundance plot
abundance.plot_order<-ggplot(insects_order, aes(x =Trap, y = abundance, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Abundance (log10)")+
  scale_y_continuous(trans="log10")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=abun.cld_order, aes(y = 600, label = .group))
abundance.plot_order

#order richness plot
richness.plot_order<-ggplot(insects_order, aes(x =Trap, y = richness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Richness")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=rich.cld_order, aes(y = 25, label = .group))
richness.plot_order

#order diversity plot
diversity.plot_order<-ggplot(insects_order, aes(x =Trap, y = diversity, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Shannon Diversity")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=div.cld_order, aes(y = 2.5, label = .group))
diversity.plot_order

#order evenness plot
evenness.plot_order<-ggplot(insects_order, aes(x =Trap, y = evenness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="Evenness")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=even.cld_order, aes(y = 1.2, label = .group))
evenness.plot_order

#Mush order plots together
library(ggpubr) 
orderfigure <- ggarrange(richness.plot_order, abundance.plot_order, diversity.plot_order, evenness.plot_order,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "bottom")
pdf("order.pdf", height=6, width=8) #height and width in inches
orderfigure
dev.off()

orderfigure
#####
#Cannot do functional group abundance by trap type for order
######

#species accumulation for order
library (BiodiversityR)
library(ggplot2)

#individual curves for each trap type
pitfall.com.matrix_order<-pitfall_order[c(4:15)]
pitfall_curve_order<-accumresult(pitfall.com.matrix_order, method = "exact", permutations = 1000)

jar.com.matrix_order<-jar_order[c(4:15)]
jar_curve_order<-accumresult(jar.com.matrix_order, method = "exact", permutations = 1000)

ramp.com.matrix_order<-ramp_order[c(4:15)]
ramp_curve_order<-accumresult(ramp.com.matrix_order, method = "exact", permutations = 1000)

sticky.com.matrix_order<-sticky_order[c(4:15)]
sticky_curve_order<-accumresult(sticky.com.matrix_order, method = "exact", permutations = 1000)

#first-order jackknife estimates are based on the number of singletons
#second-order jackknife estimates are based on the number of singletons and doubletons

#calculates order richness for each sample
specnumber(com.matrix_order) #ranges from 1 to 10

#calculates order richness by treatment (trap)
specnumber(com.matrix_order, groups = insects_order$Trap) #jar=12; pitfall=9; ramp=12; sticky=10

#total richness and jackknife
rich <- diversityresult(com.matrix_order, y=NULL, index = "richness")
rich # 12
j1 <- diversityresult(com.matrix_order, y=NULL, index = "jack1")
j1 # 12
#100%
j2 <- diversityresult(com.matrix_order, y=NULL, index = "jack2")
j2 # 12
#100%

#jar jackknife; richness = 12
j1.j <- diversityresult(jar.com.matrix_order, y=NULL, index = "jack1")
j1.j # 13.952381
#86%
j2.j <- diversityresult(jar.com.matrix_order, y=NULL, index = "jack2")
j2.j # 14.927991
#80%

#pitfall jackknife; richness = 9
j1.p <- diversityresult(pitfall.com.matrix_order, y=NULL, index = "jack1")
j1.p # 9.974359
#90%
j2.p <- diversityresult(pitfall.com.matrix_order, y=NULL, index = "jack2")
j2.p # 10.923077
#82%

#ramp jackknife; richness = 12
j1.r <- diversityresult(ramp.com.matrix_order, y=NULL, index = "jack1")
j1.r # 12
#100%
j2.r <- diversityresult(ramp.com.matrix_order, y=NULL, index = "jack2")
j2.r # 11.070848
#108% --> 100%

#sticky jackknife; richness = 10
j1.s <- diversityresult(sticky.com.matrix_order, y=NULL, index = "jack1")
j1.s # 10
#100%
j2.s <- diversityresult(sticky.com.matrix_order, y=NULL, index = "jack2")
j2.s # 10
#100%

#BiodiversityR::accumcomp
Accum.1_order <- accumcomp(com.matrix_order, y=env.matrix_order, factor='Trap', 
                     method='random', conditioned=FALSE, plotit=FALSE)
Accum.1_order

#BiodiversityR::accumcomp.long
accum.long1_order <- accumcomp.long(Accum.1_order, ci=NA, label.freq=5)
head(accum.long1_order)

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

order_accum <- ggplot(data=accum.long1_order, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_color_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  scale_shape_manual(values=c(19,17,15,25))+
  geom_line(aes(colour=Grouping), size=0.1) +
  geom_ribbon(aes(colour=Grouping, fill=after_scale(alpha(colour, 0.3))), 
              show.legend=FALSE, linetype = 0) + 
  geom_point(data=subset(accum.long1_order, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=3) +
  BioR.theme +
  labs(x = "", y = "Richness", colour = "Trap", shape = "Trap")
order_accum

pdf("order_accum.pdf", height=6, width=8) #height and width in inches
order_accum
dev.off()


############################################################################
#bring in functional data sets from github

pitfall <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Insect%20ID%202020_pitfall_functional.csv",na.strings = NULL)

ramp <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Insect%20ID%202020_yellowramp_functional.csv",na.strings = NULL)

jar <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Insect%20ID%202020_jarramp_functional.csv",na.strings = NULL)

sticky <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Insect%20ID%202020_stickycard_functional.csv",na.strings = NULL)

taxa <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/Functional%20taxa.csv")

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
#NMDS of insect community by functional classification between trap types
library (vegan)

#Create matrix of environmental variables
env.matrix<-insects[c(1:3,43)]
#create matrix of community variables
com.matrix<-insects[c(4:42)]

#ordination by NMDS
NMDS<-metaMDS(com.matrix, distance="bray", k=2, autotransform=FALSE, trymax=100)
stressplot(NMDS)
#stress=0.15

#functional classification NMDS visualization 

#what taxa to display using "taxa"
flying<-as.vector(t(taxa[1,]))
flying<-flying[-1]
crawling<-as.vector(t(taxa[2,]))
crawling<-crawling[-1]
include<-as.vector(t(taxa[3,]))
include<-include[-1]

#plot functional NMDS
#10x12
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
#legend(1.0,1.51, title=NULL, pch=c(19,17,15,25), col=c("#E69F00","#009E73","#F0E442","#CC79A7"), cex=.7, legend=c("Pitfall", "Jar ramp", "Yellow ramp", "Yellow sticky card"))
#add insect taxa as text
ordilabel(NMDS, display="species", select =which (include==TRUE & crawling == TRUE), cex=0.6, col="black", fill="white")
ordilabel(NMDS, display="species", select =which (include==TRUE & flying == TRUE), cex=0.6, col="white", fill="black")

#bootstrapping and testing for differences between the groups (traps)
fit<-adonis(com.matrix ~ Trap, data = env.matrix, permutations = 999, method="bray")
fit
#P-value = 0.001

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix)
anova(betadisper(distances_data, env.matrix$Trap))
#P-value = .0001 -- cannot assume homogeneity of multivariate dispersion


################
#calculate Abundance
insects.abun <- rowSums(insects[,4:42])
insects$abundance <- insects.abun

#calculate Richness
insects.rowsums <- rowSums(insects[,4:42]>0)
insects$richness <- insects.rowsums

#calculate Shannon diversity
diversity <-diversity(insects[,4:42])
insects$diversity <-diversity

#calculate Evenness
evenness <-diversity/log(specnumber(insects[,4:42]))
insects$evenness <- evenness

#######
#Mixed effects models
library(lme4)
library(lmerTest) #to obtain p values
library (emmeans) #for pairwise comparisons
library (multcompView) #to view letters

#richness
##AIC 718
richness.model<-lmer(richness ~ Trap + Date + (1 | Site), data=insects)
summary(richness.model)
anova(richness.model)
AIC(richness.model)
#pairwise comparison 
rich.emm<-emmeans(richness.model,pairwise~Trap)
rich.emm
#results: jar-pitfall no sig diff (0.0610), sig dif btw all others
rich.cld<-multcomp::cld(rich.emm, alpha = 0.05, Letters = LETTERS)
rich.cld

#abundance
##AIC 1795
abundance.model<-lmer(abundance ~ Trap + Date + (1 | Site), data=insects)
summary(abundance.model)
anova(abundance.model)
AIC(abundance.model)
#pairwise comparison 
abun.emm<-emmeans(abundance.model,pairwise~Trap)
abun.emm
#results: jar-pitfall no sig diff (0.8089), sig dif btw all others
abun.cld<-multcomp::cld(abun.emm, alpha = 0.05, Letters = LETTERS)
abun.cld

#diversity
##AIC 175 (152 w/o Date)
diversity.model<-lmer(diversity ~ Trap + Date + (1 | Site), data=insects)
summary(diversity.model)
anova(diversity.model)
AIC(diversity.model)
#pairwise comparison 
div.emm<-emmeans(diversity.model,pairwise~Trap)
div.emm
#results: no sig diff btw jar-pitfall (0.2016), jar-sticky (0.9540), pitfall-sticky (0.0661); sig diff btw all others 
div.cld<-multcomp::cld(div.emm, alpha = 0.05, Letters = LETTERS)
div.cld

#evenness
##AIC -172 (-193 w/o Date)
evenness.model<-lmer(evenness ~ Trap + Date + (1 | Site), data=insects)
summary(evenness.model)
anova(evenness.model)
AIC(evenness.model)
#pairwise comparison 
even.emm<-emmeans(evenness.model,pairwise~Trap)
even.emm
#results: no sig diff btw jar-pitfall (0.2851) or ramp-sticky (0.0974); sig diff btw all others
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
  labs(x="", y="")+
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
  labs(x="", y="")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=rich.cld, aes(y = 25, label = .group))
richness.plot

#diversity plot
diversity.plot<-ggplot(insects, aes(x =Trap, y = diversity, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=div.cld, aes(y = 2.5, label = .group))
diversity.plot

#evenness plot
evenness.plot<-ggplot(insects, aes(x =Trap, y = evenness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=even.cld, aes(y = 1.2, label = .group))
evenness.plot

#Mush order plots together
library(ggpubr) 
functionalfigure <- ggarrange(richness.plot, abundance.plot, diversity.plot, evenness.plot,
                         labels = c("E", "F", "G", "H"),
                         ncol = 2, nrow = 2,
                         common.legend = TRUE, legend = "bottom")
pdf("functional.pdf", height=6, width=8) #height and width in inches
functionalfigure
dev.off()

functionalfigure
###

#flying vs crawling
#input data
flying<-read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/flying.csv")
crawling<-read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/crawling.csv")

#calculating abundance for flying
flying.abun <- rowSums(flying[,2:30])
flying$abundance <- flying.abun

#calculating abundance for crawling
crawling.abun <- rowSums(crawling[,2:9])
crawling$abundance <- crawling.abun

#calculating richness for flying
flying.rich <- rowSums(flying[,2:30]>0)
flying$richness <- flying.rich

#calculating richness for crawling
crawling.rich <- rowSums(crawling[,2:9]>0)
crawling$richness <- crawling.rich

#abundance model for flying arthropods
#AIC = 1855
abundance.model_flying<-lm(abundance ~ Trap, data=flying)
summary(abundance.model_flying)
anova(abundance.model_flying)
AIC(abundance.model_flying)
#pairwise comparison
abun_f.emm<-emmeans(abundance.model_flying,pairwise~Trap)
abun_f.emm
abun_f.cld<-multcomp::cld(abun_f.emm, alpha = 0.05, Letters = LETTERS)
abun_f.cld

#abundance model for crawling arthropods
#AIC = 1443
abundance.model_crawling<-lm(abundance ~ Trap, data=crawling)
summary(abundance.model_crawling)
anova(abundance.model_crawling)
AIC(abundance.model_crawling)
#pairwise comparison
abun_c.emm<-emmeans(abundance.model_crawling,pairwise~Trap)
abun_c.emm
abun_c.cld<-multcomp::cld(abun_c.emm, alpha = 0.05, Letters = LETTERS)
abun_c.cld

#richness model for flying arthropods
#AIC = 662
richness.model_flying<-lm(richness ~ Trap, data=flying)
summary(richness.model_flying)
anova(richness.model_flying)
AIC(richness.model_flying)
#pairwise comparison
rich_f.emm<-emmeans(richness.model_flying,pairwise~Trap)
rich_f.emm
rich_f.cld<-multcomp::cld(rich_f.emm, alpha = 0.05, Letters = LETTERS)
rich_f.cld

#richness model for crawling arthropods
#AIC = 488
richness.model_crawling<-lm(richness ~ Trap, data=crawling)
summary(richness.model_crawling)
anova(richness.model_crawling)
AIC(richness.model_crawling)
#pairwise comparison
rich_c.emm<-emmeans(richness.model_crawling,pairwise~Trap)
rich_c.emm
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
figure4
pdf("Figure 4.pdf", height=6, width=8) #height and width in inches
figure4
dev.off()

#####################
#species accumulation
library (BiodiversityR)
library(ggplot2)

#individual curves for each trap type
pitfall.com.matrix<-pitfall[c(4:42)]
pitfall_curve<-accumresult(pitfall.com.matrix, method = "exact", permutations = 1000)

jar.com.matrix<-jar[c(4:42)]
jar_curve<-accumresult(jar.com.matrix, method = "exact", permutations = 1000)

ramp.com.matrix<-ramp[c(4:42)]
ramp_curve<-accumresult(ramp.com.matrix, method = "exact", permutations = 1000)

sticky.com.matrix<-sticky[c(4:42)]
sticky_curve<-accumresult(sticky.com.matrix, method = "exact", permutations = 1000)

#first-order jackknife estimates are based on the number of singletons
#second-order jackknife estimates are based on the number of singletons and doubletons

#calculates species richness for each sample
specnumber(com.matrix) #ranges from 1 to 20

#calculates species richness by treatment (trap)
specnumber(com.matrix, groups = insects$Trap) #jar=26; pitfall=21; ramp=35; sticky=31

#total richness and jackknife
rich <- diversityresult(com.matrix, y=NULL, index = "richness")
rich # 39
j1 <- diversityresult(com.matrix, y=NULL, index = "jack1")
j1 # 43.969697
#89%
j2 <- diversityresult(com.matrix, y=NULL, index = "jack2")
j2 # 44.98167
#87%

#jar jackknife; richness = 26
j1.j <- diversityresult(jar.com.matrix, y=NULL, index = "jack1")
j1.j # 32.833333
#79%
j2.j <- diversityresult(jar.com.matrix, y=NULL, index = "jack2")
j2.j # 35.783391
#73%

#pitfall jackknife; richness = 21
j1.p <- diversityresult(pitfall.com.matrix, y=NULL, index = "jack1")
j1.p # 24.897436
#84%
j2.p <- diversityresult(pitfall.com.matrix, y=NULL, index = "jack2")
j2.p # 25.921053
#81%

#ramp jackknife; richness = 35
j1.r <- diversityresult(ramp.com.matrix, y=NULL, index = "jack1")
j1.r # 41.833333
#84%
j2.r <- diversityresult(ramp.com.matrix, y=NULL, index = "jack2")
j2.r # 46.641696
#75%

#sticky jackknife; richness = 31
j1.s <- diversityresult(sticky.com.matrix, y=NULL, index = "jack1")
j1.s # 36.857143
#84%
j2.s <- diversityresult(sticky.com.matrix, y=NULL, index = "jack2")
j2.s # 39.783972
#78%

#BiodiversityR::accumcomp
Accum.1_functional <- accumcomp(com.matrix, y=env.matrix, factor='Trap', 
                     method='random', conditioned=FALSE, plotit=FALSE)
Accum.1_functional

#BiodiversityR::accumcomp.long
accum.long1_functional <- accumcomp.long(Accum.1_functional, ci=NA, label.freq=5)
head(accum.long1_functional)

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

functional_accum <- ggplot(data=accum.long1_functional, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_color_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  scale_shape_manual(values=c(19,17,15,25))+
  geom_line(aes(colour=Grouping), size=0.1) +
  geom_ribbon(aes(colour=Grouping, fill=after_scale(alpha(colour, 0.3))), 
              show.legend=FALSE, linetype = 0) + 
  geom_point(data=subset(accum.long1_functional, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=3) +
  BioR.theme +
  labs(x = "", y = "", colour = "Trap", shape = "Trap")
functional_accum

pdf("functional_accum.pdf", height=6, width=8) #height and width in inches
functional_accum
dev.off()

########################################################################
#beetles

#bring in beetle data sets from github

pitfall_beetle <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/2020%20beetles_pitfall.csv",na.strings = NULL)

ramp_beetle <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/2020%20beetles_yellowramp.csv",na.strings = NULL)

jar_beetle <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/2020%20beetles_jarramp.csv",na.strings = NULL)

sticky_beetle <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/2020%20beetles_stickycard.csv",na.strings = NULL)

taxa_beetle <- read.csv("https://raw.githubusercontent.com/katiemmanning/trapping_calibration/main/Data/beetle%20taxa.csv")

#add trap type as a column on each data file
pitfall_beetle$Trap="pitfall"
ramp_beetle$Trap="ramp"
jar_beetle$Trap="jar"
sticky_beetle$Trap="sticky"

#combine beetle data tables 
library (plyr)
pitfallramp_beetle <- rbind.fill (pitfall_beetle, ramp_beetle)
pitfallrampjar_beetle <-rbind.fill (pitfallramp_beetle, jar_beetle)
beetle <- rbind.fill (pitfallrampjar_beetle, sticky_beetle)

#############
#NMDS of beetle community between trap types
library (vegan)

#Create matrix of environmental variables
env.matrix_beetle<-beetle[c(1:3,19)]
#create matrix of community variables
com.matrix_beetle<-beetle[c(4:18)]

#ordination by NMDS
NMDS_beetle<-metaMDS(com.matrix_beetle, distance="bray", k=2, autotransform=FALSE, trymax=100)
stressplot(NMDS_beetle)
#stress= "nearly zero"
#INSUFFICENT DATA

#beetle NMDS visualization 

#what taxa to display using "taxa"
#most.abund<-as.vector(t(taxa[1,]))#greater than 100 caught
#most.abund<-most.abund[-1]
#bioind<-as.vector(t(taxa[2,]))
#bioind<-bioind[-1]
flying<-as.vector(t(taxa_beetle[1,]))
flying<-flying[-1]
crawling<-as.vector(t(taxa_beetle[2,]))
crawling<-crawling[-1]
include<-as.vector(t(taxa_beetle[3,]))
include<-include[-1]

#plot beetle NMDS
plot(NMDS_beetle, disp='sites', type="n")
#add ellipsoids with ordiellipse
ordiellipse(NMDS_beetle, env.matrix_beetle$Trap, draw="polygon", col="#F0E442",kind="sd", conf=0.95, label=FALSE, show.groups = "ramp") 
ordiellipse(NMDS_beetle, env.matrix_beetle$Trap, draw="polygon", col="#CC79A7",kind="sd", conf=0.95, label=FALSE, show.groups = "sticky")
ordiellipse(NMDS_beetle, env.matrix_beetle$Trap, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "pitfall")
ordiellipse(NMDS_beetle, env.matrix_beetle$Trap, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "jar") 
#display ground trap data as solid shapes - pitfall=circle, ramp trap=square, jar=triangle, flying trap as triangle outline
points(NMDS_beetle, display="sites", select=which(env.matrix_beetle$Trap=="pitfall"),pch=19, col="#E69F00")
points(NMDS_beetle, display="sites", select=which(env.matrix_beetle$Trap=="jar"), pch=17, col="#009E73")
points(NMDS_beetle, display="sites", select=which(env.matrix_beetle$Trap=="ramp"), pch=15, col="#F0E442")
points(NMDS_beetle, display="sites", select=which(env.matrix_beetle$Trap=="sticky"), pch=25, col="#CC79A7")
#add legend
#legend(1.5,1.5, title=NULL, pch=c(19,17,15,25), col=c("#E69F00","#009E73","#F0E442","#CC79A7"), cex=.7, legend=c("Pitfall", "Jar ramp", "Yellow ramp", "Yellow sticky card"))
#add taxa as text
ordilabel(NMDS_beetle, display="species", select =which (include==TRUE & crawling == TRUE), cex=0.6, col="black", fill="white")
ordilabel(NMDS_beetle, display="species", select =which (include==TRUE & flying == TRUE), cex=0.6, col="white", fill="black")

#bootstrapping and testing for differences between the groups (traps)
fit<-adonis(com.matrix_beetle ~ Trap, data = env.matrix_beetle, permutations = 999, method="bray")
fit
#P-value = 0.003

#check assumption of homogeneity of multivariate dispersion 
#P-value greater than 0.05 means assumption has been met
distances_data<-vegdist(com.matrix_beetle)
anova(betadisper(distances_data, env.matrix_beetle$Trap))
#P-value = 0.7987 --- assumes homogeneity


################
#calculate beetle Abundance
insects.abun_beetle <- rowSums(beetle[,4:18])
beetle$abundance <- insects.abun_beetle

#calculate beetle Richness
insects.rowsums_beetle <- rowSums(beetle[,4:18]>0)
beetle$richness <- insects.rowsums_beetle

#calculate beetle Shannon diversity
diversity_beetle <-diversity(beetle[,4:18])
beetle$diversity <-diversity_beetle

#calculate beetle Evenness
evenness_beetle <-diversity_beetle/log(specnumber(beetle[,4:18]))
beetle$evenness <- evenness_beetle

#######
#Mixed effects models
library(lme4)
library(lmerTest) #to obtain p values
library (emmeans) #for pairwise comparisons
library (multcompView) #to view letters

#beetle richness
##AIC 77
richness.model_beetle<-lmer(richness ~ Trap + Date + (1 | Site), data=beetle)
summary(richness.model_beetle)
anova(richness.model_beetle)
AIC(richness.model_beetle)
#pairwise comparison 
rich.emm_beetle<-emmeans(richness.model_beetle,pairwise~Trap)
rich.emm_beetle
#results: 
rich.cld_beetle<-multcomp::cld(rich.emm_beetle, alpha = 0.05, Letters = LETTERS)
rich.cld_beetle

#beetle abundance
##AIC 77
abundance.model_beetle<-lmer(abundance ~ Trap + Date + (1 | Site), data=beetle)
summary(abundance.model_beetle)
anova(abundance.model_beetle)
AIC(abundance.model_beetle)
#pairwise comparison 
abun.emm_beetle<-emmeans(abundance.model_beetle,pairwise~Trap)
abun.emm_beetle
#results: 
abun.cld_beetle<-multcomp::cld(abun.emm_beetle, alpha = 0.05, Letters = LETTERS)
abun.cld_beetle

#beetle diversity
##AIC 54 (40 w/o date)
diversity.model_beetle<-lmer(diversity ~ Trap + Date + (1 | Site), data=beetle)
summary(diversity.model_beetle)
anova(diversity.model_beetle)
AIC(diversity.model_beetle)
#pairwise comparison 
div.emm_beetle<-emmeans(diversity.model_beetle,pairwise~Trap)
div.emm_beetle
#results: 
div.cld_beetle<-multcomp::cld(div.emm_beetle, alpha = 0.05, Letters = LETTERS)
div.cld_beetle

#beetle evenness
##AIC -193 (-411 w/o date)
evenness.model_beetle<-lmer(evenness ~ Trap + Date + (1 | Site), data=beetle)
summary(evenness.model_beetle)
anova(evenness.model_beetle)
AIC(evenness.model_beetle)
#pairwise comparison 
even.emm_beetle<-emmeans(evenness.model_beetle,pairwise~Trap)
even.emm_beetle
#results: 
even.cld_beetle<-multcomp::cld(even.emm_beetle, alpha = 0.05, Letters = LETTERS)
even.cld_beetle

###########
library(ggplot2)
#beetle abundance plot
abundance.plot_beetle<-ggplot(beetle, aes(x =Trap, y = abundance, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="")+
  scale_y_continuous(trans="log10")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=abun.cld_beetle, aes(y = 600, label = .group))
abundance.plot_beetle

#beetle richness plot
richness.plot_beetle<-ggplot(beetle, aes(x =Trap, y = richness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=rich.cld_beetle, aes(y = 25, label = .group))
richness.plot_beetle

#beetle diversity plot
diversity.plot_beetle<-ggplot(beetle, aes(x =Trap, y = diversity, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=div.cld_beetle, aes(y = 2.5, label = .group))
diversity.plot_beetle

#beetle evenness plot
evenness.plot_beetle<-ggplot(beetle, aes(x =Trap, y = evenness, fill=Trap))+
  geom_boxplot()+
  theme_bw()+
  theme(legend.position ="NULL")+
  theme(axis.text.x=element_blank())+
  labs(x="", y="")+
  scale_fill_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  geom_text(data=even.cld_beetle, aes(y = 1.2, label = .group))
evenness.plot_beetle

#Mush order plots together
library(ggpubr) 
beetlefigure <- ggarrange(richness.plot_beetle, abundance.plot_beetle, diversity.plot_beetle, evenness.plot_beetle,
                         labels = c("A", "B", "C", "D"),
                         ncol = 2, nrow = 2,
                         common.legend = TRUE, legend = "bottom")
pdf("beetle.pdf", height=6, width=8) #height and width in inches
beetlefigure
dev.off()

beetlefigure

####
#species accumulation
library (BiodiversityR)
library(ggplot2)

#individual curves for each trap type
pitfall.com.matrix<-pitfall_beetle[c(4:18)]
pitfall_curve<-accumresult(pitfall.com.matrix, method = "exact", permutations = 1000)

jar.com.matrix<-jar_beetle[c(4:18)]
jar_curve<-accumresult(jar.com.matrix, method = "exact", permutations = 1000)

ramp.com.matrix<-ramp_beetle[c(4:18)]
ramp_curve<-accumresult(ramp.com.matrix, method = "exact", permutations = 1000)

sticky.com.matrix<-sticky_beetle[c(4:18)]
sticky_curve<-accumresult(sticky.com.matrix, method = "exact", permutations = 1000)

#first-order jackknife estimates are based on the number of singletons
#second-order jackknife estimates are based on the number of singletons and doubletons

#calculates species richness for each sample
specnumber(com.matrix_beetle) #ranges from 1 to 3

#calculates species richness by treatment (trap)
specnumber(com.matrix_beetle, groups = beetle$Trap) #jar=3; pitfall=3; ramp=8; sticky=9

#total richness and jackknife
rich <- diversityresult(com.matrix_beetle, y=NULL, index = "richness")
rich # 15
j1 <- diversityresult(com.matrix_beetle, y=NULL, index = "jack1")
j1 # 18.878788
#79%
j2 <- diversityresult(com.matrix_beetle, y=NULL, index = "jack2")
j2 # 18.996212
#79%

#jar jackknife; richness = 3
j1.j <- diversityresult(jar.com.matrix, y=NULL, index = "jack1")
j1.j # 3.8333333
#78%
j2.j <- diversityresult(jar.com.matrix, y=NULL, index = "jack2")
j2.j #3.9666667
#76%

#pitfall jackknife; richness = 3
j1.p <- diversityresult(pitfall.com.matrix, y=NULL, index = "jack1")
j1.p # 3.8333333
#78%
j2.p <- diversityresult(pitfall.com.matrix, y=NULL, index = "jack2")
j2.p # 3.9666667
#76%

#ramp jackknife; richness = 8
j1.r <- diversityresult(ramp.com.matrix, y=NULL, index = "jack1")
j1.r # 11.555556
#69%
j2.r <- diversityresult(ramp.com.matrix, y=NULL, index = "jack2")
j2.r # 13.305556
#60%

#sticky jackknife; richness = 9
j1.s <- diversityresult(sticky.com.matrix, y=NULL, index = "jack1")
j1.s # 11.75
#77%
j2.s <- diversityresult(sticky.com.matrix, y=NULL, index = "jack2")
j2.s # 11.219697
#80%

#BiodiversityR::accumcomp
Accum.1_beetle <- accumcomp(com.matrix_beetle, y=env.matrix_beetle, factor='Trap', 
                     method='random', conditioned=FALSE, plotit=FALSE)
Accum.1_beetle

#BiodiversityR::accumcomp.long
accum.long1_beetle <- accumcomp.long(Accum.1_beetle, ci=NA, label.freq=5)
head(accum.long1_beetle)

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

beetle_accum <- ggplot(data=accum.long1_beetle, aes(x = Sites, y = Richness, ymax = UPR, ymin = LWR)) + 
  scale_x_continuous(expand=c(0, 1), sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_color_manual(values=c("#009E73","#E69F00","#F0E442","#CC79A7"))+
  scale_shape_manual(values=c(19,17,15,25))+
  geom_line(aes(colour=Grouping), size=0.1) +
  geom_ribbon(aes(colour=Grouping, fill=after_scale(alpha(colour, 0.3))), 
              show.legend=FALSE, linetype = 0) + 
  geom_point(data=subset(accum.long1_beetle, labelit==TRUE), 
             aes(colour=Grouping, shape=Grouping), size=3) +
  BioR.theme +
  labs(x = "Number of samples", y = "", colour = "Trap", shape = "Trap")
beetle_accum

pdf("beetle_accum.pdf", height=6, width=8) #height and width in inches
beetle_accum
dev.off()

#######
#manuscript figures

#Figure 1 - trap photos

#Figure 2 - NMDSs

#a - order
plot(NMDS_order, disp='sites', type="n")
#add ellipsoids with ordiellipse
ordiellipse(NMDS_order, env.matrix_order$Trap, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "pitfall")
ordiellipse(NMDS_order, env.matrix_order$Trap, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "jar") 
ordiellipse(NMDS_order, env.matrix_order$Trap, draw="polygon", col="#F0E442",kind="sd", conf=0.95, label=FALSE, show.groups = "ramp") 
ordiellipse(NMDS_order, env.matrix_order$Trap, draw="polygon", col="#CC79A7",kind="sd", conf=0.95, label=FALSE, show.groups = "sticky")
#display ground trap data as solid shapes - pitfall=circle, ramp trap=square, jar=triangle, flying trap as triangle outline
points(NMDS_order, display="sites", select=which(env.matrix_order$Trap=="pitfall"),pch=19, col="#E69F00")
points(NMDS_order, display="sites", select=which(env.matrix_order$Trap=="jar"), pch=17, col="#009E73")
points(NMDS_order, display="sites", select=which(env.matrix_order$Trap=="ramp"), pch=15, col="#F0E442")
points(NMDS_order, display="sites", select=which(env.matrix_order$Trap=="sticky"), pch=25, col="#CC79A7")
#add legend
legend(1.20,1.35, title=NULL, pch=c(19,17,15,25), col=c("#E69F00","#009E73","#F0E442","#CC79A7"), cex=.7, legend=c("Pitfall", "Jar ramp", "Yellow ramp", "Yellow sticky card"))
#add insect taxa as text
ordilabel(NMDS_order, display="species", select =which (include==TRUE & crawling == TRUE), cex=0.6, col="black", fill="white")
ordilabel(NMDS_order, display="species", select =which (include==TRUE & flying == TRUE), cex=0.6, col="white", fill="black")

#b - functional
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
#legend(1.46,1.45, title=NULL, pch=c(19,17,15,25), col=c("#E69F00","#009E73","#F0E442","#CC79A7"), cex=.7, legend=c("Pitfall", "Jar ramp", "Yellow ramp", "Yellow sticky card"))
#add insect taxa as text
ordilabel(NMDS, display="species", select =which (include==TRUE & crawling == TRUE), cex=0.6, col="black", fill="white")
ordilabel(NMDS, display="species", select =which (include==TRUE & flying == TRUE), cex=0.6, col="white", fill="black")

#c - beetles
plot(NMDS_beetle, disp='sites', type="n")
#add ellipsoids with ordiellipse
ordiellipse(NMDS_beetle, env.matrix_beetle$Trap, draw="polygon", col="#E69F00",kind="sd", conf=0.95, label=FALSE, show.groups = "pitfall")
ordiellipse(NMDS_beetle, env.matrix_beetle$Trap, draw="polygon", col="#009E73",kind="sd", conf=0.95, label=FALSE, show.groups = "jar") 
ordiellipse(NMDS_beetle, env.matrix_beetle$Trap, draw="polygon", col="#F0E442",kind="sd", conf=0.95, label=FALSE, show.groups = "ramp") 
ordiellipse(NMDS_beetle, env.matrix_beetle$Trap, draw="polygon", col="#CC79A7",kind="sd", conf=0.95, label=FALSE, show.groups = "sticky")
#display ground trap data as solid shapes - pitfall=circle, ramp trap=square, jar=triangle, flying trap as triangle outline
points(NMDS_beetle, display="sites", select=which(env.matrix_beetle$Trap=="pitfall"),pch=19, col="#E69F00")
points(NMDS_beetle, display="sites", select=which(env.matrix_beetle$Trap=="jar"), pch=17, col="#009E73")
points(NMDS_beetle, display="sites", select=which(env.matrix_beetle$Trap=="ramp"), pch=15, col="#F0E442")
points(NMDS_beetle, display="sites", select=which(env.matrix_beetle$Trap=="sticky"), pch=25, col="#CC79A7")
#add legend
#legend(5,5, title=NULL, pch=c(19,17,15,25), col=c("#E69F00","#009E73","#F0E442","#CC79A7"), cex=.7, legend=c("Pitfall", "Jar ramp", "Yellow ramp", "Yellow sticky card"))
#add taxa as text
ordilabel(NMDS_beetle, display="species", select =which (include==TRUE & crawling == TRUE), cex=0.6, col="black", fill="white")
ordilabel(NMDS_beetle, display="species", select =which (include==TRUE & flying == TRUE), cex=0.6, col="white", fill="black")

#Figure 3 - trap comparison box plots

#a - order
library(ggpubr) 
orderfigure <- ggarrange(richness.plot_order, abundance.plot_order, diversity.plot_order, evenness.plot_order,
                         ncol = 4, nrow = 1)
orderfigure

#b - functional 
functionalfigure <- ggarrange(richness.plot, abundance.plot, diversity.plot, evenness.plot,
                              ncol = 4, nrow = 1)
functionalfigure

#c - beetle
beetlefigure <- ggarrange(richness.plot_beetle, abundance.plot_beetle, diversity.plot_beetle, evenness.plot_beetle,
                          ncol = 4, nrow = 1,
                          common.legend = TRUE, legend = "bottom")
beetlefigure

figure3 <- ggarrange(orderfigure, functionalfigure, beetlefigure,
                          labels = c("A", "B", "C"),
                          ncol = 1, nrow = 3,
                          common.legend = TRUE, legend = "bottom")
pdf("Figure 3.pdf", height=10, width=15) #height and width in inches
figure3
dev.off()
figure3

#Figure 4 - flying vs crawling (functional level)
##insert code here once data is updated

#Figure 5 - accumulation plots

#a - order
#b - functional 
#c - beetles

figure5 <- ggarrange(order_accum, functional_accum, beetle_accum,
                     labels = c("A", "B", "C"),
                     ncol = 1, nrow = 3,
                     common.legend = TRUE, legend = "bottom")
figure5

pdf("Figure 5.pdf", height=6, width=6) #height and width in inches
figure5
dev.off()


#Supplementary figure 1 - trap size vs mean catch
#currently on excel
