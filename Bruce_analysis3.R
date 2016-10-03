#This file analyzes Bruce Hill data on powerlines.

# In this study BH surveyed 10 study sites (2km2). 
# 5 study sites were bisected by a maintained powerline corridor, 5 were not.
# The habitat types surveyed in this study were:
# a)    Corridors;
# b)	Forest;
# c)	Forest/grassland boundary;
# d)	Non-flowering crop boundary;
# e)	Semi-natural grassland;
# f)	Maintained roadside;
# g)	Maintained drain.
# Each study site contained examples of between 4-6 habitat types.   
# Within each habitat types the most flowery 50m section was surveyed.  
# The survey was done along a 50m long by 3m wide transect in which bumblebee abundance, 
# by species, was recorded for 15 minutes.  Collection handling time was discounted and 
# when the end of a transect was reached before 15 minutes the transect was resurveyed.  
# If the bee was feeding the flower species was recorded, otherwise the bee was recorded 
# as “flying”.  Flower density, being the total percentage of the transect area covered 
# by flowers was also estimated and recorded. 
# Analysis of the combined results of all of transects in one study site is referred to 
# as the “site scale”. Analysis of each transect is referred to as the “plot scale”
# Each study plot was surveyed twice between 9th July 2014 and 25th August 2014.  
# The first round of surveying took place between 9th July- 27th July, the second 
# round took place between 29th July-25th August. 

#This file will answer the following questions:
# 1. Did study plots containing a power line corridor (corridor) 
# have more bees abundance and/or more bee species richness than study plots not containing 
# a corridor (non-corridor)? 

# 2. Did the more regularly maintained/disturbed habitat types (being corridors,
# maintained road margins, maintained drains and semi-natural grasslands) in both 
# corridor and non-corridor sites have more bees/species than the less regularly 
# maintained/disturbed habitat types (being forest, forest/grassland boundary and 
# non flowering crop edge)?

# 3. Explore the beta-diversity component

# 4. What are the influences of
# flower abundance on both bee abundance and the number of bee species in these habitats?
# Which plants sustain more bee species?


#load libraries----
library(reshape)
library(vegan)
library(bipartite)
library(nlme)


#then we read and understand the data----

d <- read.csv("data/powerlines3.csv", h = TRUE)
head(d)
str(d)

#And check for mistakes.

levels(d$Gen_sp) 
levels(d$Flower_density)
#I rename them for analysis
levels(d$Flower_density) <- c("5","20", "40", "60", "10", "80", "1")
table(d$Flower_density) 
d$Flower_density <- as.numeric(as.character(d$Flower_density))

levels(d$Flower_species) 
levels(d$Site)
#issues with swedish letters, I rename them
levels(d$Site)[1] <- "Angeby"
levels(d$Site)[4] <- "Gavsta"
levels(d$Site)[7] <- "Kottgrind"
levels(d$Site)[8] <- "Laby_Norra"
levels(d$Site)[9] <- "Laby_Sodra"

#check comments and check first site number of transects
unique(d$Notes)
d[which(d$Notes %in% unique(d$Notes)[c(2,8,9,19)]),]
#i understand all is ok.

#check for some basic things (e.g. weather):
boxplot(d$Highest_temp ~ d$Corridor)
boxplot(d$Highest_temp ~ d$Habitat) #is forest cooler b/c of shadow?
    #Nice, no visual differences!
boxplot(d$Highest_temp ~ d$Site) 
boxplot(d$Highest_temp ~ d$Round)
    #expected, second round was overall cooler.

table(d$Habitat, d$Site)
table(d$Habitat, d$Corridor)

# 1) richness and abundance at the Site level.----
# main concern is using same sampling effort, so things are comparable.
# for that we have to calculate abundance and richness at the transect level

str(d)
d2 <- cast(d, Site + Plot + Corridor + Round + Flower_density + Habitat ~ Gen_sp, length, 
           value = "Date")
head(d2)
#NA in species are sites with 0 individuals, so can be discarded.

#calculate abundance and richness
abundance <- rowSums(d2[,7:26])
d2bin <- ifelse(d2[,c(7:21,23:26)] > 0, 1, 0)
head(d2bin)
richness <- rowSums(d2bin)
#calculate chao2/ACE estimator here also.
chao2 <- estimateR(d2[,c(7:21,23:26)])[2,]
ACE <- estimateR(d2[,c(7:21,23:26)])[4,]
#plot(jitter(richness), chao2, xlim = c(0,20), ylim = c(0,20))
#plot(richness, ACE, xlim = c(0,30), ylim = c(0,30))

#cretae new dataset
dat <- data.frame(d2[,1:6], abundance, richness, chao2, ACE)
head(dat)

#some plotting for answers 1 (site) and 2 (habitat) and also flower abundance
tapply(dat$Flower_density, dat$Habitat, mean)
tapply(dat$Flower_density, dat$Habitat, summary)
boxplot(dat$abundance ~ dat$Corridor, las = 1) #same abundnce per transect in both sites...
boxplot(dat$abundance ~ dat$Habitat, las = 2) #not big diferences among habitats...

#FIG1#----
#barplots for abundance
#I use barplots for consistency with betadiversity results
ab <- tapply(dat$abundance, dat$Corridor, mean, na.rm = TRUE)
ab_n <- tapply(dat$abundance, dat$Corridor, length)
ab_sd <- tapply(dat$abundance, dat$Corridor, sd, na.rm = TRUE)/sqrt(ab_n)

names(ab)[1] <- "No corridor"
names(ab)[2] <- "Corridor"
barplot(ab, las = 1, ylim = c(0,8), ylab = "Abundance")
abline(h=0)
arrows(x0 = 0.7, y0 = ab[1]+ab_sd[1], x1 = 0.7, y1 = ab[1]-ab_sd[1], 
       angle = 90, code = 3)
arrows(x0 = 1.9, y0 = ab[2]+ab_sd[2], x1 = 1.9, y1 = ab[2]-ab_sd[2], 
       angle = 90, code = 3)


ab <- tapply(dat$abundance, dat$Habitat, mean, na.rm = TRUE)
ab_n <- tapply(dat$abundance, dat$Habitat, length)
ab_sd <- tapply(dat$abundance, dat$Habitat, sd, na.rm = TRUE)/sqrt(ab_n)

names(ab)[3] <- "Boundary"
names(ab)[4] <- "Drain"
names(ab)[5] <- "Roadside"
names(ab)[6] <- "Crop edge"
names(ab)[7] <- "Grassland"

sort_id <- order(ab)
ab <- ab[sort_id]
ab_sd <- ab_sd[sort_id]

barplot(ab, las = 2, ylim = c(0,13), ylab = "Abundance")
abline(h=0)
arrows(x0 = 0.7, y0 = ab[1]+ab_sd[1], x1 = 0.7, y1 = ab[1]-ab_sd[1], 
       angle = 90, code = 3)
arrows(x0 = 1.9, y0 = ab[2]+ab_sd[2], x1 = 1.9, y1 = ab[2]-ab_sd[2], 
       angle = 90, code = 3)
arrows(x0 = 3.1, y0 = ab[3]+ab_sd[3], x1 = 3.1, y1 = ab[3]-ab_sd[3], 
       angle = 90, code = 3)
arrows(x0 = 4.3, y0 = ab[4]+ab_sd[4], x1 = 4.3, y1 = ab[4]-ab_sd[4], 
       angle = 90, code = 3)
arrows(x0 = 5.5, y0 = ab[5]+ab_sd[5], x1 = 5.5, y1 = ab[5]-ab_sd[5], 
       angle = 90, code = 3)
arrows(x0 = 6.7, y0 = ab[6]+ab_sd[6], x1 = 6.7, y1 = ab[6]-ab_sd[6], 
       angle = 90, code = 3)
arrows(x0 = 7.9, y0 = ab[7]+ab_sd[7], x1 = 7.9, y1 = ab[7]-ab_sd[7], 
       angle = 90, code = 3)
####

boxplot(dat$richness ~ dat$Corridor, las = 1) #same pattern...
boxplot(dat$richness ~ dat$Habitat, las = 1) #not big diferences among habitats...

boxplot(dat$chao2 ~ dat$Corridor, las = 1) #same pattern...
boxplot(dat$chao2 ~ dat$Habitat, las = 1) #not big diferences among habitats...

#boxplot(dat$ACE ~ dat$Corridor, las = 1) #same pattern...
#boxplot(dat$ACE ~ dat$Habitat, las = 1) #not big diferences among habitats...

#answers 1, 2 and 3 (site and habitat and flower density) can be answered toguether.

#MODEL PRESENTED:
m <- lme(abundance ~ Flower_density + Habitat + Corridor, 
         na.action = na.omit,
         random = ~1|Site/Plot, data = dat,  weights = varIdent(form=~1|Habitat))
plot(m) #residuals are more or less ok.
summary(m)
anova(m)

#Calculate power
library(simr)
#run omly in lme4
m <- lmer(abundance ~ Flower_density + Habitat + Corridor + (1 | Site/Plot), 
         na.action = na.omit,
         data = dat) # lmer can't handle Varident #, weights = varIdent(form=~1|Habitat))

fixef(m)["HabitatSemi_natural_grasslands"] <- 1.69 #add an effect size expected of 25% difference between habitat and powerline.

#Calculate power
#powerSim(m, test = fixed("Flower_density")) 
powerSim(m, test = fixed("Habitat")) 
ps_hab <- lastResult()

#Calculate power

m <- lmer(abundance ~ Flower_density + Habitat + Corridor + (1 | Site/Plot), 
          na.action = na.omit,
          data = dat) # lmer can't handle Varident #, weights = varIdent(form=~1|Habitat))
fixef(m)["CorridorYes"] <- 1.36 #add an effect size expected.
powerSim(m, test = fixed("Corridor")) 
ps_cor <- lastResult()

#next
plot(dat$abundance ~ dat$Flower_density, las = 1)
abline(lm(dat$abundance ~ dat$Flower_density))

m <- lme(richness ~ Flower_density + Habitat + Corridor, na.action = na.omit,
         random = ~1|Site/Plot, data = dat, weights = varIdent(form=~1|Habitat))
plot(m)
summary(m)
anova(m)

#power
#run omly in lme4
m <- lmer(richness ~ Flower_density + Habitat + Corridor + (1 | Site/Plot), 
          na.action = na.omit,
          data = dat) # lmer can't handle Varident #, weights = varIdent(form=~1|Habitat))
tapply(dat$richness, dat$Habitat, mean)
fixef(m)["HabitatSemi_natural_grasslands"] <- 0.53 

#Calculate power
#powerSim(m, test = fixed("Flower_density")) 
powerSim(m, test = fixed("Habitat")) 
ps_hab <- lastResult()
m <- lmer(richness ~ Flower_density + Habitat + Corridor + (1 | Site/Plot), 
          na.action = na.omit,
          data = dat) # lmer can't handle Varident #, weights = varIdent(form=~1|Habitat))
tapply(dat$richness, dat$Corridor, mean)
fixef(m)["CorridorYes"] <- 0.48 
powerSim(m, test = fixed("Corridor")) 
ps_cor <- lastResult()


plot(dat$richness ~ dat$Flower_density, las = 1)
abline(lm(dat$richness ~ dat$Flower_density))

###Chao 2 same results
m <- lme(chao2 ~ Flower_density + Habitat + Corridor, na.action = na.omit,
         random = ~1|Site/Plot, data = dat,  weights = varIdent(form=~1|Habitat))
plot(m)
summary(m)
anova(m)


# Beta diversity analysis----

#following tylianakis et al 2005

#first: comparision at landscape scale
#H0 = beta among transects in sites with corridor should be lower than among transects in sites without
#if corridors are used as dispersal. 

d3 <- cast(d, Site + Corridor + Plot + Habitat ~ Gen_sp, length, 
           value = "Id")
head(d3)

#Calculate alpha as mean species per plot in each site
d3bin <- ifelse(d3[,c(5:19,21:24)] > 0, 1, 0)
head(d3bin)
richness <- rowSums(d3bin)

datB <- data.frame(d3[,1:2], richness)
alpha <- mean(datB$richness)
alpha_site <- tapply(datB$richness, d3$Plot, FUN = mean)
alpha_site
alpha_site_se <- tapply(datB$richness, d3$Plot, FUN = function(x) sd(x)/sqrt(length(x)))
alpha_site_se

#Calculate beta as mean of sp per site minus alpha
d_site <- cast(d, Site + Corridor ~ Gen_sp, length, 
               value = "Id")
head(d_site)

dsitebin <- ifelse(d_site[,c(3:17,19:22)] > 0, 1, 0)
head(dsitebin)
richness <- rowSums(dsitebin)

dat_site <- data.frame(d_site[,1:2], richness)
betas <- dat_site$richness - alpha_site 

#comparision corridor and non corridor sites
alpha_corr <- tapply(alpha_site, d_site$Corridor, FUN = mean)
alpha_corr
beta_corr <- tapply(betas, d_site$Corridor, FUN = mean)
beta_corr
alpha_corr_se <- tapply(alpha_site, d_site$Corridor, FUN = function(x) sd(x)/sqrt(length(x)))
beta_corr_se <- tapply(betas, d_site$Corridor, FUN = function(x) sd(x)/sqrt(length(x)))

gamma <- beta_corr+alpha_corr
gamma

#Fig 2A------
barplot(t(as.matrix(data.frame(alpha_corr, beta_corr))), ylim = c(0,15), las = 1,
        ylab = "Richness", names.arg = c("No Corridor", "Corridor"))
abline(h = 0)
arrows(x0 = 0.7 ,x1 = 0.7 ,y0 = alpha_corr[1] - alpha_corr_se[1], y1 = alpha_corr[1] + alpha_corr_se[1],
       angle = 90, code = 3)
arrows(x0 = 1.9 ,x1 = 1.9 ,y0 = alpha_corr[2] - alpha_corr_se[2], y1 = alpha_corr[2] + alpha_corr_se[2],
       angle = 90, code = 3)

arrows(x0 = 0.7, x1 = 0.7 ,y0 = gamma[1] - beta_corr_se[1], y1 = gamma[1] + beta_corr_se[1],
       angle = 90, code = 3)
arrows(x0 = 1.9 ,x1 = 1.9 ,y0 = gamma[2] - beta_corr_se[2], y1 = gamma[2] + beta_corr_se[2],
       angle = 90, code = 3)
##

#Sedond: use beta across habitats

#make a matrix pooling rounds.
d3 <- cast(d, Site + Plot + Corridor + Habitat ~ Gen_sp, length, 
           value = "Id")
head(d3)

#following tylianakis et al 2005
#alpha as mean species per plot
d3bin <- ifelse(d3[,c(5:19,21:24)] > 0, 1, 0)
head(d3bin)
richness <- rowSums(d3bin)

datB <- data.frame(d3[,1:4], richness)
alpha <- mean(datB$richness)
d3$site_hab <- paste(d3$Site,d3$Habitat, sep = "_")
alpha_hab <- tapply(datB$richness,d3$site_hab , FUN = mean)
alpha_hab

#beta as mean of sp per site minus alpha
#BUT need to correct for number of transects per site!!

#make a matrix pooling rounds.
d3 <- cast(d, Site + Plot + Corridor + Habitat ~ Gen_sp, length, 
           value = "Id")
head(d3)

#following tylianakis et al 2005
#alpha as mean species per plot
d3bin <- ifelse(d3[,c(5:19,21:24)] > 0, 1, 0)
head(d3bin)
richness <- rowSums(d3bin)

datB <- data.frame(d3[,1:4], richness)
alpha_ <- mean(datB$richness)
alpha <- tapply(datB$richness, d3$Habitat, FUN = mean)
alpha
alpha_se <- tapply(datB$richness, d3$Habitat, FUN = function(x) sd(x)/sqrt(length(x)))

#beta as mean of sp per habitat RAREFIED! to account for different sampling effort minus alpha
d_site <- cast(d, Habitat ~ Gen_sp, length, 
               value = "Id")
head(d_site)

rowSums(d_site[,c(2:17,18:21)]) #90
richness90 <- rarefy(d_site[,c(2:17,18:21)], 90)
richness70 <- rarefy(d_site[,c(2:17,18:21)], 70)
plot(richness90~richness70)
richness <- richness90

dat_site <- data.frame(Habitat = d_site[,1], richness = as.vector(richness))
betas <- dat_site$richness - alpha 

gamma <- betas + alpha
gamma

diversity <- data.frame(alpha, betas, alpha_se)
plot(diversity$alpha ~ diversity$betas)

#fig 2B----
#sort
sort_id <- order(diversity$alpha)
diversity <- diversity[sort_id,]
rownames(diversity)
barplot(t(as.matrix(diversity[,1:2])), ylim = c(0,16),
        names.arg = c("Forest", "Crop edge", "Boundary", "Drain", "Corridor", "Roadside",
                      "Grassland"), las = 2, ylab = "Richness")
abline(h = 0)
#for alpha
arrows(x0 = 0.7 ,x1 = 0.7 ,y0 = diversity$alpha[1] - diversity$alpha_se[1], y1 = diversity$alpha[1] + diversity$alpha_se[1],
       angle = 90, code = 3)
arrows(x0 = 1.9 ,x1 = 1.9 ,y0 = diversity$alpha[2] - diversity$alpha_se[2], y1 = diversity$alpha[2] + diversity$alpha_se[2],
       angle = 90, code = 3)
arrows(x0 = 3.1 ,x1 = 3.1 ,y0 = diversity$alpha[3] - diversity$alpha_se[3], y1 = diversity$alpha[3] + diversity$alpha_se[3],
       angle = 90, code = 3)
arrows(x0 = 4.3 ,x1 = 4.3 ,y0 = diversity$alpha[4] - diversity$alpha_se[4], y1 = diversity$alpha[4] + diversity$alpha_se[4],
       angle = 90, code = 3)
arrows(x0 = 5.5 ,x1 = 5.5 ,y0 = diversity$alpha[5] - diversity$alpha_se[5], y1 = diversity$alpha[5] + diversity$alpha_se[5],
       angle = 90, code = 3)
arrows(x0 = 6.7 ,x1 = 6.7 ,y0 = diversity$alpha[6] - diversity$alpha_se[6], y1 = diversity$alpha[6] + diversity$alpha_se[6],
       angle = 90, code = 3)
arrows(x0 = 7.9 ,x1 = 7.9 ,y0 = diversity$alpha[7] - diversity$alpha_se[7], y1 = diversity$alpha[7] + diversity$alpha_se[7],
       angle = 90, code = 3)
###


#analyze conservation and Ecosystem Species Provider (ESP) species separatelly----

ESP <- c("Bombus terrestis", "Bombus lapidarius", "Bombus pascuorum", "Bombus hypnorum", "B.pratorum", "B. hortorum")
    
CV <- c("Bombus muscorum", "Bombus humilis", "Bombus sorensis", "Bombus Sylvarum")

head(d)
d4 <- cast(d, Plot + Habitat + Flower_density + Site ~ Gen_sp, length, 
           value = "Date")
head(d4)

d4$abundance_ESP <- d4$Bombus_terrestris + d4$Bombus_lapidarius + d4$Bombus_pascuorum + d4$Bombus_hypnorum + 
    d4$Bombus_pratorum + d4$Bombus_hortorum
d4$abundance_CV <- d4$Bombus_muscorum + d4$Bombus_humilis + d4$Bombus_soroeensis + d4$Bombus_sylvarum

#test model for abundance of key groups.
#ESP
head(d4)
ab <- tapply(d4$abundance_ESP, d4$Habitat, mean, na.rm = TRUE)
ab_n <- tapply(d4$abundance_ESP, d4$Habitat, length)
ab_sd <- tapply(d4$abundance_ESP, d4$Habitat, sd, na.rm = TRUE)/sqrt(ab_n)

names(ab)[3] <- "Boundary"
names(ab)[4] <- "Drain"
names(ab)[5] <- "Roadside"
names(ab)[6] <- "Crop edge"
names(ab)[7] <- "Grassland"

sort_id <- order(ab)
ab <- ab[sort_id]
ab_sd <- ab_sd[sort_id]
barplot(ab, las = 2, ylim = c(0,7), ylab = "Abundance of Ecosystem services providers")
abline(h=0)
arrows(x0 = 0.7, y0 = ab[1]+ab_sd[1], x1 = 0.7, y1 = ab[1]-ab_sd[1], 
       angle = 90, code = 3)
arrows(x0 = 1.9, y0 = ab[2]+ab_sd[2], x1 = 1.9, y1 = ab[2]-ab_sd[2], 
       angle = 90, code = 3)
arrows(x0 = 3.1, y0 = ab[3]+ab_sd[3], x1 = 3.1, y1 = ab[3]-ab_sd[3], 
       angle = 90, code = 3)
arrows(x0 = 4.3, y0 = ab[4]+ab_sd[4], x1 = 4.3, y1 = ab[4]-ab_sd[4], 
       angle = 90, code = 3)
arrows(x0 = 5.5, y0 = ab[5]+ab_sd[5], x1 = 5.5, y1 = ab[5]-ab_sd[5], 
       angle = 90, code = 3)
arrows(x0 = 6.7, y0 = ab[6]+ab_sd[6], x1 = 6.7, y1 = ab[6]-ab_sd[6], 
       angle = 90, code = 3)
arrows(x0 = 7.9, y0 = ab[7]+ab_sd[7], x1 = 7.9, y1 = ab[7]-ab_sd[7], 
       angle = 90, code = 3)

#MODEL SHOWN:
m <- lme(abundance_ESP ~ Flower_density + Habitat, 
         na.action = na.omit,
         random = ~1|Site, data = d4, weights = varIdent(form=~1|Habitat))
plot(m) #residuals are more or less ok.
summary(m)
anova(m)

#CV
head(d4)
ab <- tapply(d4$abundance_CV, d4$Habitat, mean, na.rm = TRUE)
ab_n <- tapply(d4$abundance_CV, d4$Habitat, length)
ab_sd <- tapply(d4$abundance_CV, d4$Habitat, sd, na.rm = TRUE)/sqrt(ab_n)

names(ab)[3] <- "Boundary"
names(ab)[4] <- "Drain"
names(ab)[5] <- "Roadside"
names(ab)[6] <- "Crop edge"
names(ab)[7] <- "Grassland"

sort_id <- order(ab)
ab
ab <- ab[sort_id]
ab_sd <- ab_sd[sort_id]
barplot(ab, las = 2, ylim = c(0,1.5), ylab = "Abundance of conservation value species")
abline(h=0)
arrows(x0 = 0.7, y0 = ab[1]+ab_sd[1], x1 = 0.7, y1 = ab[1]-ab_sd[1], 
       angle = 90, code = 3)
arrows(x0 = 1.9, y0 = ab[2]+ab_sd[2], x1 = 1.9, y1 = ab[2]-ab_sd[2], 
       angle = 90, code = 3)
arrows(x0 = 3.1, y0 = ab[3]+ab_sd[3], x1 = 3.1, y1 = ab[3]-ab_sd[3], 
       angle = 90, code = 3)
arrows(x0 = 4.3, y0 = ab[4]+ab_sd[4], x1 = 4.3, y1 = ab[4]-ab_sd[4], 
       angle = 90, code = 3)
arrows(x0 = 5.5, y0 = ab[5]+ab_sd[5], x1 = 5.5, y1 = ab[5]-ab_sd[5], 
       angle = 90, code = 3)
arrows(x0 = 6.7, y0 = ab[6]+ab_sd[6], x1 = 6.7, y1 = ab[6]-ab_sd[6], 
       angle = 90, code = 3)
arrows(x0 = 7.9, y0 = ab[7]+ab_sd[7], x1 = 7.9, y1 = ab[7]-ab_sd[7], 
       angle = 90, code = 3)

#MODEL SHOWN:
m <- lme(abundance_CV ~ Flower_density + Habitat, 
         na.action = na.omit,
         random = ~1|Site/Plot, data = d4, weights = varIdent(form=~1|Habitat))
plot(m) #residuals are more or less ok.
summary(m)
anova(m)

#4) flower species strngth----

#for that I pool everything and create a network of plants per bees. 

str(d)
levels(d$Flower_species)
dd <- subset(d, Flower_species != "Flying")
ntw <- cast(dd, Flower_species ~ Gen_sp, fun = length)
head(ntw)
ntw <- ntw[,c(-17,-22)] 
rownames(ntw) <- ntw[,1]
ntw <- ntw[,-1]
ntw <- as.data.frame(ntw)
#Fig5 (get better names to save space)
colnames(ntw) <- gsub("Bombus_", "B.", colnames(ntw))
rownames(ntw) <- c("Arctium tomentosum"          ,"Calluna vulgaris"                            
, "Campanulaceae rapunculoides"                  ,"Campanulaceae rotundifolia"
, "Carduus arvense"                              ,"Carduus crispus"                             
, "Carduus helenioides"                          ,"Centaurea cyanus"                            
, "Centaurea jacea"                              ,"Centaurea scabiosa"                          
, "Cirsium arvense"                              ,"Crepis tectorum"               
, "Epilobium adenocaulon"                        ,"Filipendula ulmaria"                         
, "Galeopsis terrahit"                           ,"Hypericum maculatum"          
, "Lamiastrum galeobdolon"                       ,"Lamium maculatum"                            
, "Lathyrus pratensis"                           ,"Leontodon autumnalis"        
, "Lythranceae salcaria"                         ,"Malva spp"                                  
, "Melampyrum pratense"                          ,"Prunella vulgaris"                           
, "Satureja vulgaris"                            ,"Solidago virgaurea"           
, "Sonchus glabrescens"                      ,"Succisa pratensis"                           
, "Taraxacum spp"                                ,"Trifolium hybridum"             
, "Trifolium medium"                             ,"Trifolium pratense"                          
, "Trifolium repens"    ,"Vicia cracca") 
par(mar = c(13,4,4,2), xpd = TRUE)
plotweb(ntw, labsize = 0.9, text.rot = 90, low.y = 0.6)
rownames(ntw)
colnames(ntw)

#calculate plant strength
S <- specieslevel(ntw, index = "species strength", level = "lower")
S$sp <- rownames(S)
S[order(S$species.strength),]
#we can add plant abundance/numbler of links, sum of link abundance
write.table(x = S, "Strength.txt")

#subset per corridors
head(dd)
dd2 <- subset(dd, Habitat == "Corridor")
ntw <- cast(dd2, Flower_species ~ Gen_sp, fun = length)
head(ntw)
ntw <- ntw[,c(-15,-20)] 
rownames(ntw) <- ntw[,1]
ntw <- ntw[,-1]
ntw <- as.data.frame(ntw)
#calculate plant strength
S <- specieslevel(ntw, index = "species strength", level = "lower")
S$sp <- rownames(S)
S[order(S$species.strength),]

#subset per grassland
head(dd)
dd2 <- subset(dd, Habitat == "Semi_natural_grasslands")
ntw <- cast(dd2, Flower_species ~ Gen_sp, fun = length)
head(ntw)
ntw <- ntw[,c(-11)] 
rownames(ntw) <- ntw[,1]
ntw <- ntw[,-1]
ntw <- as.data.frame(ntw)
#calculate plant strength
S <- specieslevel(ntw, index = "species strength", level = "lower")
S$sp <- rownames(S)
S[order(S$species.strength),]



#create some descrpitive data (not directly used)----
dd <- cast(d,  Habitat ~ Gen_sp, length, 
           value = "Date")

barplot(as.numeric(dd[1,-1]))

ddd <- t(dd)
colnames(ddd) <- c("Corridor" ,"Forest","Boundary", "Maintained_drain",         
                   "Roadside", "Crop_edge", "Grasslands")  
ddd  

#list species per habitat? 
colnames(d)
list <- unique(d[,c("Habitat", "Gen_sp")])
list[order(list$Habitat),]
table(list$Habitat)
unique(d$Gen_sp) #wow, all species are found in corridors!!


