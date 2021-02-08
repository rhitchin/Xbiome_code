setwd("~/Desktop/ClinStat/BoneDensity")
bones <- read.csv("SUMMER_2017_FINALEXAM_DATASET.CSV", header = T)
bones <- as.data.frame(bones)
bonez <- bones

library(ggplot2)
library(e1071)
library(plyr)
#library(RStata)
#library(nortest)
#library(MASS)
#library(Hmisc)
#library(Deducer)
library(corrplot)
library(pspearman)
#library(BSDA)
#library(stats)
#library(questionr)
#library(compareGroups)
#library(pwr)
#library(rms)
library(epitools)

### Check everything for nomality, just in case
ggplot(bonez, aes(x = dexa)) + 
  geom_histogram(aes(y = ..density..), fill = "light blue", bins = 50) + 
  stat_function( 
    fun = dnorm, 
    args = with(bonez, c(mean = mean(dexa), sd = sd(dexa))) 
  ) + 
  scale_x_continuous("Dexa") +
  ggtitle("Dexa Scores")

ggplot(bonez, aes(x = vitd)) + 
  geom_histogram(aes(y = ..density..), fill = "light blue", bins = 50) + 
  stat_function( 
    fun = dnorm, 
    args = with(bonez, c(mean = mean(vitd), sd = sd(vitd))) 
  ) + 
  scale_x_continuous("Vitd") +
  ggtitle("Vitamin D")

ggplot(bonez, aes(x = ses)) + 
  geom_histogram(aes(y = ..density..), fill = "light blue", bins = 25) + 
  stat_function( 
    fun = dnorm, 
    args = with(bonez, c(mean = mean(ses), sd = sd(ses))) 
  ) + 
  scale_x_continuous("SES") +
  ggtitle("Socio-Economic Status")

ggplot(bonez, aes(x = age)) + 
  geom_histogram(aes(y = ..density..), fill = "light blue", bins = 7) + 
  stat_function( 
    fun = dnorm, 
    args = with(bonez, c(mean = mean(age), sd = sd(age))) 
  ) + 
  scale_x_continuous("years") +
  ggtitle("Age")

ggplot(bonez, aes(x = bmi)) + 
  geom_histogram(aes(y = ..density..), fill = "light blue", bins = 25) +
  stat_function( 
    fun = dnorm, 
    args = with(bonez, c(mean = mean(bmi), sd = sd(bmi))) 
  ) + 
  scale_x_continuous("BMI") +
  ggtitle("Body-Mass Index")

ggplot(bonez, aes(x = calcium)) + 
  geom_histogram(aes(y = ..density..), fill = "light blue", bins = 20) + 
  stat_function( 
    fun = dnorm, 
    args = with(bonez, c(mean = mean(calcium), sd = sd(calcium))) 
  ) + 
  scale_x_continuous("Ca") +
  ggtitle("Blood Calcium Levels")

#nhanes1$bmicat <- ifelse(nhanes1$bmi <= 29.999,
                         #ifelse(nhanes1$bmi <= 25, "Normal", "Overweight"), "Obese")

bonez$vitdcat <- ifelse(bonez$vitd < 20.0, "Low", "Normal")
bonez$vitdcat <- as.factor(bonez$vitdcat)
lowvit <- subset(bonez, vitd < 20)
normvit <- subset(bonez, vitd >= 20)
blacklow <- subset(bonez, race =="black" & vitd < 20)
blackhi <- subset(bonez, race =="black" & vitd >= 20)
realrace <- subset(bonez, race =="black" | race =="hispanic")
hisplow <- subset(bonez, race == "hispanic" & vitd < 20)
hisphi <- subset(bonez, race == "hispanic" & vitd >= 20)
table(realrace$race, realrace$vitdcat)
table(bonez$multivitamin, bonez$vitdcat)
prop.test(c(17,53),c(39,83))
table(bonez$multivitamin, bonez$vitdcat)
prop.test(c(28,28+33),c(32,32+50))
table(bonez$season, bonez$vitdcat)
prop.test(c(1,1+17),c(59,59+66))
#do all interquartile ranges with IQR()
IQR(bonez$ses)
IQR(bonez$soda)
IQR(lowvit$soda)
IQR(normvit$soda)
IQR(bonez$age)
IQR(lowvit$age)
IQR(normvit$age)

t.test(lowvit$ses, normvit$ses, var.equal=F, paired=F)
t.test(lowvit$soda, normvit$soda, var.equal=F, paired=F)
t.test(lowvit$vitd, normvit$vitd, var.equal=T, paired=F)
t.test(lowvit$age, normvit$age, var.equal=F, paired=F)
t.test(lowvit$dexa, normvit$dexa, var.equal=T, paired=F)
t.test(lowvit$bmi, normvit$bmi, var.equal=T, paired=F)
t.test(lowvit$phys_act, normvit$phys_act, var.equal=T, paired=F)
t.test(lowvit$calcium, normvit$calcium, var.equal=T, paired=F)
t.test(lowvit$dexa, normvit$dexa, var.equal=T, paired=F)
t.test(lowvit$dexa, normvit$dexa, var.equal=T, paired=F)


### Let's build some early single-variable regression models to look for correlations
simp <- lm(vitd ~ dexa, data = bonez)
plot(bones$dexa, bones$vitd, pch = 1, cex = 0.5, col = "green", main = "Simple Linear Regression", xlab = "DEXA", ylab = "Vitamin D")
abline(simp)
summary(simp)
## Statistically significant correlation, B-value of 1.87, Rsquare 0.14

calc <- lm(calcium ~ dexa, data = bonez)
plot(bones$dexa, bones$calcium, pch = 1, cex = 0.5, col = "green", main = "Simple Linear Regression", xlab = "DEXA", ylab = "Calcium")
abline(calc)
summary(calc)
## absolutely no correlation, B-value of 0.01, Rsq 0.001
simp <- lm(vitd ~ dexa, data = bonez)
plot(bones$dexa, bones$vitd, pch = 1, cex = 0.5, col = "green", main = "Simple Linear Regression", xlab = "DEXA", ylab = "Vitamin D")
abline(simp)
summary(simp)
###moderate, significant correlation

###Test for interaction between BMI and vitamin D
simp2 <- lm(vitd ~ dexa + bmi, data = bonez)
summary(simp2)

#visually check for interactions using plot(), coplot(), and termplot()
plot(simp2)
termplot(simp2)
coplot(vitd ~ dexa|bmi, data = bonez, panel = panel.smooth)
### No obvious interaction, lets check the residuals of the multiplication term
m.simp <- lm(vitd ~ dexa:bmi, data = bonez)
plot(m.simp)
summary(m.simp)
### Again no substantial evidence of interaction, let's build a model containing a few of the other variables and see if adding interaction terms changes B significantly
multi.plot <- lm(vitd ~ dexa + bmi + sex + age, data = bonez)
plot(multi.plot)
termplot(multi.plot)
summary()
summary(multi.plot)
###Plot regression with 95% CI
ggplot(bonez, aes(x=dexa, y=vitd))+
  geom_point()+
  geom_smooth(method=lm, se=TRUE)


#summary(lm(dexa ~ vitd, data = bonez))
plot(lm(dexa ~ vitd + dexa:bmi, data = bonez))
summary(lm(dexa ~ vitd, data = bonez))
summary(lm(dexa ~ vitd + season + age, data = bonez))
summary(lm(dexa ~ vitd + bmi, data = bonez))
summary(lm(dexa ~ vitd + bmi + dexa:bmi, data = bonez))
summary(lm(dexa ~ vitd + dexa:bmi, data = bonez))


#Dichotomized LM's
skinny <- subset(bonez, bmi < 85)
fat <- subset(bonez, bmi >= 85)
summary(lm(dexa ~ vitd, data = skinny))
summary(lm(dexa ~ vitd, data = fat))

#Odds ratio to appease my statistical overlords
bonez$dexcat <- ifelse(bonez$dexa < -2.0, "Low", "Normal")
bonez$dexcat <- as.factor(bonez$dexcat)
table(bonez$dexcat, bonez$vitdcat)
table(skinny$dexcat, skinny$vitdcat)
table(fat$dexcat, fat$vitdcat)


### Wow season seems to have an effect, let's graph!
plot(bones$season, bones$dexa, pch = 1, cex = 0.5, col = "brown", main = "Dexa by Season", xlab = "Season", ylab = "Dexa")
warm <- subset(bonez, season == "Summer/Spring")
cold <- subset(bonez, season == "Winter/Fall")
#warmline <- lm(vitd ~ dexa, data = warm)
#coldline <- lm(vitd ~ dexa, data = cold)
#plot(bonez$vitd, bonez$dexa, pch = 1, cex = 0.5, col = bonez$season, main = "By Season", ylab = "DEXA", xlab = "Vitamin D")
#abline(warmline)
###Trying to add both regression lines to the same plot, don't know if this is worth it.
ggplot(bonez, aes(x=dexa, y=vitd, color = season))+
  geom_point()+
  geom_smooth(method='lm',formula=y~x)+
  ggtitle("Vitamin D and Bone Density")+
  xlab("DEXA Z-Score")+
  ylab("Vitamin D (ng/mL)")
t.test(warm$dexa, cold$dexa, paired = F, var.equal = T)
t.test(warm$vitd, cold$vitd, paired = F, var.equal = T)
t.test(bonez$dexa, mu = 0)

ggplot(bonez, aes(x=dexa, y=vitd, color = race))+
  geom_point()+
  geom_smooth(method='lm',formula=y~x)+
  ggtitle("Vitamin D and Bone Density")+
  xlab("DEXA Z-Score")+
  ylab("Vitamin D (ng/mL)")
t.test(warm$dexa, cold$dexa, paired = F, var.equal = T)
t.test(warm$vitd, cold$vitd, paired = F, var.equal = T)
t.test(bonez$dexa, mu = 0)

#ggplot(bonez, aes(x=vitd, y=dexa, color = season))+
  #geom_point()+
  #geom_smooth(method='lm',formula=y~x)

### Cool! definitely Season plays a role in changing total Vitamin D in children. It seems like kids who show up to the clinic in summer/spring are just a different population to kids in the winter/fall.

### Calculate Power for this study
### Let's calculate for an effect size of one SD in DEXA change, because that correlates to an 80% (meaningful) change in fracture chance
library(pwr)
#pwr.t.test(d = #ONE SD(1), sig.level = 0.05, power = 0.8)
## Using Rsquare value from model of dexa ~ vitd + season
pwr.r.test(power = NULL, n = 143, r = 0.376)
#power = 0.9968629


###BMI Calculations

bonez$bmi_dich <- ifelse(bonez$bmi >= 85, "Over Mean BMI", "Under Mean BMI")
ggplot(bonez, aes(x=dexa, y=vitd, color = bmi_dich))+
  geom_point()+
  geom_smooth(method='lm',formula=y~x)+
  ggtitle("Vitamin D and Bone Density")+
  xlab("DEXA Z-Score")+
  ylab("Vitamin D (ng/mL)")
t.test(warm$dexa, cold$dexa, paired = F, var.equal = T)
t.test(warm$vitd, cold$vitd, paired = F, var.equal = T)
t.test(bonez$dexa, mu = 0)

fat <- subset(bonez, bmi >= 85)
skinny <- subset(bonez, bmi < 85)
ggplot(fat, aes(x=dexa, y=vitd))+
  geom_point()+
  geom_smooth(method='lm',formula=y~x)+
  ggtitle("Vitamin D and Bone Density")+
  xlab("DEXA Z-Score")+
  ylab("Vitamin D (ng/mL)")
t.test(warm$dexa, cold$dexa, paired = F, var.equal = T)
t.test(warm$vitd, cold$vitd, paired = F, var.equal = T)
t.test(bonez$dexa, mu = 0)

ggplot(skinny, aes(x=dexa, y=vitd))+
  geom_point()+
  geom_smooth(method='lm',formula=y~x)+
  ggtitle("Vitamin D and Bone Density")+
  xlab("DEXA Z-Score")+
  ylab("Vitamin D (ng/mL)")

t.test(warm$dexa, cold$dexa, paired = F, var.equal = T)
t.test(warm$vitd, cold$vitd, paired = F, var.equal = T)
t.test(bonez$dexa, mu = 0)

###Correlation Matrix
corrsub <- bonez[c("ses", "vitd", "age", "dexa", "bmi", "phys_act", "calcium")]
colnames(corrsub) <- c("SES", "Vitamin D", "Age", "DEXA", "BMI", "Physical Activity", "Calcium")
bonez.spearman <- cor(corrsub, method = "spearman", use = "complete.obs")
bonez.pearson <- cor(corrsub, method = "pearson", use = "complete.obs")

corrplot.mixed(bonez.spearman, title = "Spearman Correlation")

#png(filename="~/Desktop/CRTP/Exam/Pearson_correlation_matrix")
corrplot.mixed(bonez.pearson)
#dev.off()


cor.mtest <- function(bonez.pearson, conf.level = 0.95){
  mat <- as.matrix(bonez.pearson)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

res1 <- cor.mtest(corrsub,0.95)
res2 <- cor.mtest(corrsub,0.99)
## specialized the insignificant value according to the significant level
corrplot(bonez.pearson, p.mat = res1[[1]], sig.level=0.05)
corrplot(bonez.pearson, p.mat = res1[[1]], insig = "p-value")

png(filename="~/Desktop/CRTP/Exam/Spearman_correlation_FIXED")
corrplot(bonez.spearman, p.mat = res1[[1]], insig = "p-value")
dev.off()

cor(x=bonez$dexa, y=bonez$vitd, method = "spearman", use = "complete.obs")

