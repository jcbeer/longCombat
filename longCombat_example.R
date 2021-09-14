###########################################################
# longCombat package examples
# JCBeer joanne.beer@pennmedicine.upenn.edu
# 13 Sept 2021
###########################################################

#################################
# install longCombat package
#################################
# install.packages('devtools')
# devtools::install_github("jcbeer/longCombat")

#################################
# load longCombat package
#################################
library(longCombat)
# check documentation
?longCombat

#################################
# install and load invgamma & lmer package
#################################
# install.packages('invgamma')
library(invgamma)
library(lme4)

#################################
# simulate data to run the functions
#################################
# 100 subjects with 5 time points each
# 8 scanners / batches
# 20 features
#################################
# set random seed
set.seed(1)
# simulate the covariates
simdata <- data.frame(
  subid=rep(1:100, each=5),
  age=rep(sample(c(20:60), 100, replace=TRUE), each=5),
  diagnosis=rep(c(0,1), each=250),
  time=rep(0:4, times=100)
  )
# define 4 batch patterns 
# each column of this matrix represents a batch pattern over time
batch.patterns <- matrix(c(1,1,1,2,2,
                           3,3,4,4,5,
                           6,6,6,6,6,
                           7,7,8,8,8),
                         ncol=4)
# randomly sample 100 batch patterns
batch.pattern.sample <- sample(1:4, 100, replace=TRUE)
simdata$batch <- as.vector(batch.patterns[,batch.pattern.sample])
# simulate the brain features 
features <- matrix(rnorm(100*5*20), nrow=500)
# simulate additive batch effects (normally distributed)
gamma <- runif(n=8, min=-5, max=5)
tau <- runif(n=8, min=0.1, max=0.3)
batch.add <- matrix(c(
  rnorm(mean=gamma[1], sd=tau[1], n=20),
  rnorm(mean=gamma[2], sd=tau[2], n=20),
  rnorm(mean=gamma[3], sd=tau[3], n=20),
  rnorm(mean=gamma[4], sd=tau[4], n=20),
  rnorm(mean=gamma[5], sd=tau[5], n=20),
  rnorm(mean=gamma[6], sd=tau[6], n=20),
  rnorm(mean=gamma[7], sd=tau[7], n=20),
  rnorm(mean=gamma[8], sd=tau[8], n=20)),
  ncol=8) 
# simulate multiplicative batch effects (inverse gamma distributed)
lambda <- sample(c(2, 3), 8, replace=TRUE)
theta <- sample(c(0.5, 1), 8, replace=TRUE)
batch.mult <- matrix(c(
  rinvgamma(n=20, shape=lambda[1], scale=theta[1]),
  rinvgamma(n=20, shape=lambda[2], scale=theta[2]),
  rinvgamma(n=20, shape=lambda[3], scale=theta[3]),
  rinvgamma(n=20, shape=lambda[4], scale=theta[4]),
  rinvgamma(n=20, shape=lambda[5], scale=theta[5]),
  rinvgamma(n=20, shape=lambda[6], scale=theta[6]),
  rinvgamma(n=20, shape=lambda[7], scale=theta[7]),
  rinvgamma(n=20, shape=lambda[8], scale=theta[8])),
  ncol=8)
# add / multiply batch effects to the features
for(i in 1:500){
  features[i,] <- features[i,]*batch.mult[,simdata$batch[i]] + batch.add[,simdata$batch[i]]
}
# add covariate effects to the features
features <- features - 0.1*simdata$age + simdata$diagnosis - 0.5*simdata$time - 2*simdata$diagnosis*simdata$time 
# add subject random effect to the features (will be the same across features in this case)
features <- features + rep(rnorm(n=100), each=5)
# save feature names
featurenames <- paste0('feature', 1:20)
colnames(features) <- featurenames
# combine into one data frame
simdata <- data.frame(simdata, features)
# remove some stuff no longer needed
rm('features', 'batch.patterns', 'batch.pattern.sample', 'i')

#################################
# longCombat functions:
#################################
# batchTimeViz() -- visualize change in batch over time
# batchBoxplot() -- to visualize residuals across batches
# trajPlot() -- visualize trajectories
# addTest() -- test for additive scanner effects
# multTest() -- test for multiplicative scanner effects
# longCombat() -- apply longitudinal ComBat
#################################
# type e.g. ?longCombat to get further documentation
#################################

# these examples (ranef='(1|subid)') are for random subject intercept
# use ranef='(1 + time|subid)' to add random slope

#################################
# batchTimeViz() -- visualize change in batch over time
#################################
batchTimeViz(batchvar='batch',
             timevar='time',
             data=simdata)

#################################
# batchBoxplot() -- to visualize residuals across batches
# can do for each feature you are interested in
#################################
# make batch boxplot for feature1, do not adjust for batch 
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature1', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=1:8)

# make batch boxplot for feature2, do not adjust for batch 
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature2', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',   
             data=simdata,
             colors=1:8)

# make batch boxplot for feature2, DO adjust for batch 
# order by increasing batch variance
# (centers boxplot means on the zero line)
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature2', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             adjustBatch=TRUE,
             orderby='var',
             colors=1:8)

#################################
# trajPlot() -- visualize trajectories
#################################
# for everyone
trajPlot(idvar='subid', 
         timevar='time',
         feature='feature2', 
         batchvar='batch',  
         data=simdata,
         point.col=simdata$batch,
         line.col=simdata$diagnosis[!duplicated(simdata$subid)]+1)
# for only diagnosis=0
trajPlot(idvar='subid', 
         timevar='time',
         feature='feature2', 
         batchvar='batch',  
         data=simdata[simdata$diagnosis==0,],
         point.col=simdata$batch[simdata$diagnosis==0])
# for only diagnosis=1
trajPlot(idvar='subid', 
         timevar='time',
         feature='feature2', 
         batchvar='batch',  
         data=simdata[simdata$diagnosis==1,],
         point.col=simdata$batch[simdata$diagnosis==1],
         line.col=rep(2,250))

#################################
# addTest() -- test for additive scanner effects
#################################
addTestTable <- addTest(idvar='subid', 
        batchvar='batch', 
        features=featurenames, 
        formula='age + diagnosis*time',
        ranef='(1|subid)',
        data=simdata)

# when we generate data with set.seed(1)
# feature5 has largest additive scanner effects
# (since it is the first row in the addTestTable)
# check boxplot to see this
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature5', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=1:8,
             title='Feature 5')

# compare with feature 2
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature2', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=1:8,
             title='Feature 2')

#################################
# multTest() -- test for multiplicative scanner effects
#################################
multTestTable <- multTest(idvar='subid', 
                        batchvar='batch', 
                        features=featurenames, 
                        formula='age + diagnosis*time',
                        ranef='(1|subid)',
                        data=simdata)

# when we generate data with set.seed(1)
# feature3 has largest multiplicative scanner effects
# (since it is the first row in the multTestTable)
# check boxplot to see this
# (we will adjust for batch and order by variance
# to best see the multiplicative batch effects)
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature3', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=1:8,
             adjustBatch=TRUE,
             orderby='var',
             title='Feature 3')

# compare with feature1
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature1', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=1:8,
             adjustBatch=TRUE,
             orderby='var',
             title='Feature 1')

#################################
# longCombat() -- apply longitudinal ComBat
#################################
simdata_combat <- longCombat(idvar='subid', 
                             timevar='time',
                             batchvar='batch', 
                             features=featurenames, 
                             formula='age + diagnosis*time',
                             ranef='(1|subid)',
                             data=simdata)

#################################
# get the harmonized data
simdata_harmonized <- simdata_combat$data_combat
# save combat feature names
featurenames.combat <- names(simdata_harmonized)[4:23]
# merge with original dataframe
simdata <- merge(simdata, simdata_harmonized[,c(1,2,4:23)], by=c('subid', 'time'))

#################################
# test for additive scanner effects in combatted data
#################################
addTestTableCombat <- addTest(idvar='subid', 
                              batchvar='batch', 
                              features=featurenames.combat, 
                              formula='age + diagnosis*time',
                              ranef='(1|subid)',
                              data=simdata)

# there are still some significant additive batch effects (p<0.05)
# but p-values tend to be larger (-log10(p-values are smaller)) in overall distribution
boxplot(-log(as.numeric(addTestTable$`KR p-value`), base=10),
        -log(as.numeric(addTestTableCombat$`KR p-value`), base=10),
        ylim=c(0, 8),
        las=1,
        ylab='additive batch effect -log10(p-value)',
        names=c('before ComBat', 'after ComBat'))

# check feature 5 boxplot before combat
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature5', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=1:8,
             title='feature 5 before combat')

# check feature 5 boxplot after combat
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature5.combat', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=1:8,
             title='feature 5 after combat')

#################################
# test for multiplicative scanner effects in combatted data
#################################
multTestTableCombat <- multTest(idvar='subid', 
                                batchvar='batch', 
                                features=featurenames.combat, 
                                formula='age + diagnosis*time',
                                ranef='(1|subid)',
                                data=simdata)

# there are still some significant multiplicative batch effects (p<0.05)
# but p-values tend to be larger (-log10(p-values are smaller)) in overall distribution
boxplot(-log(as.numeric(multTestTable$`p-value`), base=10),
        -log(as.numeric(multTestTableCombat$`p-value`), base=10),
        las=1,
        ylab='multiplicative batch effect -log10(p-value)',
        names=c('before ComBat', 'after ComBat'))

# check feature3 boxplot before combat
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature3', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=1:8,
             adjustBatch=TRUE,
             orderby='var',
             title='Feature 3 before ComBat')

# check feature3 boxplot after combat
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature3.combat', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=1:8,
             adjustBatch=TRUE,
             orderby='var',
             title='Feature 3 after ComBat')

#################################
# plot trajectories before and after combat
#################################
par(mfrow=c(1,2))
trajPlot(idvar='subid', 
         timevar='time',
         feature='feature2', 
         batchvar='batch',  
         data=simdata,
         ylimits=c(-30, 12),
         title='feature 2 before combat',
         point.col=simdata$batch,
         line.col=simdata$diagnosis[!duplicated(simdata$subid)]+1)

trajPlot(idvar='subid', 
         timevar='time',
         feature='feature2.combat', 
         batchvar='batch',  
         data=simdata,
         ylimits=c(-30,12),
         title='feature 2 after combat',
         point.col=simdata$batch,
         line.col=simdata$diagnosis[!duplicated(simdata$subid)]+1)

#################################
# fit LME models before / after ComBat
#################################
# before ComBat
feature5.fit <- lmer(feature5 ~ age + diagnosis*time + (1|subid), data=simdata)
feature5.batch.fit <- lmer(feature5 ~ age + diagnosis*time + as.factor(batch) + (1|subid), data=simdata)
# after ComBat
feature5combat.fit <- lmer(feature5.combat ~ age + diagnosis*time + (1|subid), data=simdata)
feature5combat.batch.fit <- lmer(feature5.combat ~ age + diagnosis*time + as.factor(batch) + (1|subid), data=simdata)

summary(feature5.fit)
summary(feature5.batch.fit)

summary(feature5combat.fit)
summary(feature5combat.batch.fit)
