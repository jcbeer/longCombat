###########################################################
# longCombat package examples
# JCBeer joanne.beer@pennmedicine.upenn.edu
# 11 Feb 2021
###########################################################

#################################
# install longCombat package
#################################
devtools::install_github("jcbeer/longCombat")

#################################
# load longCombat package
#################################
library(longCombat)

#################################
# simulate data to run the functions
#################################
# 10 subjects with 5 time points each
# 3 scanners / batches
# 20 features
#################################
# set random seed
set.seed(1)
# simulate the covariates
simdata <- data.frame(
  subid=rep(1:10, each=5),
  age=rep(sample(c(20:60), 10), each=5),
  diagnosis=rep(c(0,1), each=25),
  time=rep(0:4, times=10),
  batch=sample(1:3, 50, replace=TRUE)
  )
# simulate the brain features 
features <- matrix(rnorm(50*20), nrow=50)
# add a batch effect to the features
features <- (features*simdata$batch) + 2*simdata$batch
# add covariate effects to the features
features <- features + simdata$subid -0.01*simdata$age + simdata$diagnosis - 0.5*simdata$time - simdata$diagnosis*simdata$time 
# save feature names
featurenames <- paste0('feature', 1:20)
colnames(features) <- featurenames
# combine into one data frame
simdata <- data.frame(simdata, features)
rm('features')

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
# (NOTE: for the simulated data
# each batch is distributed over all time points
# so the plot is not that interesting)
#################################
batchTimeViz(batchvar='batch',
             timevar='time',
             data=simdata)

#################################
# batchBoxplot() -- to visualize residuals across batches
# can do for each feature you are interested in
#################################
# color by batch
cols <- c('lightsalmon', 'lightblue', 'darkseagreen1')

# make batch boxplot for feature1, do not adjust for batch 
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature1', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=cols)

# make batch boxplot for feature2, do not adjust for batch 
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature2', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',   
             data=simdata,
             colors=cols)

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
             colors=cols)

#################################
# trajPlot() -- visualize trajectories
#################################
# for everyone
trajPlot(idvar='subid', 
         timevar='time',
         feature='feature2', 
         batchvar='batch',  
         data=simdata)
# for only diagnosis=0
trajPlot(idvar='subid', 
         timevar='time',
         feature='feature2', 
         batchvar='batch',  
         data=simdata[simdata$diagnosis==0,])
# for only diagnosis=1
trajPlot(idvar='subid', 
         timevar='time',
         feature='feature2', 
         batchvar='batch',  
         data=simdata[simdata$diagnosis==1,])

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
# feature6 has largest additive scanner effects
# (since it is the first row in the addTestTable)
# check boxplot to see this
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature6', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=cols)

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
# feature8 has largest multiplicative scanner effects
# (since it is the first row in the multTestTable)
# check boxplot to see this
# (we will adjust for batch and order by variance
# to best see the multiplicative batch effects)
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature8', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=cols,
             adjustBatch=TRUE,
             orderby='var')

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
# but greatly reduced

# check feature 6 boxplot before combat
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature6', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=cols,
             title='feature 6 before combat')

# check feature 6 boxplot after combat
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature6.combat', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=cols,
             title='feature6 after combat')

# check feature 19 boxplot before combat
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature19', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=cols,
             title='feature 19 before combat')

# check feature 19 boxplot after combat
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature19.combat', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=cols,
             title='feature19 after combat')

#################################
# test for multiplicative scanner effects in combatted data
#################################
multTestTableCombat <- multTest(idvar='subid', 
                                batchvar='batch', 
                                features=featurenames.combat, 
                                formula='age + diagnosis*time',
                                ranef='(1|subid)',
                                data=simdata)

# there is still one significant multiplicative batch effect (p<0.05)

# check feature8 boxplot before combat
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature8', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=cols,
             adjustBatch=TRUE,
             orderby='var')

# check feature8 boxplot after combat
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature8.combat', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=cols,
             adjustBatch=TRUE,
             orderby='var')

# check feature2 boxplot before combat
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature2', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=cols,
             adjustBatch=TRUE,
             orderby='var')

# check feature2 boxplot after combat
batchBoxplot(idvar='subid', 
             batchvar='batch', 
             feature='feature2.combat', 
             formula='age + diagnosis*time',
             ranef='(1|subid)',
             data=simdata,
             colors=cols,
             adjustBatch=TRUE,
             orderby='var')

#################################
# plot trajectories before and after combat
#################################
par(mfrow=c(1,2))
trajPlot(idvar='subid', 
         timevar='time',
         feature='feature2', 
         batchvar='batch',  
         data=simdata,
         ylimits=c(0,20),
         title='feature 2 before combat')

trajPlot(idvar='subid', 
         timevar='time',
         feature='feature2.combat', 
         batchvar='batch',  
         data=simdata,
         ylimits=c(0,20),
         title='feature 2 after combat')
