#Scripts to runABCRF from the simulations run with ABCtoolbox + fastsimcoal

#!/usr/bin/Rscript

setwd("")
library(abcrf)

# chose one 3gp_apricot.obs
statobs <- read.table("3gp_apricot.obs",header = TRUE)
# read all 30 simulations' results

### local computer ###

# 30m without gene flow
ss_SET1_30msc1 <- read.table("30m_res_sc1.txt", header=TRUE, dec=".",stringsAsFactors=F,sep="\t", colClasses="numeric" )
ss_SET1_30msc2 <- read.table("30m_res_sc2.txt", header=TRUE, dec=".",stringsAsFactors=F,sep="\t", colClasses="numeric" )
# with gene flow
ss_SET1_30msc1_05 <- read.table("30m_res_sc1_05.txt", header=TRUE, dec=".",stringsAsFactors=F,sep="\t", colClasses="numeric" )
ss_SET1_30msc2_05 <- read.table("30m_res_sc2_05.txt", header=TRUE, dec=".",stringsAsFactors=F,sep="\t", colClasses="numeric" )

ss_SET1_30msc1_005 <- read.table("30m_res_sc1_005.txt", header=TRUE, dec=".",stringsAsFactors=F,sep="\t", colClasses="numeric" )
ss_SET1_30msc2_005 <- read.table("30m_res_sc2_005.txt", header=TRUE, dec=".",stringsAsFactors=F,sep="\t", colClasses="numeric" )

#delete column that are not ss

ss_SET1_30msc1 <- ss_SET1_30msc1[1:7500,-(1:6)]
ss_SET1_30msc2 <- ss_SET1_30msc2[1:7500,-(1:6)]
ss_SET1_30msc1_05 <- ss_SET1_30msc1_05[1:7500,-(1:8)]
ss_SET1_30msc2_05 <- ss_SET1_30msc2_05[1:7500,-(1:8)]
ss_SET1_30msc1_005 <- ss_SET1_30msc1_005[1:7500,-(1:8)]
ss_SET1_30msc2_005 <- ss_SET1_30msc2_005[1:7500,-(1:8)]

#### keep FST_2_1, PI_2_1, Sx1_1_0, Sx0_1_0, Sf_1_0

statobs <- statobs[,c(5,6,7,8,10)]
ss_SET1_30msc1 <- ss_SET1_30msc1[1:7500,c(5,6,7,8,10)]
ss_SET1_30msc2 <- ss_SET1_30msc2[1:7500,c(5,6,7,8,10)]
ss_SET1_30msc1_05 <- ss_SET1_30msc1_05[1:7500,c(5,6,7,8,10)]
ss_SET1_30msc2_05 <- ss_SET1_30msc2_05[1:7500,c(5,6,7,8,10)]
ss_SET1_30msc1_005 <- ss_SET1_30msc1_005[1:7500,c(5,6,7,8,10)]
ss_SET1_30msc2_005 <- ss_SET1_30msc2_005[1:7500,c(5,6,7,8,10)]

#create index of simulations

ss_SET1_30msc1$index <- rep("ss_SET1_30msc1", length(ss_SET1_30msc1$FST_2_1))
ss_SET1_30msc2$index <- rep("ss_SET1_30msc2", length(ss_SET1_30msc2$FST_2_1))
ss_SET1_30msc1_05$index <- rep("ss_SET1_30msc1_05", length(ss_SET1_30msc1_05$FST_2_1))
ss_SET1_30msc2_05$index <- rep("ss_SET1_30msc2_05", length(ss_SET1_30msc2_05$FST_2_1))
ss_SET1_30msc1_005$index <- rep("ss_SET1_30msc1_005", length(ss_SET1_30msc1_005$FST_2_1))
ss_SET1_30msc2_005$index <- rep("ss_SET1_30msc2_005", length(ss_SET1_30msc2_005$FST_2_1))

ss<-rbind(ss_SET1_30msc1,ss_SET1_30msc2,
          ss_SET1_30msc1_05,ss_SET1_30msc2_05,
          ss_SET1_30msc1_005,ss_SET1_30msc2_005)

grouping_scenarios=list(c('ss_SET1_30msc1', 'ss_SET1_30msc2'),
                        c('ss_SET1_30msc1_05','ss_SET1_30msc2_05'),
                        c('ss_SET1_30msc1_005','ss_SET1_30msc2_005'))

head(ss)

#SUMSTA

sumsta <-data.frame(ss)
nSS=as.numeric(ncol(ss))
ref_table<-ss[,c(1:(ncol(ss) - 1))]
head(ref_table)
ss$index <- as.factor(ss$index)

#index <- ss$index
##################################################################################
# Contributions to the random forests and LDA projections of the reference table #
##################################################################################
# Two graphics/figures providing:
# (i) the contributions of the 30 most important statistics to the RF (file named graph_varImpPlot.pdf).
# (ii) LDA projections of the reference table for the different scenarios plus the observed dataset (cf. black star in the figure; file named graph_lda.pdf).
# e.g. Fig. S6 and Fig. S7 in Pudlo et al. 2016.
mc.rf <- abcrf(formula=index~., data=ss,ntree=100) # Computation related to the forest classification.
plot(mc.rf, ss, obs=statobs[1,], pdf=TRUE)
plot(x=mc.rf, obs=statobs, training=ss, pdf=TRUE, n.var=13) # "pdf": a boolean that 
mc.rf
dev.off()

## change the LDA axe ##
#COLORS LDA
projections <- predict(mc.rf$model.lda, ss)$x
mf <- match.call(expand.dots=FALSE)
mf <- mf[1]
mf$formula <- mc.rf$formula
mf$data <- ss
mf[[1L]] <- as.name("model.frame")
mf <- eval(mf, parent.frame() )
mt <- attr(mf, "terms")
modindex <- model.response(mf)
nmod <- length(mc.rf$model.rf$forest$levels)
coloris <- rainbow(nmod)
colo <- coloris[modindex]

#OBS
projobs <- predict(mc.rf$model.lda, statobs)$x

#projections 1 3 
pdf("graph_LDA.pdf", width=6, height=6)
plot(x=projections[,1], y=projections[,3], xlab="LD1", ylab="LD3" , col=colo, pch=3, xlim = c(-10, 5), ylim = c(-2, 7))
legend("topright", legend = as.character(mc.rf$model.rf$forest$levels), col = coloris, 
       pch = 15, bty = "o", pt.cex = 1, cex = .8, horiz = FALSE, 
       inset = c(0, 0), ncol = 1, title = "Models", bg = "white")
points(projobs[1],projobs[3],pch="*",cex=5.3)
dev.off()

Nsimu_per_model=10000
ntrees_in_forest=500 	# Number of trees in the random forest (default=500; Pudlo et al. 2016).
nrep=10	# Number of replicates analyses used to calculated the average and the standard deviation (sd) of posterior probability (default=10, Fraimout et al. 2017). 
ncores_on_speed=20
nscenarios=3	# Number of compared scenarios (models) in the reference table.
nSS=as.numeric(ncol(ss))	# Size of the reference table (i.e. number of simulated datasets) that will be used to do the RF analysis.
Nref=nscenarios*Nsimu_per_model # NB OF SCENARIOS * NB OF SIMULATIONS #
##############################
# Random forest computations #
##############################
TotOut=list()		# Create a new empty list (necessary for the next loop).
TotOut=list()		# Create a new empty list (necessary for the next loop).
TotVote=list()		# Create a new empty list (necessary for the next loop).
TotBest=list()		# Create a new empty list (necessary for the next loop).
TotProb=list()		# Create a new empty list (necessary for the next loop).
T1 <- Sys.time()	# T1: starting time of the analysis (including reading the file reftable.txt).

## 1 Compare 30 scenarios together
for (i in 1:nrep) {
  mc.rf <- abcrf(formula=index~., data=ss, paral=TRUE, ntree=ntrees_in_forest, ncores=ncores_on_speed,grouping_scenarios)
  Out <- capture.output(mc.rf) 		# Send the console output (i.e. numerical results for the global prior error rates and the matrix of confusion (i.e. prior error rates detailed for each scenario) to a object "Out". Estimation from 10,000 pseudo-observed datasets from the reference table.
  TotOut[[i]] <- Out		# Add the numerical results for the global prior error rates and the matrix of confusion to a result list.
  Best <- predict(object=mc.rf, obs=statobs, training=ss, paral=TRUE, ntree=ntrees_in_forest)
  TotVote[[i]] <- Best$vote		# Add the pproportion of vote per scenario to a result list.
  TotBest[[i]] <- Best$allocation	# Add the best scenario to a result list.
  TotProb[[i]] <- Best$post.prob	# Add the posterior probability of the best scenario to a result list.
  T2 <- Sys.time()			# T2: ending time of the analysis.
  duration = difftime(T2, T1)		# Duration of the analysis (i.e. T2-T1).
  print(duration)				# 6.436006 mins
}

##################
# Number of vote #
##################

TOTVOTE <- do.call(rbind.data.frame, TotVote) # Convert the list of best scenario into a data frame (useful for the exportation).
write.csv(TOTVOTE, "PropVotes.csv") # Export results into a csv file.
NBVOTE <- TOTVOTE*500     # Calculate the number of vote per scenario by multiplying the proportion of vote to 500 (i.e. the total number of vote).
for (i in 1:nscenarios) { # Calculate the mean number of vote per scenario and the corresponding standard deviation (sd).
  NBVOTE[nrep+1, i] <- mean(as.numeric(NBVOTE[1:nrep, i]))
  NBVOTE[nrep+2, i] <- sd(as.numeric(NBVOTE[1:nrep, i]))
}
NBVOTE$Analysis <- c(1:nrep, "Mean", "sd") # Add a new column indicating (for each line) the corresponding analysis (e.g. "Mean"...).
write.csv(NBVOTE, "NbVotes.csv") # Export results into a csv file.

#################
# Best scenario #
#################

TOTBEST <- do.call(rbind.data.frame, TotBest) ; colnames(TOTBEST) <- "Best_Scen" # Convert the list of best scenario into a data frame (useful for the exportation).
write.csv(TOTBEST, "BestModel.csv") # Export results into a csv file.

######################################################################################
# Mean and standard deviation (sd) of the posterior probability of the best scenario #
######################################################################################

TOTPROB <- do.call(rbind.data.frame, TotProb) ; colnames(TOTPROB) <- "Post_Prob" # Convert the list of Post Prob into a data frame (useful for the exportation).
MEAN <- mean(TOTPROB[1:nrep, 1]) ; TOTPROB[nrep+1, 1]=MEAN # Computes the averaged posterior probability and add it to the data frame.
SD <- sd(TOTPROB[1:nrep, 1]) ; TOTPROB[nrep+2, 1]=SD       # Computes the standard deviation (sd) of the posterior probability and add it to the data frame.
TOTPROB$Analysis <- c(1:nrep, "Mean", "sd")    # Add a new column indicating (for each line) the corresponding analysis (e.g. "Mean"...).
TOTPROB <- TOTPROB[c("Analysis", "Post_Prob")] # Reorder the column (i.e. 1st: "Analysis"; 2nd: "Post_Prob").
write.csv(TOTPROB, "PostProb.csv")             # Export results into a csv file.

############################################################
# Mean and standard deviation (sd) of the prior error rate #
############################################################

TOTOUT_A <- do.call(rbind.data.frame, TotOut) # Convert the list of results into a data frame (useful for the exportation).
TOTOUT_B <- TOTOUT_A[, 7]                     # Export the column containing the prior error rate (i.e. the 7th column).
TOTOUT_C <- gsub(x=TOTOUT_B, pattern="Out-of-bag prior error rate: ", replacement=""); TOTOUT_C <- gsub(x=TOTOUT_C, pattern="%", replacement="") # Remove the text into the column (keep only the numerical values).
PRIOR_ERROR_RATE <- data.frame(as.numeric(TOTOUT_C)) ; colnames(PRIOR_ERROR_RATE) <- "Prior_Error_Rate" # Create a data frame containing the prior error rate for each analysis. 
PRIOR_ERROR_RATE_MEAN <- mean(PRIOR_ERROR_RATE[0:nrep, 1]) ; PRIOR_ERROR_RATE[nrep+1, 1]=PRIOR_ERROR_RATE_MEAN # Computes the averaged prior error rate and add it to the data frame.
PRIOR_ERROR_RATE_SD <- sd(PRIOR_ERROR_RATE[0:nrep, 1]) ; PRIOR_ERROR_RATE[nrep+2, 1]=PRIOR_ERROR_RATE_SD       # Computes the standard deviation (sd) of the prior error rate and add it to the data frame.
PRIOR_ERROR_RATE$Analysis <- c(1:nrep, "Mean", "sd") # Add a new column indicating (for each line) the corresponding analysis (e.g. "Mean"...).
PRIOR_ERROR_RATE <- PRIOR_ERROR_RATE[c("Analysis", "Prior_Error_Rate")] # Reorder the column (i.e. 1st: "Analysis"; 2nd: "Prior_Error_Rate").
write.csv(PRIOR_ERROR_RATE, "PriorErrorRate.csv") # Export results into a csv file.

###############################
# Combine all tables into one #
###############################

ALL <- NBVOTE                                             # Create a new dataframe based on the dataframe containing the number of vote per scenario.
ALL$Post_Prob <- TOTPROB$Post_Prob                        # Add a new column indicating the posterior probability.
ALL$Prior_Error_Rate <- PRIOR_ERROR_RATE$Prior_Error_Rate # Add a new column indicating the posterior probability.
write.csv(ALL, "ABCRF_Results_Step1.csv")                 # Export results into a csv file.

