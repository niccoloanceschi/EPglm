rm(list=ls())
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #set the directory to the path of the current script
library("rstan")
library("bpr")
source("EPGLM_fcts.R")

######################################################
# DATA LOADING
######################################################
load(file="Data-and-Results/football_large_p.RData")
X <- as.matrix(data2[,-1])
y <- data2[,1]
p <- ncol(X)

# Hyper-parameters -------------------------------------------------------------

seed = 123          # seed for random number generation 

tolerance = 1e-3          # tolerance for stopping iterative updates

nSample = 1e4       # number of MCMC samples

om2p     =  25*10/p  # marginal prior variances
Omega0   = diag(om2p,p,p)
beta0    = rep(0.,p)


# Split train-test -------------------------------------------------------------
set.seed(seed)
n        = length(y)
nTest    = 60
nTrain   = n-nTest

labTest = sample.int(n,nTest)

X_test = X[labTest,]
X      = X[-labTest,]
y_test = y[labTest]
y      = y[-labTest]

######################################################
# DEFINE POISSON MODEL IN STAN
######################################################
burn = 1e3
nChains_mcmc = 5
nSample_stan = as.integer(nSample/nChains_mcmc+burn) # so we have nSample samples after burnin

poissonGLM_diagPrior = 'data{
                      	  int<lower=0> K;
                      	  int<lower=0> N;
                      	  int<lower=0> Y[N];
                      	  matrix[N,K] X;
                      	  real<lower=0> devStd;
                        }
                        parameters {
                        	vector[K] beta;
                        }
                        model {
                          for(j in 1:K)
                        	  beta[j] ~ normal(0,devStd);
                        	for(i in 1:N)
                        	  Y[i] ~ poisson_log(X[i]*beta);
                        }'


######################################################
# GET HMC SAMPLES
######################################################  
data_stan = list(N=nTrain,K=p,Y=as.vector(y),X=X,devStd=sqrt(om2p))

# If RUN=T run the code
RUN=F
if(RUN){
  startTime = Sys.time()
  HMC_Samples = stan(model_code = poissonGLM_diagPrior, data = data_stan, iter = nSample_stan,
                     warmup=burn, chains = nChains_mcmc, init="0", algorithm="NUTS", seed=seed, cores = nChains_mcmc)
  
  betaHMC <- t(extract(HMC_Samples)$beta)
  meanBetaHMC = apply(betaHMC,1,mean)
  sdBetaHMC = apply(betaHMC,1,sd)
  
  timeHMC = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  # prediction
  predMeanHMC = rowMeans(exp(X_test%*%betaHMC))
  
  timeHMCpredMean = difftime(Sys.time(), startTime, units=("secs"))[[1]] # timeMCMC_algo[i] + difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  
  ######################################################
  # GET EP POSTERIOR MOMENTS
  ######################################################
  startTime = Sys.time()
  paramsEP = getParamsEP(X,y,family='poisson',Omega0=Omega0,beta0=beta0,approxHybrid=T,
                         fullVar=TRUE,damp=1.,nPrint=100,tolerance=tolerance,maxIter=1e5)
  timeEP = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  # prediction
  m_pred <- X_test%*%paramsEP$meanBeta
  s2_pred <- rowSums(X_test*(X_test%*%paramsEP$Omega))
  predMeanEP    = exp(m_pred+s2_pred/2)
  timeEPpredMean = difftime(Sys.time(), startTime, units=("secs"))[[1]]

  # store quantities
  meanBetaEP = paramsEP$meanBeta
  sdBetaEP = sqrt(paramsEP$diagOmega)

  # save(meanBetaHMC,meanBetaEP,
  #      sdBetaHMC,sdBetaEP,
  #      predMeanHMC,predMeanEP,
  #      timeHMC,timeEP,
  #      timeHMCpredMean,timeEPpredMean,
  #      file = "Data-and-Results/outputIllustrationPoisson_large_p_small_n.RData")
} else {
  load("Data-and-Results/outputIllustrationPoisson_large_p_small_n.RData")
}

library(ggplot2)
library(reshape2)
library(scales)


colnames(meanBetaEP) = c("EP")
meanData = melt(meanBetaEP)[,-1]
meanData$group1 = c("Expectations")
meanData$y = c(meanBetaHMC)

sdBetaEP = as.matrix(sdBetaEP, ncols = 1)
colnames(sdBetaEP) = c("EP")
sdData = melt(sdBetaEP)[,-1]
sdData$group1 = c("Standard Deviations")
sdData$y = c(sdBetaHMC)

colnames(predMeanEP) = c("EP")
predData = melt(predMeanEP)[,-1]
predData$group1 = c("Predictive Means")
predData$y = c(predMeanHMC)

data_points = rbind(meanData,sdData,predData)
data_points$group1 = factor(data_points$group1, levels=c("Expectations", "Standard Deviations","Predictive Means" ))
P_Pois = ggplot(data_points, aes(y=value, x=y,color=c("red"))) +
  geom_point(shape=21,fill="red",alpha=0.5)+facet_wrap(group1~., scales="free")+
  geom_abline(slope = 1, intercept = 0, alpha = 0.7, lty = 2)+theme_bw()+
  labs(x="Quantities computed from the exact posterior", y = "Quantities computed from the EP posterior")+
  theme(strip.text.x = element_text(size = 15),axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15), axis.text.y = element_text(size=15),
        plot.margin = margin(0.1, 0.1, 0.05, 0.2, "cm"),legend.position = "none")+
  scale_color_manual(values=c("red"))+
  scale_shape_manual(values=c(1))+
  scale_fill_manual(values=c("red"))

show(P_Pois)


