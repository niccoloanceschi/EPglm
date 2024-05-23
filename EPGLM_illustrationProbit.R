rm(list=ls())
library(rstudioapi)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #set the directory to the path of the current script
source("EPGLM_fcts.R")
library("TruncatedNormal")
library("truncnorm")
library("transport")


######################################################
# LOAD THE DATA
######################################################
load("Data-and-Results/Alzheimer_Interactions.RData")

seed = 1
set.seed(seed)
n = 300 # training set size

# split training and test sets
X_Test = X[-trainingSet,]
yTest = y[-trainingSet]
X = X[trainingSet,]
y = y[trainingSet]


# get number of covariates
p = dim(X)[2]

# prior variance
nu2 = 25

# fix number of samples 
nSample = 1e4

# If RUN=T run the code
RUN=F
if(RUN){
  # initialize vectors for predictive probabilities
  
  predProbPFM    = double(length = length(yTest))
  predProbEP     = double(length = length(yTest))
  predProbMC     = double(length = length(yTest))
  
  ######################################################
  # PFM
  ######################################################
  tolerance = 1e-3 # tolerance to establish ELBO convergence
  startTime = Sys.time()
  paramsPFM = getParamsPFM(X=X,y=y,nu2=nu2,moments=TRUE,predictive=TRUE,tolerance=tolerance,maxIter=1e4)
  timePFM = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  # posterior moments
  meanPFM = paramsPFM$postMoments.meanBeta
  sdPFM = sqrt(paramsPFM$postMoments.varBeta)
  
  # predictive probabilities
  # obtain a sample from q^*(z) to be used for the predictive probabilities of PFM
  nSampleZ = 1e4
  muTN = paramsPFM$mu
  muTN[y==0] = -muTN[y==0]
  sampleTruncNorm = matrix(rtruncnorm(n*nSampleZ, a = 0, b = Inf, mean = muTN, sd = sqrt(paramsPFM$sigma2)), nrow = n, ncol = nSampleZ, byrow = F )
  sampleTruncNorm[y==0,] = -sampleTruncNorm[y==0,] # need to adjust the sign of the variables for which y_i is 0
  
  for(j in 1:length(yTest)){
    xNew = matrix(X_Test[j,],ncol = 1)
    Xx = X%*%xNew
    sd = as.double(sqrt(1+nu2*(sum(xNew^2)-nu2*t(Xx)%*%paramsPFM$predQuant.invIXXt%*%Xx)))
    
    predProbPFM[j] = mean(pnorm((t(xNew)%*%paramsPFM$predQuant.VXt%*%sampleTruncNorm)/sd))
  }
  timePFMpredProb = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  ######################################################
  # EP
  ######################################################
  Omega0   = diag(nu2,p,p)
  beta0    = rep(0.,p)
  
  startTime = Sys.time()
  tolerance = 1e-3 # tolerance to establish convergence
  paramsEP = getParamsEPprobit(X=X,y=y,nu2=nu2,tolerance=tolerance,predictive=T,fullVar = F,maxIter=1e4) 
  # paramsEP = getParamsEPprobit(X,y,family='probit',Omega0=Omega0,beta0=beta0,approxHybrid=T,fullVar=TRUE,damp=1.,nPrint=100,tolerance=tolerance,maxIter=1e5)
  timeEP = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  # posterior moments
  meanEP = paramsEP$meanBeta
  sdEP = sqrt(paramsEP$diagOmega)
  
  # predictive probabilities
  for(j in 1:length(yTest)){
    xNew = matrix(X_Test[j,],ncol = 1)
    predProbEP[j] = as.double(predictEPprobit(paramsEP,xNew,nu2))
  }
  timeEPpredProb = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  ######################################################
  # EXACT SAMPLING
  ######################################################
  startTime = Sys.time()
  betaMC = rSUNpost(X=X,y=y,nu2=nu2,nSample=nSample)
  
  # posterior moments
  meanMC = apply(betaMC,1,mean)
  sdMC = apply(betaMC,1,sd)
  timeMC = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  # predictive probabilities
  for(j in 1:length(yTest)){
    xNew = matrix(X_Test[j,],ncol = 1)
    predProbMC[j] = mean(pnorm(t(xNew)%*%betaMC))
  }
  timeMCpredProb = difftime(Sys.time(), startTime, units=("secs"))[[1]]
  
  ######################################################
  # CHECK RUNNING TIMES AND NUMBER OF ITERATIONS
  ######################################################
  timeMC # in second
  timeMCpredProb # in seconds
  timePFM # in seconds
  timePFMpredProb # in seconds
  timeEP # in seconds
  timeEPpredProb # in seconds
  
  paramsPFM$nIter
  paramsEP$nIter
  
  ######################################################
  # COMPUTE AND SAVE THE NEEDED SUMMARY STATISTICS
  ######################################################
  
  # save(meanMC,meanPFM,meanEP,
  #      predProbMC,predProbPFM,predProbEP,
  #      sdMC,sdPFM,sdEP,
  #      timeMC,timePFM,timeEP,
  #      timeMCpredProb,timePFMpredProb,timeEPpredProb,
  #      file = "Data-and-Results/outputIllustrationProbit.RData")
} else {
  load("Data-and-Results/outputIllustrationProbit.RData")
} 

library(ggplot2)
library(reshape)
library(RColorBrewer)

means = cbind(meanEP,meanPFM)
colnames(means) = c("EP","PFM-VB")
meanData = melt(means)[,-1]
meanData$group1 = c("Expectations")
meanData$y = c(meanMC,meanMC)

sds = cbind(sdEP,sdPFM)
colnames(sds) = c("EP","PFM-VB")
sdData = melt(sds)[,-1]
sdData$group1 = c("Standard Deviations")
sdData$y = c(sdMC,sdMC)

pred = cbind(predProbEP,predProbPFM)
colnames(pred) = c("EP","PFM-VB")
predData = melt(pred)[,-1]
predData$group1 = c("Predictive Probabilities")
predData$y = c(predProbMC,predProbMC)

######################################################################
# BOXPLOTS
######################################################################
### Comparison different methods via box-plot
EP_Diff_mean_Beta  = meanMC - meanEP
PFM_Diff_mean_Beta = meanMC - meanPFM

EP_Diff_logsd_Beta  = log(sdMC) - log(sdEP)
PFM_Diff_logsd_Beta = log(sdMC) - log(sdPFM)

EP_Diff_sd_Beta  = sdMC - sdEP
PFM_Diff_sd_Beta = sdMC - sdPFM

EP_Diff_predProb  = predProbMC - predProbEP
PFM_Diff_predProb = predProbMC - predProbPFM


Diff_Mean  = c(EP_Diff_mean_Beta, PFM_Diff_mean_Beta)
Diff_logsd = c(EP_Diff_logsd_Beta, PFM_Diff_logsd_Beta)
Diff_predProb = c(EP_Diff_predProb, PFM_Diff_predProb)


#Par        = as.factor(rep(rep(c(rep(1,n),rep(2,n)),3),2))
Method     = as.factor(
  c(rep(rep(c("EP","PFM-VB"),each=length(EP_Diff_mean_Beta)),2),
    rep(c("EP","PFM-VB"),each=length(EP_Diff_predProb)))
)
Functional = as.factor(c(rep(1:2,each=2*length(EP_Diff_mean_Beta)),
                         rep(3,2*length(EP_Diff_predProb)))
)

Diff       = c(Diff_Mean,Diff_logsd,Diff_predProb)
library("viridis")
library("latex2exp")
Data_boxplot = data.frame(Diff,Method,Functional)
scale_color_viridis(discrete=TRUE)
col_meth   = c("red","blue")
levels(Data_boxplot$Functional) = 
  c("1"=TeX("$\\Delta E(\\beta_{j} | y)$"), 
    "2"=TeX("$\\Delta \\log(\\sqrt{var(\\beta_{j} | y)})$"),
    "3"=TeX("$\\Delta P(y_{n+1}=1 | y)$")
  )
Data_boxplot$Method = factor(Data_boxplot$Method, c("EP","PFM-VB"))

Plot_Diff = ggplot(Data_boxplot, aes(y=Diff, x=Method,col=Method)) +
  geom_boxplot()+
  scale_colour_manual(values = alpha(col_meth,c(1,1)))+
  theme_bw()+ geom_hline(yintercept=0, linetype="dashed", linewidth=2)+
  theme(axis.title=element_blank(),axis.text=element_text(size=15) ,
        strip.text = element_text(size=15),legend.position = "none")+
  facet_wrap(.~Functional, scales = "free", labeller=label_parsed)

Plot_Diff





