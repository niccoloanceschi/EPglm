getParamsPFM = function(X , y ,nu2, moments = TRUE, predictive = FALSE, tolerance = 1e-2, maxIter = 1e4) {
  ###################################################### -
  # PRECOMPUTATION
  ###################################################### -
  # define model dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  
  # compute H = X%*%V%*%t(X) and Omega_z directly or with Woodbury
  if(p<=n) {
    # define prior covariance matrix and its inverse
    # Omega = diag(rep(nu2,p),p,p)
    invOmega = diag(rep(1/nu2,p),p,p)
    V = solve(crossprod(X)+invOmega)
    H = X%*%V%*%t(X)
    invOmZ = diag(1,nrow=n,ncol=n) - H # needed for ELBO
  } else{
    XXt = tcrossprod(X)
    invOmZ = solve(diag(1,nrow=n,ncol=n)+nu2*XXt) # needed for ELBO
    H = nu2*XXt%*%invOmZ
  }
  
  # compute optimal sigma2
  # h = diag(diag(H))
  sigma2 = as.double(1/(1-diag(H)), ncol = 1)
  sigma = sqrt(sigma2)
  
  # compute matrix to write the CAVI update in a vectorized form
  A = sigma2*H
  A[cbind(1:n,1:n)] = 0
  
  # other useful quantities needed for ELBO
  diagInvOmZ = diag(invOmZ)
  
  # initialization of variables
  mu = matrix(0,n,1)
  
  # inizialize coherently the vector of expectations meanZ
  musiRatio = as.double(mu/sigma)
  phiPhiRatio = exp(dnorm(musiRatio,log = T)-pnorm((2*y-1)*musiRatio,log.p = T))
  meanZ = mu + (2*y-1)*sigma*phiPhiRatio
  
  elbo = -Inf
  diff = 1
  nIter=0
  
  ###################################################### -
  # CAVI ALGORITHM
  ###################################################### -
  
  while(diff > tolerance & nIter < maxIter) {
    elboOld = elbo
    sumLogPhi = 0
    
    for(i in 1:n) {
      mu[i] = A[i,]%*%meanZ
      
      # compute first (needed for algorithm) and second (needed for ELBO) moments
      musiRatio = mu[i]/sigma[i]
      phiPhiRatio = exp(dnorm(musiRatio, log = T) - pnorm((2*y[i]-1)*musiRatio, log.p = T))
      meanZ[i] = mu[i] + (2*y[i]-1)*sigma[i]*phiPhiRatio
      sumLogPhi = sumLogPhi + pnorm((2*y[i]-1)*musiRatio, log.p = T)
    }
    
    # computation of ELBO (up to an additive constant not depending on mu)
    elbo = -(t(meanZ)%*%invOmZ%*%meanZ -
               sum((meanZ^2)*diagInvOmZ))/2 -
      sum(meanZ*mu/sigma2) + sum((mu^2)/sigma2)/2 + sumLogPhi
    
    diff = abs(elbo-elboOld)
    nIter = nIter+1
    
    if(nIter%%100==0) {print(paste0("iter: ", nIter, ", ELBO: ", elbo))}
  }
  
  # get the optimal parameters of the normals before truncation, now that convergence has been reached
  mu = A%*%meanZ
  
  results = list(mu = mu, sigma2 = sigma2, nIter = nIter)
  
  ###################################################### -
  # (OPTIONAL) CLOSED-FORM MOMENTS' COMPUTATION
  ###################################################### -
  
  if(moments == TRUE) {
    # compute V and V%*%t(X), directly or with Woodbury
    if(p<=n) {
      diagV = diag(V) # V already computed
      VXt = V%*%t(X)
    } else{ # use Woodbury
      VXt = t(nu2*X)%*%invOmZ
      #V = Omega - VXt%*%(nu2*X)
      diagV = nu2*(1-colSums(t(VXt) * X))
    }
    
    musiRatio = mu/sigma
    phiPhiRatio = exp(dnorm(musiRatio, log = T) - pnorm((2*y-1)*musiRatio, log.p = T))
    
    meanZ = mu + (2*y-1)*sigma*phiPhiRatio
    postVarZ = as.double(sigma2*(1-(2*y-1)*musiRatio*phiPhiRatio - phiPhiRatio^2))
    
    W = apply(VXt,1,function(x) sum(x*x*postVarZ))
    
    meanBeta = VXt%*%meanZ
    varBeta = diagV + W
    
    moments_PFM = list(meanBeta=meanBeta,varBeta=matrix(varBeta,ncol = 1))
    
    results = c(results,postMoments=moments_PFM)
  }
  
  if(predictive==T){
    # return also VXt and InvOmZ which are needed for the computation of predictive probabilities
    if(moments==F){
      # we weed to compute the quantities
      if(p<=n) {
        VXt = V%*%t(X)
      } else{
        VXt = t(nu2*X)%*%invOmZ
      }
    }
    # now we have quantities we need
    predQuant = list(VXt=VXt,invIXXt=invOmZ)
    results = c(results,predQuant=predQuant)
  }
  
  return(results)
}

rSUNpost = function(X,y,nu2,nSample) {
  # define model dimensions
  n = dim(X)[1]
  p = dim(X)[2]
  
  # get parameters useful for sampling
  Omega = diag(rep(nu2,p),p,p)
  invOmega = diag(rep(1/nu2,p),p,p)

  signX = X
  signX[y==0,] = -X[y==0,]
  Gamma_post_unnormalized = diag(1,n,n)+(nu2*signX)%*%t(signX)
  inv_Gamma_post_unnormalized = solve(Gamma_post_unnormalized)
  s = diag(sqrt(Gamma_post_unnormalized[cbind(1:n,1:n)]),n,n)
  s_1 = diag(1/s[cbind(1:n,1:n)],n,n)
  gamma_post = matrix(0,n,1) # because prior mean is set to 0
  Gamma_post = s_1%*%Gamma_post_unnormalized%*%s_1
  
  V = Omega-t(nu2*signX)%*%inv_Gamma_post_unnormalized%*%(signX*nu2)
  V = 0.5*(V+t(V))
  L = t(chol(V))
  
  # compute multiplicative coefficients for the truncated multivariate normal component
  coefTruncNorm = t(nu2*signX)%*%inv_Gamma_post_unnormalized%*%s
  
  # sample the multivariate normal component
  sampleMultNorm = matrix(rnorm(nSample*p),p,nSample)
  
  # sample the truncated multivariate normal component
  # if(n == 1) {
  #   sampleTruncNorm = matrix(rtruncnorm(n = nSample, a = -gamma_post, b = Inf, mean = 0, sd = 1), nrow = 1, ncol = nSample)
  # } else{
  #   sampleTruncNorm = mvrandn(l = -gamma_post, u = rep(Inf,n), Sig = Gamma_post, n = nSample)
  # } # this part is old now, after the update of TruncatedNormal to version 2.0
  sampleTruncNorm = t(rtmvnorm(n = nSample, mu = rep(0,n), sigma = Gamma_post, lb = -gamma_post, u = rep(Inf,n)))
  
  # combine the multivariate normal and truncated normal components
  sampleSUN = L%*%sampleMultNorm+coefTruncNorm%*%sampleTruncNorm
}



# ------------------------------- EFFICIENT EP FOR PROBIT ---------------------------------
zeta1 = function(x){exp(dnorm(x, log = T) - pnorm(x,log.p = T))}
zeta2 = function(x,z1){-z1^2-x*z1}

getParamsEPprobit = function(X,y,nu2,
                             tolerance=1e-3,maxIter=1e3,nPrint=100,
                             fullVar=FALSE,predictive=FALSE,
                             force_highdim_algo=FALSE){

  n = dim(X)[1]
  p = dim(X)[2]
  
  ### Initialization
  if(p<n && force_highdim_algo == F){
    invQ = diag(nu2,p,p)
  }else{
    V = nu2*t(X)
  }

  r = rep(0,p)

  k = double(length = n)
  m = double(length = n)

  diff = 1
  nIter = 0

  ### Iterations

  if(p<n && force_highdim_algo == F){

    while(diff > tolerance && nIter < maxIter){

      diff = 0.
      count = 0

      for(i in 1:n){

        xi = X[i,]

        r_i = r - m[i]*xi

        Oxi = invQ%*%xi

        Oi = invQ + tcrossprod(Oxi)*k[i]/as.double(1.-k[i]*crossprod(xi,Oxi))

        Oixi = Oi%*%xi
        xiOixi = as.double(crossprod(xi,Oixi))

        if(xiOixi>0){

          r_iOixi = as.double(crossprod(r_i,Oixi))

          s = (2.*y[i]-1.)/sqrt(1.+xiOixi)
          tau = s*r_iOixi

          z1 = zeta1(tau)
          z2 = zeta2(tau,z1)

          kNew = - z2/(1.+xiOixi+z2*xiOixi)
          mNew = s*z1 + kNew*r_iOixi + kNew*s*z1*xiOixi

          maxDiff = max(abs(c(kNew - k[i], mNew - m[i])))
          if(maxDiff>diff){diff = maxDiff} #???

          k[i] = kNew
          m[i] = mNew

          r = r_i + m[i]*xi

          invQ = Oi - tcrossprod(Oixi)*k[i]/(1.+k[i]*xiOixi)
        }else{
          count = count+1
          print(paste0(count," units skipped"))
        }

      }

      nIter = nIter + 1
      if(nIter %% nPrint == 0) {print(paste0("iteration ",nIter))}
    }

  } else {

    while(diff > tolerance && nIter < maxIter){

      diff = 0.
      count = 0

      for(i in 1:n){

        v = V[,i]
        xi = X[i,]
        xTv = crossprod(xi,v)[1]

        d = 1-k[i]*xTv
        w = v/d
        xTw = xTv/d

        if(xTw>0){

          r_iTw = crossprod(r,w) - m[i]*xTw

          s = (2*y[i]-1)/sqrt(1+xTw)
          tau = s*r_iTw

          z1 = zeta1(tau)
          z2 = zeta2(tau,z1)

          kNew = as.double(-z2/(1 + xTw + z2*xTw))
          mNew = as.double(z1*s + kNew*r_iTw + kNew*z1*s*xTw)

          r = r + (mNew - m[i])*xi

          ratio = (k[i]-kNew)/(1.+(kNew-k[i])*xTv)
          V = V + (v*ratio)%*%crossprod(xi,V)

          maxDiff = max(abs(c(kNew - k[i], mNew - m[i])))
          if(maxDiff>diff){diff = maxDiff}

          k[i] = kNew
          m[i] = mNew

        }else{
          count = count+1
          print(paste0(count," units skipped"))
        }
      }

      nIter = nIter + 1
      if(nIter %% nPrint == 0) {print(paste0("iteration ",nIter))}
    }

  }

  ### Posterior Approximate Moments

  if(p<n && force_highdim_algo == F){

    meanBeta = invQ%*%r
    diagOmega = diag(invQ)

  }else{

    diagOmega = rep(nu2,p)*(rep(1,p) - rowSums(V*t(k*X)))
    meanBeta = nu2*(r - V%*%(k*X%*%r))
  }

  results = list(meanBeta = meanBeta, diagOmega = diagOmega,
                 nIter = nIter, kEP = k, mEP = m)

  if(fullVar==TRUE){
    if(p>=n || force_highdim_algo == T) {
      invQ = nu2*(diag(1,p,p) - V%*%(k*X))
    }
    results = append(list(Omega=invQ),results)
  }

  if(predictive==TRUE){
    if(p>=n || force_highdim_algo == T){
      results = append(list(V=V),results)
    } else{
      if(fullVar==FALSE){
        results = append(list(Omega=invQ),results)
      }
    }
  }

  return(results)
}

predictEPprobit = function(paramsEP,xNew,nu2,force_highdim_algo=FALSE){
  
  p = length(paramsEP$meanBeta)
  n = length(paramsEP$kEP)

  if(p>=n || force_highdim_algo==T){
    Xx = X%*%xNew
    KVtx = paramsEP$k*(t(paramsEP$V)%*%xNew)
    sd = as.double(sqrt(1+nu2*(sum(xNew^2)-sum(KVtx*Xx))))
  }else{
    sd = as.double(sqrt(1+t(xNew)%*%paramsEP$Omega%*%xNew))
  }
  predProb = as.double(pnorm(t(xNew)%*%paramsEP$meanBeta/sd))
  return(predProb)
}

# ------------------------------- EFFICIENT EP FOR GLM ---------------------------------

# -------------------------- PROBIT moment matching ----------------------------
momentMatchingPROBIT <- function(y_i,m_i,s2_i){
  
  sP = (2.*y_i-1.)/sqrt(1.+s2_i)
  tau = sP*m_i
  
  logZnew = pnorm(tau,log.p = T)
  
  zeta1 = exp(dnorm(tau, log = T) - pnorm(tau,log.p = T))
  zeta2 = -zeta1^2-tau*zeta1
  
  precNew = -zeta2/(1.+s2_i+zeta2*s2_i)
  locNew = sP*zeta1 + precNew*m_i + precNew*sP*zeta1*s2_i
  
  list(precNew=precNew,locNew=locNew,logZnew=logZnew)
}


# ----------------------- LOGIT approx moment matching -------------------------
momentMatchingLOGIT <- function(y_i,m_i,s2_i){
  
  sP = 1./sqrt(1. + s2_i*pi/8.)
  Znew = plogis((2.*y_i-1.)*m_i*sP)
  
  dZdm_Z = (1.-Znew)*(2.*y_i-1.)*sP
  dZds2_Z = -0.5*(pi/8.)*m_i*sP^2*dZdm_Z
  
  Delta = dZdm_Z^2 - 2.*dZds2_Z
  
  precNew = Delta / (1.-s2_i*Delta)
  locNew = (dZdm_Z + m_i*Delta) / (1.-s2_i*Delta)
  
  list(precNew=precNew,locNew=locNew,logZnew=log(Znew))
}


# ---------------------- POISSON approx moment matching ------------------------

Hermite06 <- function(x){
  c(1,2*x,4*x^2-2,8*x^3-12*x,16*x^4-48*x^2+12,
    32*x^5-160*x^3+120*x,64*x^6-480*x^4+720*x^2-120)
}
erf  <- function(x) { 2 * pnorm(x * sqrt(2)) - 1 }
erfc <- function(x) { 1 - erf(x) }

momentMatchingPOISSONapprox_Rossberg <- function(y_i,m_i,s2_i){
  
  s_i = sqrt(s2_i)
  tau_i = -(m_i+y_i*s2_i)/s_i
  
  K0 = 6
  alpha = 1-(2^(-c(1:K0)-1))*c(1.6911373, 0.0875520, 1.4803348, 0.5847012, 1.1523156, 0.8769133)
  
  vecH <- Hermite06(tau_i/sqrt(2))
  vecC <- (alpha[K0] - alpha[-K0]) / ((s_i * sqrt(2))^(1:(K0-1)))
  
  G_i <- 0.5 + 0.5 * erf(tau_i/sqrt(2)) + exp(-tau_i^2/2)*sum(vecC*vecH[-c(K0,K0+1)])/sqrt(pi) -
    alpha[K0] * 0.5 * exp(s2_i/2 - s_i*tau_i) * erfc((s_i - tau_i)/sqrt(2))
  
  dGdy <- exp(-tau_i^2/2)/sqrt(2*pi) +
    alpha[K0]*0.5*exp(-tau_i^2/2)*(exp((s_i - tau_i)^2/2)*s_i*erfc((s_i - tau_i)/sqrt(2)) - sqrt(2/pi)) -
    exp(-tau_i^2/2)*tau_i*vecC[1]/sqrt(pi) -
    exp(-tau_i^2/2)*sum(vecC[-1]*vecH[-c(1,2,K0)])/sqrt(2*pi) 
  
  d2Gdy2 <- -tau_i*exp(-tau_i^2/2)/sqrt(2*pi) -
    alpha[K0]*0.5*exp(-tau_i^2/2)*(exp((s_i - tau_i)^2/2)*s_i^2*erfc((s_i - tau_i)/sqrt(2)) - (s_i + tau_i)*sqrt(2/pi)) - 
    (1-tau_i^2)*exp(-tau_i^2/2)*(alpha[K0] - alpha[1])/(sqrt(pi)*(s_i*sqrt(2))) +
    exp(-tau_i^2/2)*sum(vecC[-1]*vecH[-c(1,2,3)])/sqrt(2*pi) 
  
  Theta_i <- y_i - (dGdy/G_i)/s_i;
  Delta_i <- (y_i^2 - 2*y_i*(dGdy/G_i)/s_i + (d2Gdy2/G_i)/s2_i) - Theta_i^2 ;
  
  # varNew = s2_i + s2_i*Delta_i*s2_i
  # meanNew = m_i + s2_i*Theta_i
  
  logZnew <- y_i*(m_i + y_i*s2_i/2) - log(factorial(y_i)) + log(G_i) 
  
  precNew = - Delta_i / (1.+ s2_i*Delta_i)
  locNew = precNew * (m_i - Theta_i / Delta_i) 
  
  list(precNew=precNew,locNew=locNew,logZnew=logZnew)
}

momentMatchingPOISSONapprox_Asmussen <- function(y_i,m_i,s2_i){
  
  tau = m_i+y_i*s2_i
  
  if(tau>10){
    lWf = log(s2_i)+tau - log(log(s2_i)+tau) + 
      log(log(s2_i)+tau)/(log(s2_i)+tau) 
  } else{
    lWf = lamW::lambertW0(s2_i*exp(tau))    ### requires lamW package
  }
  
  logZnew = y_i*(m_i + y_i*s2_i/2) - 0.5*lWf^2/s2_i - lWf/s2_i - 
    0.5*log(1+lWf) - log(factorial(y_i))
  
  dZdm_Z = y_i - ((1+lWf)/s2_i+0.5/(1+lWf))*lWf/(1+lWf)
  dZds2_Z = 0.5*y_i^2 + (lWf/s2_i)*(1+0.5*lWf)/s2_i - 
    ((1+lWf)/s2_i+0.5/(1+lWf))*(lWf/s2_i)*(1+y_i*s2_i)/(1+lWf)
  
  Delta_ds2 = dZdm_Z^2 - 2.*dZds2_Z
  Delta_dm2 = (1/s2_i + (0.5-lWf)/(1+lWf)^3)*lWf/(1+lWf) # = dZdm_Z^2 - dZdm2_Z
  Delta = 0.5*(Delta_ds2+Delta_dm2)
  
  # varNew = s2_i - s2_i*Delta*s2_i
  # meanNew = m_i + s2_i*dZdm_Z
  
  precNew = Delta / (1.-s2_i*Delta)
  locNew = (dZdm_Z + m_i*Delta) / (1.-s2_i*Delta)
  
  list(precNew=precNew,locNew=locNew,logZnew=logZnew)
}

momentMatchingPOISSONapprox <- function(y_i,m_i,s2_i){
  
  # Empirical region of good approximation via Rossberg et al. (2007) for y_i=0
  y1 = 2; x1 = -20 + (3/9)*40; y2 = 10; x2 = -20 + (6/9)*40
  m = (y1 - y2)/(x1 - x2); q = (x1*y2 - x2*y1)/(x1 - x2)
  
  if(y_i==0 & (s2_i > q + m*m_i)) {
    momentMatchingPOISSONapprox_Rossberg(y_i,m_i,s2_i)
  } else {
    momentMatchingPOISSONapprox_Asmussen(y_i,m_i,s2_i)
  }
}

# ----------------------- GAMMA approx moment matching -------------------------
momentMatchingGAMMA <- function(y_i,m_i,s2_i,nu_g){
  
  tau = nu_g*s2_i-m_i
  
  if(tau>10){
    lWf = log(nu_g*y_i)+tau - log(log(nu_g*y_i)+tau) + 
      log(log(nu_g*y_i)+tau)/(log(nu_g*y_i)+tau) 
  } else{
    lWf = lambertW0(nu_g*y_i*exp(tau))    ### requires lamW package
  }
  
  logZnew = dgamma(y_i,nu_g,rate=nu_g,log=T) - 0.5*log(1+lWf) + nu_g*y_i - 
    nu_g*(m_i-nu_g*s2_i/2) - 0.5*lWf^2/s2_i - lWf/s2_i 
  
  dZdm_Z = - nu_g + ((1+lWf)/s2_i+0.5/(1+lWf))*lWf/(1+lWf)
  dZds2_Z = 0.5*nu_g^2 + (1+0.5*lWf)*lWf/s2_i^2 -
    ((1+lWf)/s2_i+0.5/(1+lWf))*(lWfp1/(1.+lWf))*(1.+nu_g*s2_i)/s2_i
  
  Delta = dZdm_Z^2 - 2.*dZds2_Z
  
  # varNew = s2_i - s2_i*Delta*s2_i
  # meanNew = m_i + s2_i*dZdm_Z
  
  precNew = Delta / (1.-s2_i*Delta)
  locNew = (dZdm_Z + m_i*Delta) / (1.-s2_i*Delta)
  
  list(precNew=precNew,locNew=locNew,logZnew=logZnew)
}


# ------------------------ GLM exact moment matching ---------------------------
likelihoodGLM <- function(family){
  if(family=="poisson"){
    function(x,y_i,nu_g=NULL){dpois(y_i,exp(x))}
  } else if(family=="gamma"){
    function(x,y_i,nu_g){dgamma(y_i,nu_g,scale=exp(x))}
  } else if(family=="negbin"){
    function(x,y_i,nu_g){dnbinom(y_i,size=nu_g,mu=exp(x))}
  }
}

momentMatchingGLMexact <- function(glmLikelihood,y_i,m_i,s2_i,nu_g=NULL){
  
  # integration limit
  xLim = max(abs(floor(m_i-7*sqrt(s2_i))%/%10),abs(ceiling(m_i+7*sqrt(s2_i))%/%10))*10
  
  # numerical evaluation of hybrid moments
  risZ = Rmpfr::integrateR(function(x) glmLikelihood(x,y_i,nu_g)*dnorm(x,m_i,sqrt(s2_i)),
                           lower=-xLim,upper=xLim,ord = NULL,abs.tol = 1e-4, rel.tol = 1e-4)
  risMean = Rmpfr::integrateR(function(x) x*glmLikelihood(x,y_i,nu_g)*dnorm(x,m_i,sqrt(s2_i)),
                              lower=-xLim,upper=xLim,ord = NULL,abs.tol = 1e-4, rel.tol = 1e-4)
  risVar = Rmpfr::integrateR(function(x) x^2*glmLikelihood(x,y_i,nu_g)*dnorm(x,m_i,sqrt(s2_i)),
                             lower=-xLim,upper=xLim,ord = NULL,abs.tol = 1e-4, rel.tol = 1e-4)
  
  # extracting the moments of interest
  Zexact = risZ$value
  meanExact = risMean$value/Zexact
  varExact = risVar$value/Zexact - meanExact^2
  
  # i-th site updated parameters
  precNew = 1./varExact - 1./s2_i
  locNew = meanExact/varExact - m_i/s2_i
  
  list(precNew=precNew,locNew=locNew,logZnew=log(Zexact))
  
}


# ---------------------- MODEL-SPECIFIC moment matching ------------------------
momentMatching <- function(family,approxHybrid=TRUE){
  
  if(family=="probit"){
    function(y_i,m_i,s2_i,nu_g=NULL){momentMatchingPROBIT(y_i,m_i,s2_i)}
    
  } else if(family=="logit") {
    function(y_i,m_i,s2_i,nu_g=NULL){momentMatchingLOGIT(y_i,m_i,s2_i)}
    
  } else if(family=="poisson") {
    function(y_i,m_i,s2_i,nu_g=NULL){
      if(approxHybrid) {
        momentMatchingPOISSONapprox(y_i,m_i,s2_i)
      } else {
        likelihoodPoiss <- likelihoodGLM("poisson")
        momentMatchingGLMexact(likelihoodPoiss,y_i,m_i,s2_i)
      }
    }
    
  } else if(family=="gamma") {
    if(!is.null(nu)){stop("Invalid input: provide a positive parameter 'nu' for gamma glm")}
    function(y_i,m_i,s2_i,nu_g){
      if(approxHybrid){
        momentMatchingGAMMA(y_i,m_i,s2_i,nu_g)
      } else {
        likelihoodGamma <- likelihoodGLM("gamma")
        momentMatchingGLMexact(likelihoodGamma,y_i,m_i,s2_i,nu_g)
      }
    }
    
  } else if (family=='negbin'){  
    function(y_i,m_i,s2_i,nu_g){
      likelihoodGamma <- likelihoodGLM("negbin")
      momentMatchingGLMexact(likelihoodGamma,y_i,m_i,s2_i,nu_g)
    }
    
  } else {
    stop("Invalid input: 'family' must be on of 'probit', 'logit', 'poisson', 'gamma' or 'negbin'")
  }
}


# ------------------------------ EFFICIENT EP ----------------------------------
getParamsEP = function(X,y,family,Omega0,beta0,nu_g=NULL,
                       tolerance=1e-3,damp=1.,maxIter=1e3,nPrint=100,
                       fullVar=FALSE,predictive=FALSE,approxHybrid=TRUE){
  
  MomMatch <- momentMatching(family,approxHybrid)
  
  if(family=='poisson' | family=='negbin'){
    ySort = order(y,decreasing = TRUE)
    y = y[ySort]
    X = X[ySort,]
  }
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  ### Pre-Computations
  
  if(diag_prior){
    Q0 = diag(1/diag(Omega0),p,p)
    r = beta0/diag(Omega0)
    if(p<n){
      Omega = Omega0
      Q = Q0
      logDetQ0 = sum(log(diag(Q0)))
    }else{
      U = diag(Omega0)*t(X)
      U0 = U
    }
  } else {
    Q0 = solve(Omega0)
    r = Q0%*%beta0
    if(p<n){
      Omega = Omega0
      Q = Q0
      logDetQ0 = determinant(Q0, logarithm = TRUE)
    }else{
      U = tcrossprod(Omega0,X)
      U0 = U
    }
  }
  
  logZ0 = 0.5*crossprod(r,beta0)
  
  ### Initialization
  
  logZ = double(length = n)
  k = double(length = n)
  m = double(length = n)
  
  skipCounts = double(length = maxIter)
  
  diff = Inf
  nIter = 0
  
  ### Iterations
  
  while(diff > tolerance && nIter < maxIter){
    
    countSkips = 0
    
    kOld = k
    mOld = m
    
    if(p<n){
      
      for(i in c(1:n)){
        
        xi = X[i,]
        r_i = r - m[i]*xi
        Q_i = Q - k[i]*tcrossprod(xi)
        
        Oxi = Omega%*%xi
        Oi = Omega + tcrossprod(Oxi)*k[i]/as.double(1.-k[i]*crossprod(xi,Oxi))
        Oixi = Oi%*%xi
        xiOixi = as.double(crossprod(xi,Oixi))
        
        
        r_iOixi = as.double(crossprod(r_i,Oixi))
        
        newCoeff = MomMatch(y[i],r_iOixi,xiOixi,nu_g)
        
        if(newCoeff$precNew>=0){
          
          k[i] = damp*newCoeff$precNew + (1.-damp)*kOld[i]
          m[i] = damp*newCoeff$locNew + (1.-damp)*mOld[i]
          
          logZ[i] = newCoeff$logZnew + 0.5*log(1.+k[i]*xiOixi) +
            0.5*((r_iOixi)^2)/xiOixi - 0.5*((m[i] + r_iOixi/xiOixi)^2)*xiOixi/(1.+k[i]*xiOixi)
          
          r = r_i + m[i]*xi
          Q = Q_i + k[i]*tcrossprod(xi)
          
          Omega = Oi - tcrossprod(Oixi)*k[i]/(1.+k[i]*xiOixi)
          
        }else{countSkips = countSkips+1}
      }
    } else {
      
      for(i in c(1:n)){
        
        u = U[,i]
        xi = X[i,]
        xTu = crossprod(xi,u)[1]
        
        d = 1-k[i]*xTu
        w = u/d
        xTw = xTu/d
        
        r_iTw = crossprod(r,w) - m[i]*xTw
        
        newCoeff = MomMatch(y[i],r_iTw,xTw,nu_g)
        
        if(newCoeff$precNew>=0){
          
          k[i] = damp*newCoeff$precNew + (1.-damp)*kOld[i]
          m[i] = damp*newCoeff$locNew + (1.-damp)*mOld[i]
          
          r = r + (m[i]-mOld[i])*xi
          
          ratio = (kOld[i]-k[i])/(1.+(k[i]-kOld[i])*xTu)
          U = U + (u*ratio)%*%crossprod(xi,U)
          
          logZ[i] = newCoeff$logZnew + 0.5*log(1.+k[i]*xTw) +
            0.5*((r_iTw)^2)/xTw - 0.5*((m[i]+r_iTw/xTw)^2)*xTw/(1.+k[i]*xTw)
          
        }else{countSkips = countSkips+1}
      }
    }
    
    diff = max(abs(c(kOld - k, mOld - m)))
    nIter = nIter + 1
    skipCounts[nIter] = countSkips
    
    if(nIter %% nPrint == 0) {print(paste0("---iteration ",nIter))}
  }
  
  ### Posterior Approximate Moments
  
  if(p<n){
    
    meanBeta = Omega%*%r
    diagOmega = diag(Omega)
    
    logDetQ = determinant(Q, logarithm = TRUE)
    logML = sum(logZ) + 0.5*crossprod(r,meanBeta) - logZ0 + 
      0.5*logDetQ0$modulus[1] - 0.5*logDetQ$modulus[1] 
    
  }else{
    
    diagOmega = diag(Omega0) - rowSums(U*t(k*t(U0)))
    meanBeta = beta0 + U0%*%m - U%*%(k*crossprod(U0,r))
    
    logDet = determinant(diag(1,n,n) + k*X%*%U0, logarithm = TRUE)
    logML = sum(logZ) + 0.5*crossprod(r,meanBeta) - logZ0 - 0.5*logDet$modulus[1]
  }
  
  results = list(meanBeta = meanBeta, diagOmega = diagOmega, logML = logML, 
                 nIter = nIter, kEP = k, mEP = m)
  
  if(fullVar==TRUE){
    if(p>=n) {
      Omega = Omega0 - U%*%(k*t(U0))
    }
    results = append(list(Omega=Omega),results)
  }
  
  if(predictive==TRUE){
    if(p>=n) {
      results = c(results,postPredictive=list(U=U,U0=U0))
    }
  }
  
  skipCounts = skipCounts[1:nIter]
  if(max(skipCounts)>0){
    idx = which(skipCounts > 0)
    print("Iterations:")
    print(paste(as.character(idx), collapse=", "), sep=" ")
    print("Number of skipped units:")
    print(paste(as.character(skipCounts[idx]), collapse=", "), sep=" ")
  }
  
  return(results)
}
