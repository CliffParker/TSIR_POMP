library(pomp)
library(plyr)
library(reshape2)
library(magrittr)
library(ggplot2)
library(foreach)
library(doParallel)

registerDoParallel()
set.seed(998468235L,kind="L'Ecuyer")
stopifnot(packageVersion("pomp")>="1.4.8")

Biweekly=function(Data){
  n=nrow(Data)/2
  m=ncol(Data) 
  mat = matrix( ,n,m)
  for (i in 0:n - 1 ){
    mat[i + 1,]= rep(0, m)
    for (j in 1:2) {
      x = (2*i)+j
      mat[i+1,] = c(mat[i+1,]) + c(Data[x,])
    }
  }
  return(mat)
}

Cumulative = function(Data){
  n=nrow(Data)
  m=ncol(Data)
  Dta=matrix( ,n,m)
  Dta[1,] = Data[1,]
  for (i in 2:n){
    Dta[i,] = Dta[i-1,] + Data[i,]
  }
  
  return(Dta)
}

"Loading dataset"
load("twentycities.rda")
demog$town = factor(demog$town)
measles$town = factor(measles$town)

"Creating Bi-weekly incidence data"
for (names in c("London")) {
  tmp <- subset(measles, town == names)
  tmp %>% 
    dcast(date~"cases", fun.aggregate = sum) %>%
    mutate(year=as.integer(format(date,"%Y"))) %>%
    subset(year>=1944 & year<1965) %>%
    mutate(time=(julian(date,origin=as.Date("1944-01-01")))/365.25+1944) %>%
    subset(time>1944 & time<1965, select=c(time,cases)) -> tmp
  
  ##################
  
  NewData<- as.matrix(tmp)
  Fdat = Biweekly(NewData)# Biweekly data
  Fdat.d = cbind(Fdat[1:539,2]) # Adjusting for the  delay caused by maternal immunity
  NewData1 = cbind(times=seq(from = 1944, by = 1/26.01786, length.out = 539),cases=Fdat.d[,1])
  NewData2 = as.data.frame(NewData1)
  
  London_cases<-NewData2
  
}



# Please load the attacted london data, it has the reporting rates as 1/b1
# and also the dynamics of the reconstructed susceptibles as residuals. and needed for the codes below
load("ZTSIR_POMP_SRA.rda")
names(London_SRA)<- c("time","b1","Z")
Covar = cbind(London_cases,London_SRA)
Covar$Tcases = Covar$cases*Covar$b1 
Covar = Covar[,-c(1,2)]

"POMP Parts"
toEst <- Csnippet("
                  Tm = log(m);
                  TB1 = log(B1);
                  TB2 = log(B2);
                  TB3 = log(B3);
                  TB4 = log(B4);
                  TB5 = log(B5);
                  TB6 = log(B6);
                  TB7 = log(B7);
                  TB8 = log(B8);
                  TB9 = log(B9);
                  TB10 = log(B10);
                  TB11 = log(B11);
                  TB12 = log(B12);
                  TB13 = log(B13);
                  TB14 = log(B14);
                  TB15 = log(B15);
                  TB16 = log(B16);
                  TB17 = log(B17);
                  TB18 = log(B18);
                  TB19 = log(B19);
                  TB20 = log(B20);
                  TB21 = log(B21);
                  TB22 = log(B22);
                  TB23 = log(B23);
                  TB24 = log(B24);
                  TB25 = log(B25);
                  TB26 = log(B26);
                  
                  TI0 = log(I0);
                  TS = log(S);
                  ")

fromEst <- Csnippet("
                    Tm = exp(m);
                    TB1 = exp(B1);
                    TB2 = exp(B2);
                    TB3 = exp(B3);
                    TB4 = exp(B4);
                    TB5 = exp(B5);
                    TB6 = exp(B6);
                    TB7 = exp(B7);
                    TB8 = exp(B8);
                    TB9 = exp(B9);
                    TB10 = exp(B10);
                    TB11 = exp(B11);
                    TB12 = exp(B12);
                    TB13 = exp(B13);
                    TB14 = exp(B14);
                    TB15 = exp(B15);
                    TB16 = exp(B16);
                    TB17 = exp(B17);
                    TB18 = exp(B18);
                    TB19 = exp(B19);
                    TB20 = exp(B20);
                    TB21 = exp(B21);
                    TB22 = exp(B22);
                    TB23 = exp(B23);
                    TB24 = exp(B24);
                    TB25 = exp(B25);
                    TB26 = exp(B26);
                    
                    
                    TI0 = exp(I0);
                    TS = exp(S);
                    ")

initlz <- Csnippet("
                   I = nearbyint(I0);
                   lambda = nearbyint(lambda0);
                   theta = nearbyint(theta0);
                   
                   ")

stochStep <- Csnippet("
                      double B;
                      
                      
                      // By-week time  seasonality
                      int  tstar = (t - 1944)*26.01786;
                      
                      
                      if (tstar  % 26 == 1)
                      B = B1;
                      else if (tstar  % 26 == 2)
                      B = B2;
                      else if (tstar  % 26 == 3)
                      B = B3;
                      else if (tstar  % 26 == 4)
                      B = B4;
                      else if (tstar  % 26 == 5)
                      B = B5;
                      else if (tstar  % 26 == 6)
                      B = B6;
                      else if (tstar  % 26 == 7)
                      B = B7;
                      else if (tstar  % 26 == 8)
                      B = B8;
                      else if (tstar  % 26 == 9)
                      B = B9;
                      else if (tstar  % 26 == 10)
                      B = B10;
                      else if (tstar  % 26 == 11)
                      B = B11;
                      else if (tstar  % 26 == 12)
                      B = B12;
                      else if (tstar  % 26 == 13)
                      B = B13;
                      else if (tstar  % 26 == 14)
                      B = B14;
                      else if (tstar  % 26 == 15)
                      B = B15;
                      else if (tstar  % 26 == 16)
                      B = B16;
                      else if (tstar  % 26 == 17)
                      B = B17;
                      else if (tstar  % 26 == 18)
                      B = B18;
                      else if (tstar  % 26 == 19)
                      B = B19;
                      else if (tstar  % 26 == 20)
                      B = B20;
                      else if (tstar  % 26 == 21)
                      B = B21;
                      else if (tstar  % 26 == 22)
                      B = B22;
                      else if (tstar  % 26 == 23)
                      B = B23;
                      else if (tstar  % 26 == 24)
                      B = B24;
                      else if (tstar  % 26 == 25)
                      B = B25;
                      else  
                      B = B26;
                      
                      
                      theta = rpois(m);
                      lambda = B*(S + Z)*pow(Tcases + theta,alpha);  // working great lambda was I
                      
                      // In literature the approximation below was used for estimation, this isnt neccesary in pomp
                      //log(lambda) = log(B) + alpha * log(lambda) + (alpha * m )/ lambda + log(S + Z);

                      I = rnbinom_mu( I, lambda);                   
                      ")



#A overdispersed binomial mesurement is used here.

dmeas <- Csnippet("
                  double psi = 0.116;
                  double mn = I/b1;
                  double v = mn*(1.0-1/b1+psi*psi*mn);
                  double tol = 1.0e-18;
                  if (cases > 0.0) {
                  lik = pnorm(cases+0.5,mn,sqrt(v)+tol,1,0)-pnorm(cases-0.5,mn,sqrt(v)+tol,1,0)+tol;
                  } else {
                  lik = pnorm(cases+0.5,mn,sqrt(v)+tol,1,0)+tol;
                  }
                  ")

rmeas <- Csnippet("
                  double psi = 0.116;
                  double mn = I/b1;
                  double v = mn*(1.0-1/b1+psi*psi*mn);
                  double tol = 1.0e-18;
                  cases = rnorm(mn,sqrt(v)+tol);
                  if (cases > 0.0) {
                  cases = nearbyint(cases);
                  } else {
                  cases = 0.0;
                  }
                  ")



##############################################################################################################################################
pomp(London_cases,times="times",t0=1944,
     rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1/26.01786),
     paramnames=c("m","B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11",
                  "B12","B13","B14","B15","B16","B17","B18","B19","B20","B21",
                  "B22","B23","B24","B25","B26","I0","theta0","lambda0",
                  "alpha","S"),
     initializer=initlz,
     params = c(m=1, B1 = 2.695622e-06, B2 = 2.695622e-06, B3 = 2.695622e-06, B4 = 2.695622e-06, B5 = 2.695622e-06, B6 = 2.695622e-06,
                B7 = 2.695622e-06, B8 = 2.695622e-06, B9 = 2.695622e-06, B10 = 2.695622e-06,
                B11 = 2.695622e-06, B12 = 2.695622e-06, B13 = 2.695622e-06, B14 = 2.695622e-06, B15 = 2.695622e-06, B16 = 2.695622e-06,
                B17 = 2.695622e-06, B18 = 2.695622e-06, B19 = 2.695622e-06, B20 =2.695622e-06,
                B21 = 2.695622e-06, B22 = 2.695622e-06, B23 = 2.695622e-06, B24 = 2.695622e-06, B25 = 2.695622e-06, B26 = 2.695622e-06,
                alpha = 1 , S = 481608, I0 = 10, lambda0 =40, theta0 =1),
     statenames=c("I","lambda","theta"),
     rmeasure = rmeas,
     dmeasure = dmeas,
     toEstimationScale = toEst,
     fromEstimationScale = fromEst,
     covar = Covar,
     tcovar= "time") -> TSIR


plot(simulate(TSIR))
firstFit <- mif2(TSIR, Nmif = 150, Np = 1000, #10 was 10000
                 rw.sd = rw.sd(
                   m=0.03, alpha=0.03, S=0.03,
                   B1 = .10, B2 = .10, B3 = .10, B4 = .10, B5 = .10, B6 = .10, B7 = .10, B8 = .10, B9 = .10, B10 = .10,
                   B11 = .10, B12 = .10, B13 = .10, B14 = .10, B15 = .10, B16 = .10, B17 = .10, B18 = .10, B19 = .10, B20 = .10,
                   B21 = .10, B22 = .10, B23 = .10, B24 = .10, B25 = .10, B26 = .10),
                 transform = T,
                 cooling.type = "geometric", cooling.fraction.50 = .05,
                 tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"))

coef(TSIR)<-coef(firstFit)

TSIR %>%
  simulate(params=coef(firstFit),nsim=9,as.data.frame=TRUE,include.data=TRUE) %>%
  ggplot(aes(x=time,y=cases,group=sim,color=(sim=="data")))+
  guides(color=FALSE)+
  geom_line()+facet_wrap(~sim,ncol=2)


TSIR %>%
  simulate(params=coef(firstFit),nsim=100,as.data.frame=TRUE,include.data=TRUE) %>%
  subset(select=c(time,sim,cases)) %>%
  mutate(data=sim=="data") %>%
  ddply(~time+data,summarize,
        p=c(0.05,0.5,0.95),q=quantile(cases,prob=p,names=FALSE)) %>%
  mutate(p=mapvalues(p,from=c(0.05,0.5,0.95),to=c("lo","med","hi")),
         data=mapvalues(data,from=c(TRUE,FALSE),to=c("data","simulation"))) %>%
  dcast(time+data~p,value.var='q') %>%
  ggplot(aes(x=time,y=med,color=data,fill=data,ymin=lo,ymax=hi))+
  geom_ribbon(alpha=0.2)+ theme_bw() + ylab("cases") + ggtitle("TSIR cases with Data")


############################################################################################################################################

stew(file="MyStochTSIR_London_pomp_results.rda",{
  
  registerDoParallel()
  #####################################################################################################'
  
  #' Parameters to be estimated
  #names<-levels(demog$town)
  name<-names <-c("London")
  for (name in names) {
    
    
    #######################################################################'
    #                            POMP BODY
    #######################################################################'
    #######################################################################'
    #######################################################################
   
    get(paste0(name,"_cases")) %>%
      pomp(times="times",t0=1944,
           rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1/26.01786),
           paramnames=c("m","B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11",
                        "B12","B13","B14","B15","B16","B17","B18","B19","B20","B21",
                        "B22","B23","B24","B25","B26","I0","theta0","lambda0",
                        "alpha","S"),
           initializer=initlz,
           params = c(m=1, B1 = 2.695622e-06, B2 = 2.695622e-06, B3 = 2.695622e-06, B4 = 2.695622e-06, B5 = 2.695622e-06, B6 = 2.695622e-06, 
                      B7 = 2.695622e-06, B8 = 2.695622e-06, B9 = 2.695622e-06, B10 = 2.695622e-06,
                      B11 = 2.695622e-06, B12 = 2.695622e-06, B13 = 2.695622e-06, B14 = 2.695622e-06, B15 = 2.695622e-06, B16 = 2.695622e-06, 
                      B17 = 2.695622e-06, B18 = 2.695622e-06, B19 = 2.695622e-06, B20 =2.695622e-06,
                      B21 = 2.695622e-06, B22 = 2.695622e-06, B23 = 2.695622e-06, B24 = 2.695622e-06, B25 = 2.695622e-06, B26 = 2.695622e-06,
                      alpha = 0.99 , S = 481608, I0 = 10, lambda0 =40, theta0 =1),
           statenames=c("I","lambda","theta"),
           rmeasure = rmeas,
           dmeasure = dmeas,
           toEstimationScale = toEst,
           fromEstimationScale = fromEst,
           covar = Covar,
           tcovar= "time") -> m1
    
    ####################################################'
    ####################################################'
    # m1 %>% 
    #   simulate(params=coef(secondFit),nsim=9,as.data.frame=TRUE,include.data=TRUE) %>%
    #   ggplot(aes(x=time,y=cases,group=sim,color=(sim=="data")))+
    #   guides(color=FALSE)+
    #   geom_line()+facet_wrap(~sim,ncol=2)
    # 
    # 
    # m1 %>% 
    #   simulate(params=coef(firstFit),nsim=100,as.data.frame=TRUE,include.data=TRUE) %>%
    #   subset(select=c(time,sim,cases)) %>%
    #   mutate(data=sim=="data") %>%
    #   ddply(~time+data,summarize,
    #         p=c(0.05,0.5,0.95),q=quantile(cases,prob=p,names=FALSE)) %>%
    #   mutate(p=mapvalues(p,from=c(0.05,0.5,0.95),to=c("lo","med","hi")),
    #          data=mapvalues(data,from=c(TRUE,FALSE),to=c("data","simulation"))) %>%
    #   dcast(time+data~p,value.var='q') %>%
    #   ggplot(aes(x=time,y=med,color=data,fill=data,ymin=lo,ymax=hi))+
    #   geom_ribbon(alpha=0.2)+ theme_bw() + ylab("cases") + ggtitle("TSIR cases with Data")
    # # 
    ####################################################'
    
    #first fit with larger sd's of rw
    firstFit <- mif2(m1, Nmif = 150, Np = 2000, #10 was 10000
                     rw.sd = rw.sd(
                       m=0.03, alpha=0.03, S=0.03,
                       B1 = .10, B2 = .10, B3 = .10, B4 = .10, B5 = .10, B6 = .10, B7 = .10, B8 = .10, B9 = .10, B10 = .10,
                       B11 = .10, B12 = .10, B13 = .10, B14 = .10, B15 = .10, B16 = .10, B17 = .10, B18 = .10, B19 = .10, B20 = .10,
                       B21 = .10, B22 = .10, B23 = .10, B24 = .10, B25 = .10, B26 = .10,
                       I0=0.03, theta0=0.03, lambda0=0.03),
                     transform = T,
                     cooling.type = "geometric", cooling.fraction.50 = .05,
                     tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"))
    secondFit<-continue(firstFit, Nmif = 150, Np = 2000,
                        rw.sd = rw.sd(
                          m=0.5,
                          B1 = .10, B2 = .10, B3 = .10, B4 = .10, B5 = .10, B6 = .10, B7 = .10, B8 = .10, B9 = .10, B10 = .10,
                          B11 = .10, B12 = .10, B13 = .10, B14 = .10, B15 = .10, B16 = .10, B17 = .10, B18 = .10, B19 = .10, B20 = .10,
                          B21 = .10, B22 = .10, B23 = .10, B24 = .10, B25 = .10, B26 = .10,
                          alpha=0.5, S=0.5,
                          I0=0.5, theta0=0.5, lambda0=0.5))

    
     theta<-coef(m1)<-coef(secondFit)
    #Saving model
    assign(paste0(name,"_TSIR_Model"),m1)
    
    
    # Point estimate is in secondFit.
    ###' Confidence intervals
    #parnames <- names(coef(m1))
    parr<-parnames<-c("S")
    
    
    #For loop for Confidence intervals
    for (parr in parnames) {
      estpars <- setdiff(names(theta),c(paste0(parr)))# parameter space to be profiled on
      
      
      
      theta.t <- partrans(m1,theta,"toEstimationScale")
      theta.t.hi <- theta.t.lo <- theta.t
      #parspace
      theta.t.lo[estpars] <- theta.t[estpars]-log(2) #Lower bound for parspace to be profiled on
      theta.t.hi[estpars] <- theta.t[estpars]+log(2) #Upper bound for parspace to be profiled on
      #estspace
      FROM <- theta.t[paste0(parr)]-log(2) #Lower bound for parspace to be profiled on
      TO <- theta.t[paste0(parr)]+log(2) #Upper bound for parspace to be profiled on
      
      
      
      profileDesign(
        assign(paste0(parr),seq(from=FROM ,to=TO ,length=20)),# 2 was 20
        lower=theta.t.lo,upper=theta.t.hi,nprof=40           # 4 was 40
      ) -> pd 
      names(pd)[1]<-paste0(parr)
      
      pd <- as.data.frame(t(partrans(m1,t(pd),"fromEstimationScale")))
      
      
      ########### par_rw.sd
      for (par in names(coef(m1))) {
        assign(paste0(par,"_rw.sd"),0.02)
      }
      assign(paste0(parr,"_rw.sd"),0)
      
      ######################################################################################################'
      ######################################################################################################'
      #########################  starters   created   for maximisation ######################################'
      ######################################################################################################'
      
      ###########################################################################################'
      ######################################################################################################'
      
      dtaf<-data.frame()# to store the info
      
      foreach (p=iter(pd,"row"),
               .combine=rbind,
               .errorhandling="remove",
               .inorder=FALSE,
               .options.mpi=list(chunkSize=1,seed=1598260027L,info=TRUE)
      ) %dopar% {
        tryCatch({
          p <- unlist(p)
          
          tic <- Sys.time()
          
          library(magrittr)
          library(plyr)
          library(reshape2)
          library(pomp)
          
          options(stringsAsFactors=FALSE)
          
          ##
          m1 %>% 
            mif2(start = unlist(p),
                 Nmif = 50, Np = 2000, #10 was 10000
                 rw.sd = rw.sd(
                   m=m_rw.sd, alpha=alpha_rw.sd, S=S_rw.sd,
                   B1 = B1_rw.sd, B2 = B2_rw.sd, B3 = B3_rw.sd, B4 = B4_rw.sd, B5 = B5_rw.sd, B6 = B6_rw.sd, B7 = B7_rw.sd,
                   B8 = B8_rw.sd, B9 = B9_rw.sd, B10 = B10_rw.sd, B11 = B11_rw.sd, B12 = B12_rw.sd,
                   B13 = B13_rw.sd, B14 = B14_rw.sd, B15 = B15_rw.sd, B16 = B16_rw.sd, B17 = B17_rw.sd,
                   B18 = B18_rw.sd, B19 = B19_rw.sd, B20 = B20_rw.sd, B21 = B21_rw.sd, B22 = B22_rw.sd,
                   B23 = B23_rw.sd, B24 = B24_rw.sd, B25 = B25_rw.sd, B26 = B26_rw.sd,
                   I0=I0_rw.sd, theta0=theta0_rw.sd, lambda0=lambda0_rw.sd),
                 transform = T,
                 cooling.type = "geometric", cooling.fraction.50 = .1,
                 tol = 1e-17, max.fail = Inf, verbose = getOption("verbose")) %>%
            mif2() -> mf
          
          
          
          
          ## Runs 10 particle filters to assess Monte Carlo error in likelihood
          ##################################################################################################'
          foreach(i=1:10, 
                  .packages="pomp",
                  .options.multicore=list(set.seed=TRUE)
          ) %dopar% {
            pfilter(mf, Np = 2000)
          } -> pf
          ##################################################################################################'       
          
          
          
          
          ll <- sapply(pf,logLik)
          ll <- logmeanexp(ll, se = TRUE)
          nfail <- sapply(pf,getElement,"nfail")
          
          toc <- Sys.time()
          etime <- toc-tic
          units(etime) <- "hours"
          
          data.frame(as.list(coef(mf)),
                     loglik = ll[1],
                     loglik.se = ll[2],
                     nfail.min = min(nfail),
                     nfail.max = max(nfail),
                     etime = as.numeric(etime))
        }, error=function(e){})
        
      }->dtat
      
      dtaf <- rbind(dtaf,dtat)
      # Table of info
      
      assign(paste0(name,"_",parr,"-profile1"),dtaf)
      
      
      # Confidence interval 95%
      maxloglik <- max(dtaf[["loglik"]])
      cutoff <- maxloglik-qchisq(p=0.95,df=1)/2
      subset(dtaf,loglik>cutoff)-> CI
      assign(paste0(name,parr,"_CI"),CI)
      
    }
    
  }
  
  
  #################################################################################################'
  #
  
})
################################################################################################'
#################################################################################################'






















########################################################################################################################################

fit = mif2(TSIR, Nmif = 100, Np = 10000, #10 was 10000
           rw.sd = rw.sd(m=0),
           transform = T,
           cooling.type = "hyperbolic", cooling.fraction.50 = .05,
           tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"))

firstFit <- mif2(TSIR, Nmif = 50, Np = 2000, #10 was 10000
                 rw.sd = rw.sd(
                   m=0.03, alpha=0.03, S=0.03,
                   B1 = .10, B2 = .10, B3 = .10, B4 = .10, B5 = .10, B6 = .10, B7 = .10, B8 = .10, B9 = .10, B10 = .10,
                   B11 = .10, B12 = .10, B13 = .10, B14 = .10, B15 = .10, B16 = .10, B17 = .10, B18 = .10, B19 = .10, B20 = .10,
                   B21 = .10, B22 = .10, B23 = .10, B24 = .10, B25 = .10, B26 = .10,
                   I0=0.03, theta0=0.03, lambda0=0.03),
                 transform = T,
                 cooling.type = "hyperbolic", cooling.fraction.50 = .05,
                 tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"))
secondFit<-continue(firstFit, Nmif = 50, Np = 10000,
                    rw.sd = rw.sd(
                      m=0.5, 
                      B1 = .10, B2 = .10, B3 = .10, B4 = .10, B5 = .10, B6 = .10, B7 = .10, B8 = .10, B9 = .10, B10 = .10,
                      B11 = .10, B12 = .10, B13 = .10, B14 = .10, B15 = .10, B16 = .10, B17 = .10, B18 = .10, B19 = .10, B20 = .10,
                      B21 = .10, B22 = .10, B23 = .10, B24 = .10, B25 = .10, B26 = .10,
                      alpha=0.5, S=0.5, 
                      I0=0.5, theta0=0.5, lambda0=0.5))


coef(TSIR)<-coef(secondFit)





























































