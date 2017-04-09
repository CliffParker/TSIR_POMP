#' # Tue Mar 21 16:42:09 2017 ------------------------------
#' 
#' library(pomp)
#' library(plyr)
#' library(reshape2)
#' library(magrittr)
#' library(ggplot2)
#' library(foreach)
#' library(doParallel)
#' 
#' registerDoParallel()
#' set.seed(998468235L,kind="L'Ecuyer")
#' stopifnot(packageVersion("pomp")>="1.4.8")
#' 
#' 
#' daturl <- "http://kingaa.github.io/pomp/vignettes/twentycities.rda"
#' datfile <- file.path(tempdir(),"twentycities.rda")
#' download.file(daturl,destfile=datfile,mode="wb")
#' load(datfile)
#' 
#' demog$town = factor(demog$town)
#' measles$town = factor(measles$town)
#' 
#' "creating City datasets"
#' "Cases"
#' for (names in levels(measles$town)) {
#'   tmp <- subset(measles, town == names)
#'   tmp %>% 
#'     dcast(date~"cases", fun.aggregate = sum) %>%
#'     mutate(year=as.integer(format(date,"%Y"))) %>%
#'     subset(year>=1944 & year<1965) %>%
#'     mutate(time=(julian(date,origin=as.Date("1944-01-01")))/365.25+1944) %>%
#'     subset(time>1944 & time<1965, select=c(time,cases)) -> tmp
#'   assign(paste0(names,"_cases"),tmp)
#' }
#' 
#' #########################################
#' # Data compactibility
#' 
#' 
#' # Please load the attacted london data, it has the reporting rates as 1/b1
#' # and also the dynamics of the reconstructed susceptibles as residuals. and needed for the codes below
#' names(London_SRA)<- c("time","b1","Z")
#' 
#' toEst <- Csnippet("
#'                   Tm = log(m);
#'                   TB1 = log(B1);
#'                   TB2 = log(B2);
#'                   TB3 = log(B3);
#'                   TB4 = log(B4);
#'                   TB5 = log(B5);
#'                   TB6 = log(B6);
#'                   TB7 = log(B7);
#'                   TB8 = log(B8);
#'                   TB9 = log(B9);
#'                   TB10 = log(B10);
#'                   TB11 = log(B11);
#'                   TB12 = log(B12);
#'                   TB13 = log(B13);
#'                   TB14 = log(B14);
#'                   TB15 = log(B15);
#'                   TB16 = log(B16);
#'                   TB17 = log(B17);
#'                   TB18 = log(B18);
#'                   TB19 = log(B19);
#'                   TB20 = log(B20);
#'                   TB21 = log(B21);
#'                   TB22 = log(B22);
#'                   TB23 = log(B23);
#'                   TB24 = log(B24);
#'                   TB25 = log(B25);
#'                   TB26 = log(B26);
#' 
#'                   TI0 = log(I0);
#'                   TS = log(S);
#'                   ")
#' 
#' fromEst <- Csnippet("
#'                     Tm = exp(m);
#'                     TB1 = exp(B1);
#'                     TB2 = exp(B2);
#'                     TB3 = exp(B3);
#'                     TB4 = exp(B4);
#'                     TB5 = exp(B5);
#'                     TB6 = exp(B6);
#'                     TB7 = exp(B7);
#'                     TB8 = exp(B8);
#'                     TB9 = exp(B9);
#'                     TB10 = exp(B10);
#'                     TB11 = exp(B11);
#'                     TB12 = exp(B12);
#'                     TB13 = exp(B13);
#'                     TB14 = exp(B14);
#'                     TB15 = exp(B15);
#'                     TB16 = exp(B16);
#'                     TB17 = exp(B17);
#'                     TB18 = exp(B18);
#'                     TB19 = exp(B19);
#'                     TB20 = exp(B20);
#'                     TB21 = exp(B21);
#'                     TB22 = exp(B22);
#'                     TB23 = exp(B23);
#'                     TB24 = exp(B24);
#'                     TB25 = exp(B25);
#'                     TB26 = exp(B26);
#'                    
#' 
#'                     TI0 = exp(I0);
#'                     TS = exp(S);
#'                     ")
#' 
#' initlz <- Csnippet("
#'                    I = nearbyint(I0);
#'                    lambda = nearbyint(lambda0);
#'                    theta = nearbyint(theta0);
#'                   
#'                    ")
#' 
#' stochStep <- Csnippet("
#'                       double B;
#'                       //double P;
#'                      
#' 
#'                       
#'                       // By-week time  seasonality
#'                       int  tstar = (t - 1944)*26.01786;
#' 
#' 
#'                       if (tstar  % 26 == 1)
#'                       B = B1;
#'                       else if (tstar  % 26 == 2)
#'                       B = B2;
#'                       else if (tstar  % 26 == 3)
#'                       B = B3;
#'                       else if (tstar  % 26 == 4)
#'                       B = B4;
#'                       else if (tstar  % 26 == 5)
#'                       B = B5;
#'                       else if (tstar  % 26 == 6)
#'                       B = B6;
#'                       else if (tstar  % 26 == 7)
#'                       B = B7;
#'                       else if (tstar  % 26 == 8)
#'                       B = B8;
#'                       else if (tstar  % 26 == 9)
#'                       B = B9;
#'                       else if (tstar  % 26 == 10)
#'                       B = B10;
#'                       else if (tstar  % 26 == 11)
#'                       B = B11;
#'                       else if (tstar  % 26 == 12)
#'                       B = B12;
#'                       else if (tstar  % 26 == 13)
#'                       B = B13;
#'                       else if (tstar  % 26 == 14)
#'                       B = B14;
#'                       else if (tstar  % 26 == 15)
#'                       B = B15;
#'                       else if (tstar  % 26 == 16)
#'                       B = B16;
#'                       else if (tstar  % 26 == 17)
#'                       B = B17;
#'                       else if (tstar  % 26 == 18)
#'                       B = B18;
#'                       else if (tstar  % 26 == 19)
#'                       B = B19;
#'                       else if (tstar  % 26 == 20)
#'                       B = B20;
#'                       else if (tstar  % 26 == 21)
#'                       B = B21;
#'                       else if (tstar  % 26 == 22)
#'                       B = B22;
#'                       else if (tstar  % 26 == 23)
#'                       B = B23;
#'                       else if (tstar  % 26 == 24)
#'                       B = B24;
#'                       else if (tstar  % 26 == 25)
#'                       B = B25;
#'                       else  
#'                       B = B26;
#'                       
#'                       
#'                       
#'                       theta = rpois(m);
#'                       lambda = B*(S + Z)*pow(I + theta,alpha); // working great
#'                       //P = 1/(I*lambda + 1);                     // the prob 
#'                       I = rnbinom_mu( 1/I, lambda);                   
#'                       ")
#' 
#' 
#' rmeas <- Csnippet("cases = rpois(I/b1);")  
#' # the number of cases is the number of susceptibles by the rate of reporting 1/bi
#' dmeas <- Csnippet("lik = dpois(cases,I/b1,give_log);")
#' 
#' 
#' pomp(London_cases,times="time",t0=1944,
#'      rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1/26.01786),
#'      paramnames=c("m","B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11",
#'                   "B12","B13","B14","B15","B16","B17","B18","B19","B20","B21",
#'                   "B22","B23","B24","B25","B26","I0","theta0","lambda0",
#'                   "alpha","S"),
#'      initializer=initlz,
#'      params = c(m=5, B1 = 1e-3, B2 = 1e-3, B3 = 1e-3, B4 = 1e-3, B5 = 1e-3, B6 = 1e-3, 
#'                 B7 = 1e-3, B8 = 1e-3, B9 = 1e-3, B10 = 1e-3,
#'                 B11 = 1e-3, B12 = 1e-3, B13 = 1e-3, B14 = 1e-3, B15 = 1e-3, B16 = 1e-3, 
#'                 B17 = 1e-3, B18 = 1e-3, B19 = 1e-3, B20 = 1e-3,
#'                 B21 = 1e-3, B22 = 1e-3, B23 = 1e-3, B24 = 1e-3, B25 = 1e-3, B26 = 1e-3,
#'                 alpha = .2 , S = 59558.00, I0 = 10, lambda0 =40, theta0 =5),
#'      statenames=c("I","lambda","theta"),
#'      rmeasure = rmeas,
#'      dmeasure = dmeas,
#'      toEstimationScale = toEst,
#'      fromEstimationScale = fromEst,
#'      covar = London_SRA,
#'      tcovar= "time") -> TSIR
#' 
#' 
#' 
#' 
#' plot(simulate(TSIR))
#' simulate(TSIR, as=T)
#' 
#' #' The mean of I(t+1) is lambda(T), from the simulation its evident we are getting the periodic mean, but the problem is with the draws
#' #' from the negative binomial distribution, the model specifies that I(t+1) is negative binomially distributed with expectation lambda(t+1) and clumping parameter I(t)
#' #' This basically translate to size = 1/ I(t) and probabilty = 1/(lambda(t+1)*I(t) + 1)
#' #' but im not able to get the right draws , if you could help id be grateful. Thank you. 
#' 
#' firstFit <- mif2(TSIR, Nmif = 50, Np = 2000, #10 was 10000
#'                  rw.sd = rw.sd(
#'                    m=0.03, alpha=0.03, S=0.03,
#'                    B1 = .10, B2 = .10, B3 = .10, B4 = .10, B5 = .10, B6 = .10, B7 = .10, B8 = .10, B9 = .10, B10 = .10,
#'                    B11 = .10, B12 = .10, B13 = .10, B14 = .10, B15 = .10, B16 = .10, B17 = .10, B18 = .10, B19 = .10, B20 = .10,
#'                    B21 = .10, B22 = .10, B23 = .10, B24 = .10, B25 = .10, B26 = .10,
#'                    I0=0.03, theta0=0.03, lambda0=0.03),
#'                  transform = T,
#'                  cooling.type = "hyperbolic", cooling.fraction.50 = .05,
#'                  tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"))
#' secondFit<-continue(firstFit, Nmif = 50, Np = 10000,
#'                     rw.sd = rw.sd(
#'                       m=0.5, 
#'                       B1 = .10, B2 = .10, B3 = .10, B4 = .10, B5 = .10, B6 = .10, B7 = .10, B8 = .10, B9 = .10, B10 = .10,
#'                       B11 = .10, B12 = .10, B13 = .10, B14 = .10, B15 = .10, B16 = .10, B17 = .10, B18 = .10, B19 = .10, B20 = .10,
#'                       B21 = .10, B22 = .10, B23 = .10, B24 = .10, B25 = .10, B26 = .10,
#'                       alpha=0.5, S=0.5, 
#'                       I0=0.5, theta0=0.5, lambda0=0.5))
#' 
#' 
#' coef(TSIR)<-coef(secondFit)
#' 

# Tue Mar 28 21:19:40 2017 ------------------------------


#Packages to install
install.packages("pomp")
install.packages("plyr")
install.packages("reshape2")
install.packages("magrittr")
install.packages("ggplot2")
install.packages("foreach")
install.packages("doParallel")


#Load packages

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


daturl <- "http://kingaa.github.io/pomp/vignettes/twentycities.rda"
datfile <- file.path(tempdir(),"twentycities.rda")
download.file(daturl,destfile=datfile,mode="wb")
load(datfile)

demog$town = factor(demog$town)
measles$town = factor(measles$town)

"creating City datasets"
"Cases"
for (names in c("London")) {
  tmp <- subset(measles, town == names)
  tmp %>% 
    dcast(date~"cases", fun.aggregate = sum) %>%
    mutate(year=as.integer(format(date,"%Y"))) %>%
    subset(year>=1944 & year<1965) %>%
    mutate(time=(julian(date,origin=as.Date("1944-01-01")))/365.25+1944) %>%
    subset(time>1944 & time<1965, select=c(time,cases)) -> tmp
  assign(paste0(names,"_cases"),tmp)
}


-----------------------------------------------------------------
# created the Datasets required for this work
###############################################################################################################'
###############################################################################################################'
########################       Creating pomp Body parts      ##################################################'
###############################################################################################################'
###############################################################################################################'
# Csnippet for model structure
load("~/GitHub/TSIR_POMP/ZTSIR_POMP_SRA.rda")
  
names(London_SRA)<- c("time","b1","Z")

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
                      //double P;
                      
                      
                      
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
                      lambda = B*(S + Z)*pow(I + theta,alpha); // working great
                      //P = 1/(I*lambda + 1);                     // the prob 
                      I = rnbinom_mu( 1/I, lambda);                   
                      ")


rmeas <- Csnippet("cases = rpois(I/b1);")  
# the number of cases is the number of susceptibles by the rate of reporting 1/bi
dmeas <- Csnippet("lik = dpois(cases,I/b1,give_log);")

#################################################################################################'
#################################################################################################'
#################################################################################################'
##### Main ######################################################################################'
######################################## Object #################################################'
#################################################################################################'
########################### to ##################################################################'
#################################################################################################'
#################################################################################################'
#################################################################### run ########################'
#################################################################################################'
#################################################################################################'
stew(file="MyStochTSIR_pomp_results.rda",{
  
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
    pomp(times="time",t0=1944,
         rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1/26.01786),
         paramnames=c("m","B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11",
                      "B12","B13","B14","B15","B16","B17","B18","B19","B20","B21",
                      "B22","B23","B24","B25","B26","I0","theta0","lambda0",
                      "alpha","S"),
         initializer=initlz,
         params = c(m=5, B1 = 1e-3, B2 = 1e-3, B3 = 1e-3, B4 = 1e-3, B5 = 1e-3, B6 = 1e-3, 
                    B7 = 1e-3, B8 = 1e-3, B9 = 1e-3, B10 = 1e-3,
                    B11 = 1e-3, B12 = 1e-3, B13 = 1e-3, B14 = 1e-3, B15 = 1e-3, B16 = 1e-3, 
                    B17 = 1e-3, B18 = 1e-3, B19 = 1e-3, B20 = 1e-3,
                    B21 = 1e-3, B22 = 1e-3, B23 = 1e-3, B24 = 1e-3, B25 = 1e-3, B26 = 1e-3,
                    alpha = .2 , S = 59558.00, I0 = 10, lambda0 =40, theta0 =5),
         statenames=c("I","lambda","theta"),
         rmeasure = rmeas,
         dmeasure = dmeas,
         toEstimationScale = toEst,
         fromEstimationScale = fromEst,
         covar=get(paste0(name,"_SRA")),
         tcovar= "time") -> m1
    
 ####################################################'
    ####################################################'
    
    #first fit with larger sd's of rw
    firstFit <- mif2(m1, Nmif = 50, Np = 2000, #10 was 10000
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
    
    
    theta<-coef(m1) <- coef(secondFit)
    #Saving model
    assign(paste0(name,"_Model"),m1)
    
    
    # Point estimate is in secondFit.
    ###' Confidence intervals
    #parnames <- names(coef(m1))
     parr<-parnames<-c("R0")
    
    
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
        assign(paste0(parr),seq(from=FROM ,to=TO ,length=30)),# 2 was 20
        lower=theta.t.lo,upper=theta.t.hi,nprof=400            # 4 was 40
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
                 m=0.03, alpha=0.03, S=0.03,
                 B1 = .10, B2 = .10, B3 = .10, B4 = .10, B5 = .10, B6 = .10, B7 = .10, B8 = .10, B9 = .10, B10 = .10,
                 B11 = .10, B12 = .10, B13 = .10, B14 = .10, B15 = .10, B16 = .10, B17 = .10, B18 = .10, B19 = .10, B20 = .10,
                 B21 = .10, B22 = .10, B23 = .10, B24 = .10, B25 = .10, B26 = .10,
                 I0=0.03, theta0=0.03, lambda0=0.03),
               transform = T,
               cooling.type = "geometric", cooling.fraction.50 = .1,
               tol = 1e-17, max.fail = Inf, verbose = getOption("verbose")) %>%
               mif2() -> mf
          
          

          
          ## Runs 10 particle filters to assess Monte Carlo error in likelihood
          ##################################################################################################'
          foreach(i=1:30, 
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


