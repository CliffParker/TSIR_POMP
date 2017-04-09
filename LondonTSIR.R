# Tue Mar 21 16:42:09 2017 ------------------------------

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

"creating City datasets"
"Cases"

"Loading dataset"
load("~/GitHub/Measles/twentycities.rda")
demog$town = factor(demog$town)
measles$town = factor(measles$town)



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


#########################################
# Data compactibility


# Please load the attacted london data, it has the reporting rates as 1/b1
# and also the dynamics of the reconstructed susceptibles as residuals. and needed for the codes below
load("~/GitHub/TSIR_POMP/ZTSIR_POMP_SRA.rda")



names(London_SRA)<- c("time","b1","Z")
Covar = cbind(London_cases,London_SRA)
Covar$Tcases = Covar$cases*Covar$b1 
Covar = Covar[,-c(1,2)]

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



                      //log(lambda) = log(B) + alpha * log(lambda) + (alpha * m )/ lambda + log(S + Z);
                      I = rnbinom_mu( 1/I, lambda);                   
                      ")


#rmeas <- Csnippet("cases = rpois(I/b1);")  
rmeas <- Csnippet("cases = rpois(lambda/b1);")
# the number of cases is the number of susceptibles by the rate of reporting 1/bi
dmeas <- Csnippet("lik = dpois(cases,lambda/b1,give_log);")


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



TSIR %>% 
  simulate(params=coef(fit),nsim=9,as.data.frame=TRUE,include.data=TRUE) %>%
  ggplot(aes(x=time,y=cases,group=sim,color=(sim=="data")))+
  guides(color=FALSE)+
  geom_line()+facet_wrap(~sim,ncol=2)


TSIR %>% 
  simulate(params=coef(fit),nsim=100,as.data.frame=TRUE,include.data=TRUE) %>%
  subset(select=c(time,sim,cases)) %>%
  mutate(data=sim=="data") %>%
  ddply(~time+data,summarize,
        p=c(0.05,0.5,0.95),q=quantile(cases,prob=p,names=FALSE)) %>%
  mutate(p=mapvalues(p,from=c(0.05,0.5,0.95),to=c("lo","med","hi")),
         data=mapvalues(data,from=c(TRUE,FALSE),to=c("data","simulation"))) %>%
  dcast(time+data~p,value.var='q') %>%
  ggplot(aes(x=time,y=med,color=data,fill=data,ymin=lo,ymax=hi))+
  geom_ribbon(alpha=0.2)
# 
# plot(simulate(TSIR))
# simulate(TSIR, as=T)

#' The mean of I(t+1) is lambda(T), from the simulation its evident we are getting the periodic mean, but the problem is with the draws
#' from the negative binomial distribution, the model specifies that I(t+1) is negative binomially distributed with expectation lambda(t+1) and clumping parameter I(t)
#' This basically translate to size = 1/ I(t) and probabilty = 1/(lambda(t+1)*I(t) + 1)
#' but im not able to get the right draws , if you could help id be grateful. Thank you. 
#' 
#' 
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



# c(m            B1            B2            B3            B4 
# 5.553188e+01  2.228755e-02  5.984999e-07  7.471505e-01  1.737106e-03 
# B5            B6            B7            B8            B9 
# 3.998573e-05  1.600438e-05  2.817636e-03  4.932737e+00  7.303878e-05 
# B10           B11           B12           B13           B14 
# 5.912202e-03  7.440274e-05  1.670361e-06  1.188753e-04  2.398089e-02 
# B15           B16           B17           B18           B19 
# 4.231911e-08  1.197301e-04  4.846562e-04  2.303584e-05  2.964398e+02 
# B20           B21           B22           B23           B24 
# 4.685109e+00  1.338773e-06  1.627416e-02  1.600515e-07  3.111910e-03 
# B25           B26         alpha             S            I0 
# 4.291352e-05  1.012378e-03 -4.935183e-02  3.267822e+02  8.227073e-11 
# lambda0        theta0 
# 1.499892e+02  1.530597e+01 )









subset(demog, town == "London")-> L



























































