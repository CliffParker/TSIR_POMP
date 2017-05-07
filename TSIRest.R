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
  Fdat.d = cbind(Fdat[1:547,2]) # Adjusting for the  delay caused by maternal immunity
  NewData1 = cbind(times=seq(from = 1944, by = 1/26.01786, length.out = 547),cases=Fdat.d[,1])
  NewData2 = as.data.frame(NewData1)

  London_cases<-NewData2

}



# Please load the attacted london data, it has the reporting rates as 1/b1
# and also the dynamics of the reconstructed susceptibles as residuals. and needed for the codes below
load("TSIR_SRA.rda")
names(London_SRA)<- c("time","b1","Z")
Covar = cbind(London_cases,London_SRA)
Covar$Tcases = Covar$cases*Covar$b1
Covar = Covar[,-c(1,2)]
Covar$birth = Fdat.d[,2]

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

                  Talpha = logit(alpha);
                  TSbar = log(Sbar);
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


                    Talpha = expit(alpha);
                    TSbar = exp(Sbar);
                    ")

initlz <- Csnippet("
                   I = Tcases;
                   lambda = Tcases;
                   //S =  Sbar + Z;
                   ")

stochStep <- Csnippet("
                      double B;

                      #define max(x, y) (((x) > (y)) ? (x) : (y))
                      #define min(x, y) (((x) < (y)) ? (x) : (y))

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
                      lambda = min( Sbar + Z,  B*(Sbar + Z)*pow(Tcases + theta,alpha));  // working great lambda was I


                      I = rnbinom_mu( I, lambda);
                      ")



#A overdispersed binomial mesurement is used here.

dmeas <- Csnippet("
                  double psi = 0.116;
                  double mn = lambda/b1;
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
                  double mn = lambda/b1;
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
                  "B22","B23","B24","B25","B26",
                  "alpha","Sbar"),
     initializer=initlz,
     params = c(m=13, B1 = 1.285935e-05, B2 = 1.214533e-05, B3 = 1.464063e-05, B4 = 1.277254e-05, B5 = 1.236623e-05, B6 = 1.181611e-05,
                B7 = 1.059032e-05, B8 = 1.183479e-05, B9 = 1.024478e-05, B10 = 1.099727e-05,
                B11 = 1.065480e-05, B12 = 1.136872e-05, B13 = 9.970371e-06, B14 = 9.655970e-06, B15 = 8.828046e-06, B16 = 8.999024e-06 ,
                B17 = 6.775660e-06, B18 = 6.917613e-06, B19 = 7.245048e-06, B20 = 1.110509e-05,
                B21 = 1.418629e-05, B22 = 1.017632e-05, B23 = 1.180947e-05, B24 = 1.052219e-05, B25 = 1.050035e-05 , B26 = 1.175871e-05,
                alpha = 9.70e-01 , Sbar = 114281.6),
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
                   m=0.03, alpha=0.03, Sbar=0.03,
                   B1 = .10, B2 = .10, B3 = .10, B4 = .10, B5 = .10, B6 = .10, B7 = .10, B8 = .10, B9 = .10, B10 = .10,
                   B11 = .10, B12 = .10, B13 = .10, B14 = .10, B15 = .10, B16 = .10, B17 = .10, B18 = .10, B19 = .10, B20 = .10,
                   B21 = .10, B22 = .10, B23 = .10, B24 = .10, B25 = .10, B26 = .10),
                 transform = T,
                 cooling.type = "geometric", cooling.fraction.50 = .05,
                 tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"))

coef(TSIR)<-coef(firstFit)

TSIR %>%
  simulate(params=coef(firstFit),nsim=5,as.data.frame=TRUE,include.data=TRUE) %>%
  ggplot(aes(x=time,y=cases,group=sim,color=(sim=="data")))+
  guides(color=FALSE)+
  geom_line()+facet_wrap(~sim,ncol=2)+ ggtitle("TSIR Cases Reconstruction")+theme_bw()


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


##########################################################################################'
