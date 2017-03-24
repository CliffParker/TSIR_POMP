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


daturl <- "http://kingaa.github.io/pomp/vignettes/twentycities.rda"
datfile <- file.path(tempdir(),"twentycities.rda")
download.file(daturl,destfile=datfile,mode="wb")
load(datfile)

demog$town = factor(demog$town)
measles$town = factor(measles$town)

"creating City datasets"
for (names in levels(demog$town)) {
  tmp<- subset(demog, town == names)
  tmp<-tmp[,-1]
  tmp %>% subset(year>=1944 & year<1964) %>% 
    summarize(
      time=seq(from=min(year),to=max(year),by=1/12),
      pop=predict(smooth.spline(x=year,y=pop),x=time)$y,
      birthrate=predict(smooth.spline(x=year+0.5,y=births),x=time-4)$y
    ) -> covar
  
  assign( paste0(names,"_covar"),covar)
}


"Cases"
for (names in levels(measles$town)) {
  tmp <- subset(measles, town == names)
  tmp %>% 
    dcast(date~"cases", fun.aggregate = sum) %>%
    mutate(year=as.integer(format(date,"%Y"))) %>%
    subset(year>=1944 & year<1965) %>%
    mutate(time=(julian(date,origin=as.Date("1944-01-01")))/365.25+1944) %>%
    subset(time>1944 & time<1965, select=c(time,cases)) -> tmp
  assign(paste0(names,"_cases"),tmp)
}

#########################################
# Data compactibility
names(London_SRA)<- c("time","b0","Z")

toEst <- Csnippet("
                  Tm = log(m);
                  //TB = log(B);
                  //Talpha = log(alpha);
                  TS = log(S);
                  ")

fromEst <- Csnippet("
                    Tm = exp(m);
                    //TB = exp(B);
                    //Talpha = exp(alpha);
                    TS = exp(S);
                    ")

#TSIR <- pomp(London_cases,times="time",t0=1944)


stochStep <- Csnippet("
                      double B;


                      
                      // By-week time  seasonality
                      t = (t - 1944)*26.01786;
                      
                      if (t == 1 % 26)
                      B = B1;
                      else if (t == 2 % 26)
                      B = B2;
                      else if (t == 3 % 26)
                      B = B3;
                      else if (t == 4 % 26)
                      B = B4;
                      else if (t == 5 % 26)
                      B = B5;
                      else if (t == 6 % 26)
                      B = B6;
                      else if (t == 7 % 26)
                      B = B7;
                      else if (t == 8 % 26)
                      B = B8;
                      else if (t == 9 % 26)
                      B = B9;
                      else if (t == 10 % 26)
                      B = B10;
                      else if (t == 11 % 26)
                      B = B11;
                      else if (t == 12 % 26)
                      B = B12;
                      else if (t == 13 % 26)
                      B = B13;
                      else if (t == 14 % 26)
                      B = B14;
                      else if (t == 15 % 26)
                      B = B15;
                      else if (t == 16 % 26)
                      B = B16;
                      else if (t == 17 % 26)
                      B = B17;
                      else if (t == 18 % 26)
                      B = B18;
                      else if (t == 19 % 26)
                      B = B19;
                      else if (t == 20 % 26)
                      B = B20;
                      else if (t == 21 % 26)
                      B = B21;
                      else if (t == 22 % 26)
                      B = B22;
                      else if (t == 23 % 26)
                      B = B23;
                      else if (t == 24 % 26)
                      B = B24;
                      else if (t == 25 % 26)
                      B = B25;
                      else  
                      B = B26;
                      
                      
                      
                      theta = rpois(m);
                      lambda = B*(S + Z)*pow(I + theta,alpha);
                      I = rnbinom(1/I,1/(lambda*I + 1));
                      ")


rmeas <- Csnippet("cases = rpois(I);")
dmeas <- Csnippet("lik = dpois(cases,I,give_log);")


pomp(London_cases,times="time",t0=1944,
     rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1/26.01786),
     paramnames=c("m","B1","B2","B3","B4","B5","B6","B7","B8","B9","B10","B11",
                  "B12","B13","B14","B15","B16","B17","B18","B19","B20","B21",
                  "B22","B23","B24","B25","B26",
                  "alpha","S"),
     params = c(m=5, B1 = 1e-3, B2 = 1e-3, B3 = 1e-3, B4 = 1e-3, B5 = 1e-3, B6 = 1e-3, 
                B7 = 1e-3, B8 = 1e-3, B9 = 1e-3, B10 = 1e-3,
                B11 = 1e-3, B12 = 1e-3, B13 = 1e-3, B14 = 1e-3, B15 = 1e-3, B16 = 1e-3, 
                B17 = 1e-3, B18 = 1e-3, B19 = 1e-3, B20 = 1e-3,
                B21 = 1e-3, B22 = 1e-3, B23 = 1e-3, B24 = 1e-3, B25 = 1e-3, B26 = 1e-3,
                alpha = .2 , S = 2457117, I.0 = 100, lambda.0 =50, theta.0 =5),
     statenames=c("I","lambda","theta"),
     rmeasure = rmeas,
     dmeasure = dmeas,
     toEstimationScale = toEst,
     fromEstimationScale = fromEst,
     covar = London_SRA,
     tcovar= "time") -> TSIR



plot(simulate(TSIR))
simulate(TSIR, as=T)

firstFit <- mif2(TSIR, Nmif = 50, Np = 2000, #10 was 10000
                 rw.sd = rw.sd(
                   m=0.03, alpha=0.03, S=0.03,
                   B1 = .10, B2 = .10, B3 = .10, B4 = .10, B5 = .10, B6 = .10, B7 = .10, B8 = .10, B9 = .10, B10 = .10,
                   B11 = .10, B12 = .10, B13 = .10, B14 = .10, B15 = .10, B16 = .10, B17 = .10, B18 = .10, B19 = .10, B20 = .10,
                   B21 = .10, B22 = .10, B23 = .10, B24 = .10, B25 = .10, B26 = .10,
                   I.0=0.03, theta.0=0.03, lambda.0=0.03),
                 transform = T,
                 cooling.type = "hyperbolic", cooling.fraction.50 = .05,
                 tol = 1e-17, max.fail = Inf, verbose = getOption("verbose"))
secondFit<-continue(firstFit, Nmif = 50, Np = 10000,
                    rw.sd = rw.sd(
                      m=0.5, B=0.5, alpha=0.5, S=0.5, 
                      I.0=0.5, theta.0=0.5, lambda.0=0.5))


coef(TSIR)<-coef(secondFit)
