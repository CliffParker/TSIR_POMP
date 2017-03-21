subset(demog, town == "Bedwellty", select = c(year,pop,births))->dta

subset(measles, town == "Bedwellty", select = c(date,cases))->dta_meas
###########################################################################
loc <- url("http://kingaa.github.io/sbied/intro/parus.csv")
dat <- read.csv(loc)
head(dat)
names(dat) = c("time","cases")
dat$time = dat$time - 20
names(dta) = c("time","Z","b1")
##########################################################################
library(pomp)
TSIR <- pomp(dat,times="time",t0=1940)

stochStep <- Csnippet("
                      theta = rpois(m);
                      lambda = B*(S + Z)*pow(I + theta,alpha);
                      double ll = nearbyint(lambda);
                      I = rnbinom(ll,I);
                      ")
rmeas <- Csnippet("cases = rpois(I);")
dmeas <- Csnippet("lik = dpois(cases,I,give_log);")


pomp(TSIR,
     rprocess=discrete.time.sim(step.fun=stochStep,delta.t=1),
     paramnames=c("m","B","alpha","S"),
     params = c(m=50, B = .2, alpha = .5 , S = 10, I.0 = .6, lambda.0 =50, theta.0 =5),
     statenames=c("I","lambda","theta"),
     rmeasure = rmeas,
     dmeasure = dmeas,
     covar = dta,
     tcovar= "time") -> TSIR



plot(simulate(TSIR))
simulate(TSIR, as=T)
