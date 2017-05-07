library(pomp)
library(plyr)
library(reshape2)
library(magrittr)
library(ggplot2)
library(foreach)
library(doParallel)
require(tsiR)
########### TSIR Data and Covar.
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


"Loading dataset"
load("~/GitHub/Measles/twentycities.rda")
demog$town = factor(demog$town)
measles$town = factor(measles$town)

subset(demog,town=="London")

"TSIR"
for (names in c("London")) {
  tmp <- subset(measles, town == names)
  tmp %>%
    dcast(date~"cases", fun.aggregate = sum) %>%
    mutate(year=as.integer(format(date,"%Y"))) %>%
    subset(year>=1944 & year<1965) %>%
    mutate(time=(julian(date,origin=as.Date("1944-01-01")))/365.25+1944) %>%
    subset(time>1944 & time<1965, select=c(time,cases)) -> tmp

  ##################
  tmp1<- subset(demog, town == names)
  tmp1<-tmp1[,-1]

  tmp1 %>% subset(year>=1944 & year<1964) %>%
    summarize(
      time= tmp$time,

      birthrate=(predict(smooth.spline(x=year,y=births),x=time-(16*7/365.25))$y)/52,

      pop=(predict(smooth.spline(x=year,y=pop),x=time)$y)/2) -> covar

  merge(tmp,covar, by = "time")->x

  # adding new to table
  x<- as.matrix(x)
  Bdat = Biweekly(x)# Biweekly data

  #here
  result = cbind(time=seq(from = 1944 , by = (14/365.25), length.out = 547),cases=Bdat[,2],births=Bdat[,3],pop = Bdat[,4])
  result = as.data.frame(result)



  assign( paste0(names,"_TSIR_Data"),result)

}


#############################################
#TSIR packages
London_TSIR_Data$time -> time
London_TSIR_Data$cases-> cases
London_TSIR_Data$births->births
London_TSIR_Data$pop -> pop


London <- twentymeas[["London"]]


tsirdat <- London_TSIR_Data

parest <- estpars(tsirdat, xreg = "cumcases", IP = 2, seasonality = "standard",
      regtype = "lowess", sigmamax = 5, family = "gaussian",
      link = "identity",  alpha = NULL, sbar = NULL,
      printon = F)


################################################################################
################################################################################
################################################################################
################################################################################



run <- runtsir(data=tsirdat, xreg = "cumcases", IP = 2, nsim = 10,
      sigmamax = 3, alpha = NULL, sbar = NULL, family = "gaussian",
      link = "identity", method = "negbin", inits.fit = FALSE,
      epidemics = "cont",pred = "forward", seasonality = "standard",
      add.noise.sd = 0, mul.noise.sd = 0, printon = F, fit = NULL,
      fittype = NULL)
plotres(run)



#sim <- simulatetsir(tsirdat,parms=parest,inits.fit=T)

sim <-simulatetsir(data=tsirdat, nsim = 100, IP = 2, parms=parest, method = "negbin",
      epidemics = "cont", pred = "forward",
      inits.fit = T, add.noise.sd = 0, mul.noise.sd = 0)

plotres(sim)




London <- twentymeas[["London"]]
## Not run:
plotdata(London)
res <- runtsir(data=London,method='pois',nsim=10, IP=2,inits.fit=FALSE)
plotres(res)





run1 <- run(data=tsirdat, xreg = "cumcases", IP = 2, nsim = 10,
      sigmamax = 5, alpha = NULL, sbar = NULL, family = "gaussian",
      link = "identity", method = "negbin", inits.fit = FALSE,
      epidemics = "cont",pred = "forward", seasonality = "standard",
      add.noise.sd = 0, mul.noise.sd = 0, printon = F, fit = NULL,regtype = "lm",
      fittype = NULL)


plotres(run1)
