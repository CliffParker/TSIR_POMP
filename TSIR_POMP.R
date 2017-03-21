# Sat Mar 18 16:03:04 2017 ------------------------------
# Create datasets of reported cases and biths and # later population
# Use predict to expand on births #Might not be so neccesary
#----

library(pomp)
require(locpol)
library(plyr)
library(reshape2)
library(magrittr)
library(ggplot2)
library(foreach)
library(doParallel)

registerDoParallel()
# Sat Mar 18 18:27:08 2017 ------------------------------

"creating City datasets"

stew("TSIR_POMP_SRA_DATA",{
  
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
  
  
  
  SSE2 = function(h){
    ESS = 0
    d <- data.frame(X = NewData2$CIncidence)
    d$Y <- NewData2$CBirths
    sd = sqrt(var(d$X))
    fit <- lm(CBirths ~ CIncidence, data=NewData2)
    lpest1 <- locPolSmootherC(d$X,d$Y , d$X, bw = h * sd, deg = 1, gaussK)
    lpest2 <- locpol(Y ~ X, data = d, bw = h * sd , kernel = gaussK, deg = 1, xeval = d$X)
    ESS = sum(lpest2$residuals^2)
    error = fitted.values(fit) - ( lpest1$beta0 + lpest2$lpFit$Y1 * lpest2$lpFit[,1])
    ESS1 = sum(error^2)
    return(c(ESS/1e+12,ESS1/1e+15))
  }
  
  
  "Loading dataset"
  load("~/GitHub/Measles/twentycities.rda")
  demog$town = factor(demog$town)
  measles$town = factor(measles$town)
  
  
  
  for (names in levels(measles$town)) {
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
        birthrate=predict(smooth.spline(x=year+0.5,y=births),x=time-4)$y
      ) -> covar
    merge(tmp,covar, by = "time")->x
    
    NewData<- as.matrix(x)
    Fdat = Biweekly(NewData)# Biweekly data
    Fdat.d = cbind(Fdat[-(1:8),2],Fdat[1:539,3]) # Adjusting for the  delay caused by maternal immunity
    Fdat1 = Cumulative(Fdat.d)
    NewData1 = cbind(Time=seq(from = 1944, to= 1964.983, length.out = length(Fdat1)),CIncidence=Fdat1[,1],CBirths=Fdat1[,2])
    NewData2 = as.data.frame(NewData1)
    
    h = seq(.05,5,.005)
    y = sapply(h,SSE2)
    h = h[abs(unlist(y[2,]-y[1,]))==min(abs(unlist(y[2,]-y[1,])))]
    
    d <- data.frame(X = NewData2$CIncidence)
    d$Y <- NewData2$CBirths
    plot(d$X,d$Y)
    sd = sqrt(var(d$X))
    lpest1 <- locPolSmootherC(d$X,d$Y , d$X, bw = h * sd , deg = 1, gaussK)
    lpest2 <- locpol(Y ~ X, data = d, bw = h * sd , kernel = gaussK, deg = 1, xeval = d$X)
    
    data.frame(Time = NewData2[["Time"]], beta1=lpest1[["beta1"]],residuals=lpest2[["residuals"]])->result
    
    assign( paste0(names,"_SRA"),result)
    
  }
  
})


###################################################################





###############################################################################################
##########################    DATA      CREATED    ############################################
###############################################################################################














