Incidence = read.table(file.choose())

library(pomp)
library(plyr)
library(reshape2)
library(magrittr)
library(ggplot2)
theme_set(theme_bw())
stopifnot(packageVersion("pomp")>="1.4.8")
set.seed(594709947L)

daturl <- "http://kingaa.github.io/pomp/vignettes/twentycities.rda"
datfile <- file.path(tempdir(),"twentycities.rda")
download.file(daturl,destfile=datfile,mode="wb")
load(datfile)
#measles london cases
measles %>% 
  mutate(year=as.integer(format(date,"%Y"))) %>%
  subset(town=="London" & year>=1950 & year<1964) %>%
  mutate(time=(julian(date,origin=as.Date("1950-01-01")))/365.25+1950) %>%
  subset(time>1950 & time<1964, select=c(time,cases)) -> dat
demog %>% subset(town=="London",select=-town) -> demogLondon

dat %>% ggplot(aes(x=time,y=cases))+geom_line()

#####cases with demograpgy.
m1 %>% as.data.frame() %>% 
  melt(id="time") %>%
  ggplot(aes(x=time,y=value))+
  geom_line()+
  facet_grid(variable~.,scales="free_y")