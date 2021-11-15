
library(plm)
# library("devtools")
# install_version("lfe", version = "2.8-3", repos = "http://cran.us.r-project.org")
library(lfe) # no longer available on CRAN
library(lmtest)
library(dplyr)
library(readr)
source("vcov_plm.R")

######################################################
### COMPPARISON PLM VS FELM
######################################################

df = readr::read_csv("df.csv")

# n = 41
# 9 id groups
# 5 year groups
# 3 cluster groups (gid)
# 2 x-vars (x & u)

## GROUP-FIXED EFFECTS----

plm2a=plm(y~x +u,
          data=df,
          model="within", 
          effect="individual", 
          index=c("id", "year"))

# NO CORRECTION (no HC & no df correction)
summary(plm2a) 
summary(felm(y ~ u + x| id | 0 | 0, df))

# HC CORRECTION
coeftest(plm2a, vcov=(41/(41-9-2))*vcovHC(plm2a, type="HC0", method = "white1"))
coeftest(plm2a, vcovHR(plm2a))
summary(felm(y ~ u + x| id | 0 | 0, df), robust = TRUE)

# GROUP CLUSTER
# clustering can also be done using vcovHC but without Stata-like df adjustment
vcovHC(plm2a, type="HC0", cluster="group", method="arellano") 
# default: cluster = group & method = arellano
vcovG(plm2a, type = "HC0", l = 0, inner = "cluster", cluster="group")
# this is the workhorse behind vcovHC

coeftest(plm2a, vcov=vcovG(plm2a, type = "HC0", l = 0, inner = "cluster", cluster="group") )
coeftest(plm2a, vcov=vcovCL(x=plm2a, cluster=df$id, stata = F))

coeftest(x = plm2a, vcov=vcovCL(x=plm2a, cluster=df$id))
summary(felm(y ~ u + x | id | 0 | id, df))

summary(felm(y ~ u + x + as.character(id) | 0 | 0 | id, df))
sqrt(40/(41-2-1)*9/8*0.1683416^2)

# TIME CLUSTER
coeftest(plm2a, vcov=vcovG(plm2a, type = "HC0", l = 0, inner = "cluster", cluster="time") )
coeftest(plm2a, vcov=vcovCL(x=plm2a, cluster=df$year, stata = F))

coeftest(plm2a, vcov=vcovCL(x=plm2a, cluster=df$year))
summary(felm(y ~ u + x| id | 0 | year, df))

# HIGHER-LEVEL CLUSTERING
coeftest(plm2a, vcov=vcovCL(x=plm2a, cluster=df$gid))
summary(felm(y ~ u + x| id | 0 | gid, df))

# TWOWAY CLUSTERING
coeftest(plm2a, vcov=vcovDC(plm2a, type="HC0"))
coeftest(plm2a, vcov=vcovTC(x=plm2a, df$year, df$id, stata = F))

summary(felm(y ~ u + x| id | 0 | (id+year), df))
coeftest(plm2a, vcov=vcovTC(x=plm2a, df$year, df$id, stata = T))
# felm uses the intersection of both the group and time cluster to generate
# a third cluster variable; if the id-year combos are unique then this is equivalent
# to White's heteroskedasticity-robust correction
# since felm allows non-unique year-id combos (and also allows for multi-way clustering), 
# the third variable method is implemented, whereas plm uses the White method
# the df correction for the substraction term takes the number of clusters for the
# third variable for felm whereas plm uses the standard White HC df correction
vcovCL(x=plm2a, cluster=df$idyear)
vcovHC(plm2a, type="HC0", method = "white1")

###########################

## TIME-FIXED EFFECTS----

plm2b=plm(y~x +u,
          data=df,
          model="within", 
          effect="time", 
          index=c("id", "year"))

# NO CORRECTION (no HC & no df correction)
summary(plm2b)
summary(felm(y ~ u + x| year | 0 | 0, df))

# HC CORRECTION
coeftest(plm2b, vcov=(41/(41-5-2))*vcovHC(plm2b, type="HC0", method = "white1"))
coeftest(plm2b, vcovHR(plm2b))
summary(felm(y ~ u + x| year | 0 | 0, df), robust = TRUE)

# GROUP CLUSTER
coeftest(plm2b, vcov=vcovCL(x=plm2b, cluster=df$id))
summary(felm(y ~ u + x| year | 0 | id, df))

# TIME CLUSTER
coeftest(plm2b, vcov=vcovCL(x=plm2b, cluster=df$year))
summary(felm(y ~ u + x| year | 0 | year, df))

# HIGHER-LEVEL CLUSTERING
coeftest(plm2b, vcov=vcovCL(x=plm2b, cluster=df$gid))
summary(felm(y ~ u + x| year | 0 | gid, df))

# TWOWAY CLUSTERING
coeftest(plm2b, vcov=vcovDC(plm2b, type="HC0"))
coeftest(plm2b, vcov=vcovTC(x=plm2b, df$year, df$id))
summary(felm(y ~ u + x| year | 0 | (id+year), df))

###########################

## TWOWAY FIXED EFFECTS----

plm3=plm(y~x +u,
         data=df,
         model="within", 
         effect="twoway", 
         index=c("id", "year"))

# NO CORRECTION (no HC & no df correction)
summary(plm3) 
summary(felm(y ~ u + x| id+year | 0 | 0, df))

# HC CORRECTION
coeftest(plm3, vcov=(41/(41-14+1-2))*vcovHC(plm3, type="HC0", method = "white1"))
coeftest(plm3, vcovHR(plm3))
summary(felm(y ~ u + x| id+year | 0 | 0, df), robust = TRUE)

# GROUP CLUSTER
coeftest(plm3, vcov=vcovCL(x=plm3, cluster=df$id))
summary(felm(y ~ u + x| id+year | 0 | id, df))

# TIME CLUSTER
coeftest(plm3, vcov=vcovCL(x=plm3, cluster=df$year))
summary(felm(y ~ u + x| id+year | 0 | year, df))

# HIGHER-LEVEL CLUSTERING
coeftest(plm3, vcov=vcovCL(x=plm3, cluster=df$gid, stata = F))

coeftest(plm3, vcov=vcovCL(x=plm3, cluster=df$gid))
summary(felm(y ~ u + x| id+year | 0 | gid, df))

sqrt(40/(41-2-1)*3/2*0.0518084^2)

summary(felm(y ~ u + x + as.character(id) + 
               as.character(year)| 0 | 0 | gid, df))

# TWOWAY CLUSTERING
coeftest(plm3, vcov=vcovDC(plm3, type="HC0"))
coeftest(plm3, vcov=vcovTC(x=plm3, df$year, df$id, stata = F))

coeftest(plm3, vcov=vcovTC(x=plm3, df$year, df$id))
summary(felm(y ~ u + x| id+year | 0 | (id+year), df))
