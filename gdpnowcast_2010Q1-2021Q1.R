# R codes to nowcast 2010 Q1 to 2021 Q1 in a more streamlined manner

# Packages required
library(tidyverse);library(DMwR);library(jsonlite);library(tempdisagg);library(tsbox)
library(tsoutliers);library(forecast);library(nowcasting);library(imputeTS)
options(scipen=999)

yoy_d4 <- function(x){
  (x-lag(x,4))/lag(x,4)
}

################################################################################
# 0.Notes:

# 3 eqns were used to predict SG GDP:
# EQ#1: GDP ~ GPI + SPI + D -1
# EQ#2: GPI ~ iip + pp + gpi.L1 + D -1
# EQ#3: SPI ~ fnb + svcexp + D + fnb.L1 + spi.L4 -1
# Typically, nowcasts of GDP are released two weeks after the end of a quarter.
# During which, we will observe two months of iip, pp, fnb in a given quarter.
# We will need to predict the third month, an using ARIMA model.
# Thereafter, we average the 3 months to obtain quarterly iip, pp, fnb.
# The quarterly variables are then used to predict GPI and SPI in EQ#2 and EQ#3.
# Finally, we use the predicted GPI and SPI to predict GDP in EQ#1.
# The results below reproduce GDP nowcasts for 2010 Q1 to 2021 Q1.

# Disclaimer: The results/comments made are for educational purposes only,
# and are solely my own interpretations to nowcast SG GDP.
################################################################################


#######################################################
# 1. First, we will need the GDP key data from Singstat
#######################################################
keydata <- fromJSON("https://www.tablebuilder.singstat.gov.sg/publicfacing/api/json/title/17113.json")

GDP <- keydata[1] %>% as.data.frame() %>%
  select(Level1.quarter,Level1.value) %>%
  rename(
    date = Level1.quarter,
    GDP = Level1.value
  )

GPI <- keydata[2] %>% as.data.frame() %>%
  filter(Level2.level_2 %in% "Goods Producing Industries") %>%
  select(Level2.quarter,Level2.value) %>%
  rename(
    date = Level2.quarter,
    GPI = Level2.value
  )

SPI <- keydata[2] %>% as.data.frame() %>%
  filter(Level2.level_2 %in% "Services Producing Industries") %>%
  select(Level2.quarter,Level2.value) %>%
  rename(
    date = Level2.quarter,
    SPI = Level2.value
  )

mykeydata <- data.frame(
  date = GDP$date,
  GDP = GDP$GDP,
  GPI = GPI$GPI,
  SPI = SPI$SPI,
  D = ifelse(GDP$GDP < 0 , 1, 0)
) %>%
  mutate_at(vars(-date),as.numeric) %>%
  mutate(
    GDP = GDP/100,
    GPI = GPI/100,
    SPI = SPI/100
  )

# LM - EQ#1 (using actual GPI and SPI to predict GDP)
lmmod <- lm(GDP ~ GPI + SPI + D -1, data=mykeydata)
summary(lmmod)
plot(mykeydata$GDP,type="l")
lines(lmmod$fitted.values,col="red")

# use up to 2009 Q4 to train the model and
# use a loop to predict one step ahead recursively, for 2010 Q1 to 2021 Q1:
lmmod <- list()
lmmod.pred <- list()
qtrs <- seq(136,180,by=1)
for (i in qtrs){
  lmmod[[i]] <- lm(GDP ~ GPI + SPI + D -1, data=mykeydata[1:i,])
  
  lmmod.pred[[i]] <- predict(lmmod[[i]], newdata=mykeydata[(i+1),])
}
# do.call
lmmod.pred_df <- do.call(rbind.data.frame,lmmod.pred)
colnames(lmmod.pred_df) <- "pred"

# plot
plot(mykeydata$GDP[137:181],type="o",las=1,pch=16,ylim=c(-0.15,0.2), 
     main="GDP using actual GPI and SPI")
lines(lmmod.pred_df$pred,col="red",las=1,type="o",pch=16)
regr.eval(100*(mykeydata$GDP[137:181]),100*(lmmod.pred_df$pred))

#################################
# 2. GPI equation
#################################

####################
# 2.1 Setup the data
####################
#### setup GPI and D as ts ####
gpi_qtr = ts(mykeydata$GPI, start = c(1976,1), freq=4)

D_qtr = ts(mykeydata$D, start = c(1976,1), freq=4)

#### setup IIP ####
yoy_d12 <- function(x){
  (x-lag(x,12))/lag(x,12)
}

iip <- fromJSON("https://www.tablebuilder.singstat.gov.sg/publicfacing/api/json/title/16863.json")
iip1 <-iip[1] %>% as.data.frame()
iip1 <- iip1 %>% select(Level1.month,Level1.value)
colnames(iip1) <- c("date_str","iip")
iip1$iip <- as.numeric(iip1$iip)

# make iip1 yoy %
iip1 <- iip1 %>%
  mutate_at(vars(-date_str),yoy_d12)


# cast iip as ts
iip_mth <- ts(iip1$iip, start=c(1983,1),freq=12)
iip_mth


# add NA to the end if there are no 3 iip obs in a given qtr
if (length(iip_mth) %% 3 == 2){
  iip_mth <- ts(c(iip_mth,NA),start=c(1983,1),freq=12)
} else if (length(iip_mth) %% 3 == 1){
  iip_mth <- ts(c(iip_mth,NA,NA),start=c(1983,1),freq=12)
} else {
  iip_mth <- iip_mth
}

iip_qtr <- ta(iip_mth,conversion="average",to="quarterly")
iip_qtr


### setup pp ###
pp <- fromJSON("https://www.tablebuilder.singstat.gov.sg/publicfacing/api/json/title/15304.json")
pp1<-pp[1] %>% as.data.frame() %>% 
  filter(Level1.level_1 %in% "Total Public & Private Sector") %>%
  select(Level1.month,Level1.value)
colnames(pp1) <- c("date_str","pp")
pp1$pp <- as.numeric(pp1$pp)

# make pp1 yoy %
pp1 <- pp1 %>%
  mutate_at(vars(-date_str),yoy_d12)


# cast pp as ts
pp_mth <- ts(pp1$pp, start=c(1981,1),freq=12)

# add NA to the end if there are no 3 pp obs in a given qtr
if (length(pp_mth) %% 3 == 2){
  pp_mth <- ts(c(pp_mth,NA),start=c(1981,1),freq=12)
} else if (length(pp_mth) %% 3 == 1){
  pp_mth <- ts(c(pp_mth,NA,NA),start=c(1981,1),freq=12)
} else {
  pp_mth <- pp_mth
}
pp_mth

pp_qtr <- ta(pp_mth,conversion="average",to="quarterly")
pp_qtr

### combine the gpi, iip, pp data
gpidat <- cbind(
  gpi = gpi_qtr, iip = iip_qtr, pp = pp_qtr, 
  pp.L4 = ts( c(rep(NA,4), pp_qtr[1:(NROW(pp_qtr)-4)]), start=c(1981,1),freq=4),
  D_qtr
)

gpidat <- gpidat %>% ts_df() %>% ts_wide()

gpidat <- gpidat %>%
  mutate(
    iip.L1 = c( rep(NA,1), iip[1:(nrow(gpidat)-1)]),
    iip.L4 = c( rep(NA,4), iip[1:(nrow(gpidat)-4)]),
    gpi.L1 = c( rep(NA,1), gpi[1:(nrow(gpidat)-1)]),
    gpi.L4 = c( rep(NA,4), gpi[1:(nrow(gpidat)-4)]),
  )

# set "D" for 2021 Q1 to 1
gpidat[nrow(gpidat),6] <- 1

# Predict the 2021 Q1 value for pp, which we will then subst back into gpidat
pp.arima <- na.omit(pp_qtr)
arima.mod <- tso(pp.arima,types=c("AO","LS","TC"))
arima.mod
arima.mod$fit$xreg
pp.fcst <- forecast(arima.mod$fit,h=1,xreg=cbind(TC154=0.343))
pp.fcst$mean

# Put the point fcst back into gpidat
gpidat[nrow(gpidat),4] <- as.numeric(pp.fcst$mean)

colnames(gpidat) <- c("date", "gpi", "iip", "pp", "pp.L4", "D", "iip.L1",
                      "iip.L4", "gpi.L1", "gpi.L4")


###################################################
# 2.2 use the actual iip and pp data to predict GPI
###################################################

# **Predict gpi one step ahead recursively using actual iip and pp data **
# **This serves as the "ideal" benchmark, if we fully observe iip and pp (no delays) **
# create a loop to forecast GPI one step ahead recursively using iip, pp, pp.L4 and D
gpidat.train <- list()
gpidat.test <- list()
gpi.mod <- list()
gpi.fcst <- list()

# start from 1985 Q1 in gpidat (obs 37) to 2009 Q4 (obs 136)
# predict 2010 Q1 (obs 137) to 2021 Q1 (181)
qtrs <- seq(136,180,by=1)

for (i in qtrs){
  gpidat.train[[i]] <- gpidat %>% slice(37:i)
  
  gpidat.test[[i]] <- gpidat %>% slice(i + 1)
  
  gpi.mod[[i]] <- lm(gpi ~ iip + pp + D  -1
                     ,data=gpidat.train[[i]])
  
  gpi.fcst[[i]] <- predict(gpi.mod[[i]], gpidat.test[[i]])
}

# do.call
gpi.fcst_df <- do.call(rbind.data.frame, gpi.fcst)
colnames(gpi.fcst_df) <- "pred"
gpi.fcst_df <- na.omit(gpi.fcst_df)

plot(gpidat$gpi[137:181],type="l")
lines(gpi.fcst_df$pred,col="red")
regr.eval(gpidat$gpi[137:181],gpi.fcst_df$pred)
# mae    rmse
# 0.009 0.01319


##############################################
# 2.3 Realistic nowcast of GPI
##############################################
# In reality, we only observe 2 months in a given quarter.
# We have to predict the third month, usually using ARIMA models.
# Then, we average the 3 months of iip and pp, to get quarterly iip and pp.
# I have ran ARIMA to predict the third month recursively, and averaged them.
# The results are stored in the excel files below.

#### ** Now we use the predicted values of iip and pp to predict GPI ** ####
library(readxl)
iip.pred <- read_excel("C:\\Users\\CHEONG WEI SI\\Documents\\iip_qdfdf.xlsx")
pp.pred <- read_excel("C:\\Users\\CHEONG WEI SI\\Documents\\pp_qdfdf.xlsx")

mygpi <- data.frame(
  date = seq(as.Date("2010-01-01"),length=nrow(iip.pred),by="quarter"),
  gpi = gpidat$gpi[137:nrow(gpidat)],
  gpi.L1 = gpidat$gpi[136:(nrow(gpidat)-1)],
  gpi.L4 = gpidat$gpi[133:(nrow(gpidat)-4)],
  iip = iip.pred$pred,
  pp = pp.pred$pred,
  pp.L4 = c( gpidat$pp.L4[133:136], pp.pred$pred[1:(NROW(pp.pred$pred)-4)] ),
  D = gpidat$D[137:nrow(gpidat)]
)



# **Predict GPI using Predicted iip and pp values **
# Run the reg model loop to predict GPI for 2010 Q1 to 2021 Q1
mygpi.train <- list()
mygpi.mod <- list()
mygpi.fcst <- list()

# start from 1985 Q1 in gpidat (obs 37) to 2009 Q4 (obs 136). 
# predict 2010 Q1 (obs 137) to 2021 Q1 (181)
qtrs <- seq(136,180,by=1)

for (i in qtrs){
  mygpi.train[[i]] <- gpidat %>% slice(37:i)
  
  mygpi.mod[[i]] <- lm(gpi ~ iip + pp + gpi.L1 + D  -1
                       ,data=mygpi.train[[i]])
  
  mygpi.fcst[[i]] <- predict(mygpi.mod[[i]], mygpi[(i-135),])
}

# do.call
mygpi.fcst_df <- do.call(rbind.data.frame, mygpi.fcst)
colnames(mygpi.fcst_df) <- "pred"
mygpi.fcst_df <- na.omit(mygpi.fcst_df)


plot(gpidat$gpi[137:181],type="o",pch=16,main="GPI")
lines(mygpi.fcst_df$pred,col="red",type="o",pch=16)
regr.eval(100*(gpidat$gpi[137:181]),100*(mygpi.fcst_df$pred))
# mae    rmse
# 1.534  2.263


####################################
# 3. SPI Equation
####################################

################
# 3.1 Setup data
################
#### setup SPI and D as ts ####
spi_qtr = ts(mykeydata$SPI, start = c(1976,1), freq=4)

D_qtr = ts(mykeydata$D, start = c(1976,1), freq=4)


#### setup FNB ####
yoy_d12 <- function(x){
  (x-lag(x,12))/lag(x,12)
}

fnb <- fromJSON("https://www.tablebuilder.singstat.gov.sg/publicfacing/api/json/title/17038.json")
fnb1 <-fnb[1] %>% as.data.frame()
fnb1 <- fnb1 %>% select(Level1.month,Level1.value)
colnames(fnb1) <- c("date_str","fnb")
fnb1$fnb <- as.numeric(fnb1$fnb)

# make fnb1 yoy %
fnb1 <- fnb1 %>%
  mutate_at(vars(-date_str),yoy_d12)

# cast fnb as ts
fnb_mth <- ts(fnb1$fnb, start=c(1985,1),freq=12)
fnb_mth


# add NA to the end if there are no 3 fnb obs in a given qtr
if (length(fnb_mth) %% 3 == 2){
  fnb_mth <- ts(c(fnb_mth,NA),start=c(1985,1),freq=12)
} else if (length(fnb_mth) %% 3 == 1){
  fnb_mth <- ts(c(fnb_mth,NA,NA),start=c(1985,1),freq=12)
} else {
  fnb_mth <- fnb_mth
}

fnb_qtr <- ta(fnb_mth,conversion="average",to="quarterly")
fnb_qtr


#### setup RETAIL ####
yoy_d12 <- function(x){
  (x-lag(x,12))/lag(x,12)
}

retail <- fromJSON("https://www.tablebuilder.singstat.gov.sg/publicfacing/api/json/title/16925.json")
retail1 <-retail[1] %>% as.data.frame()
retail1 <- retail1 %>% 
  filter(Level1.level_1 %in% "Total") %>%
  select(Level1.month,Level1.value)
colnames(retail1) <- c("date_str","retail")
retail1$retail <- as.numeric(retail1$retail)

# make retail1 yoy %
retail1 <- retail1 %>%
  mutate_at(vars(-date_str),yoy_d12)

# cast retail as ts
retail_mth <- ts(retail1$retail, start=c(1985,1),freq=12)
retail_mth


# add NA to the end if there are no 3 retail obs in a given qtr
if (length(retail_mth) %% 3 == 2){
  retail_mth <- ts(c(retail_mth,NA),start=c(1985,1),freq=12)
} else if (length(retail_mth) %% 3 == 1){
  retail_mth <- ts(c(retail_mth,NA,NA),start=c(1985,1),freq=12)
} else {
  retail_mth <- retail_mth
}

retail_qtr <- ta(retail_mth,conversion="average",to="quarterly")
retail_qtr


#### setup SVCEXP (Outlook for SPI) ####
outlook<-fromJSON("https://www.tablebuilder.singstat.gov.sg/publicfacing/api/json/title/15221.json")
outlook1 <-outlook[1] %>% as.data.frame()
outlook1 <- outlook1 %>% select(Level1.quarter,Level1.value) %>%
  rename(
    date_str=Level1.quarter,
    outlook=Level1.value
  ) %>%
  mutate_at(vars(-date_str),as.numeric) %>%
  mutate(
    outlook = outlook/100
  )

# convert outlook1 to a quarterly ts, insert NA for first row, and shift other rows down.
outlook_qtr <- ts( c(NA, outlook1$outlook), start=c(1995,1), freq=4)
outlook_qtr

max(time(outlook_qtr))
max(time(fnb_qtr))

if ( max(time(outlook_qtr)) > max(time(fnb_qtr))){
  # drop the last obs of outlook_qtr to be consistent with fnb_qtr's observations
  outlook_qtr <- ts( outlook_qtr[1:(NROW(outlook_qtr)-1)], start=c(1995,1), freq=4)
} else {}



#### Combine spi, fnb, retail and svcexp ####

spidat <- cbind(
  spi = spi_qtr, fnb = fnb_qtr, retail = retail_qtr, svcexp = outlook_qtr
)

spidat <- diff(spidat)

spidat_df <- spidat %>% ts_df() %>% ts_wide()

spidat_df <- spidat_df %>%
  mutate_at(vars(-time),as.numeric) %>%
  mutate(
    fnb.L1 = c(rep(NA,1),fnb[1:(nrow(spidat)-1)]),
    fnb.L4 = c(rep(NA,4),fnb[1:(nrow(spidat)-4)]),
    spi = spi,
    spi.L1 = c( rep(NA,1), spi[1:(nrow(spidat_df)-1)]),
    spi.L4 = c( rep(NA,4), spi[1:(nrow(spidat_df)-4)]),
    svcexp.L1 = c( rep(NA,1), svcexp[1:(nrow(spidat_df)-1)]),
    svcexp.L4 = c( rep(NA,4), svcexp[1:(nrow(spidat_df)-4)]),
    D = ifelse(spi < 0, 1, 0)
  )


# predict the last obs in fnbs_qtr
myfnb <- diff(fnb_qtr)
myfnb <- na.omit(myfnb)
arima.mod <- auto.arima(myfnb)
arima.pred <- forecast(arima.mod,h=1)
arima.pred$mean

which(colnames(spidat_df)=="fnb") # col 3
spidat_df[nrow(spidat_df),3] <- as.numeric(arima.pred$mean)


# predict the last obs in retail_qtr
myretail <- diff(retail_qtr)
myretail <- na.omit(myretail)
arima.mod <- auto.arima(myretail)
arima.pred <- forecast(arima.mod,h=1)
arima.pred$mean

which(colnames(spidat_df)=="retail") # col 4
spidat_df[nrow(spidat_df),4] <- as.numeric(arima.pred$mean)


#######################################################
# 3.2 use the actual fnb and svcexp data to predict SPI
#######################################################

#### Use the actual fnb, retail and svcexp to predict SPI ####
spidat.train <- list()
spidat.test <- list()
spidat.mod <- list()
spidat.pred <- list()

# obs 135 is 2009 Q4, obs 179 is 2020 Q4.
# we want to predict obs 136: 2010 Q1 to obs 180: 2021 Q1
qtrs <- seq(135,179,by=1) 

for (i in qtrs){
  
  # start from obs 78 (1995 Q3) since it is the row where all variables are non-NA
  spidat.train[[i]] <- spidat_df %>% slice(78:i)
  
  spidat.test[[i]] <- spidat_df %>% slice(i+1)
  
  spidat.mod[[i]] <- lm(spi ~  
                   fnb + svcexp + fnb.L1 + spi.L4 +
                   D  -1
                 ,data=spidat.train[[i]])
  
  # The prediction is in first difference of YoY %, so we have to convert it back later.
  spidat.pred[[i]] <- predict(spidat.mod[[i]],spidat.test[[i]])
}

# do.call
spidat.pred_df <- do.call(rbind.data.frame,spidat.pred)
colnames(spidat.pred_df) <- "pred"
spidat.pred_df <- na.omit(spidat.pred_df)

plot(spidat_df$spi[136:180],type="l",las=1)
lines(spidat.pred_df$pred,col="red",las=1)
regr.eval(spidat_df$spi[136:180],spidat.pred_df$pred)


###################################
# Convert back to actual YoY % form
###################################
ngdp <- fromJSON("https://www.tablebuilder.singstat.gov.sg/publicfacing/api/json/title/17166.json")
svc_actual <- ngdp[2] %>% as.data.frame() %>%
  filter(Level2.level_2 %in% "Services Producing Industries") %>%
  select(Level2.quarter,Level2.value) %>%
  rename(
    date_str = Level2.quarter,
    svcs = Level2.value
  ) %>%
  mutate_at(vars(-date_str),as.numeric) %>%
  filter(date_str >= "1998-Q1")

# convert back to YoY % after our prediction
svc_actual <- svc_actual %>%
  mutate(
    svcs.L1 = c(rep(NA,1),svcs[1:(nrow(svc_actual)-1)]),
    svcs.L5 = c(rep(NA,5),svcs[1:(nrow(svc_actual)-5)]),
    chg = (svcs.L1 - svcs.L5)/svcs.L5
  )

svc_actual <- svc_actual %>% filter(date_str >= "2010-Q1")

svc_actual <- svc_actual %>%
  mutate(
    pred = spidat.pred_df$pred,
    predyoy = spidat.pred_df$pred + chg
  )

plot(mykeydata$SPI[137:181],type="l")
lines(svc_actual$predyoy,col="red")
regr.eval(mykeydata$SPI[137:181], svc_actual$predyoy)
# if we use fnb + svcexp + fnb.L1 + fnb.L4 + D
# mae     rmse
# 0.00857 0.01146

# if we use fnb + svcexp + D
# mae     rmse
# 0.00989 0.01304

# if we use fnb + svcexp + fnb.L1 + spi.L4 + D
# mae     rmse
# 0.00827 0.01116


#####################################
# 3.3 Realistic nowcast of SPI
#####################################

#### ** Now we use the predicted values of fnb to predict SPI ** ####
#####################################################################
# Create a realistic exercise by making only last month of fnb as NA
# Use Arima to predict the fnb which are NA
# then average out to obtain the quarterly fnb
# Finally, use the predicted quarterly fnb to predict services output
#####################################################################
fnb_mdf <- fnb_mth %>% ts_df() %>% ts_wide()
fnb_mdf <- fnb_mdf %>% filter(time >= "1998-01-01")
colnames(fnb_mdf) <- c("time","fnb_mth")

# Create `services_mth` in order to create a dummy variable D
services <- ngdp[2] %>% as.data.frame() %>%
  filter(Level2.level_2 %in% "Services Producing Industries") %>%
  select(Level2.quarter,Level2.value) %>%
  rename(
    date_str = Level2.quarter,
    svcs = Level2.value
  ) %>%
  mutate_at(vars(-date_str),as.numeric) %>%
  mutate_at(vars(-date_str),yoy_d4) %>%
  filter(date_str >= "1998-Q1")

services_qtr <- ts(services$svcs,start=c(1998,1),freq=4)
services_qtr
services_mth <- qtr2month(services_qtr)
if (length(services_mth) %% 3 == 2){
  services_mth <- ts(c(services_mth,NA),start=c(1998,1),freq=12)
} else if (length(services_mth) %% 3 == 1){
  services_mth <- ts(c(services_mth,NA,NA),start=c(1998,1),freq=12)
} else {
  services_mth <- services_mth
}

services_mth <- services_mth %>% ts_df() %>% ts_wide()
services_mth <- services_mth %>% fill(value)
services_mth <- services_mth %>% ts_long() %>% ts_ts()
services_qtr <- ta(services_mth,conversion="average",to="quarterly")


# create fnb_mdf
fnb_mdf <- fnb_mdf %>%
  mutate(
    svcs_mth = as.numeric(services_mth)
  )

fnb_mdf <- fnb_mdf %>%
  mutate(
    D = ifelse(svcs_mth < -0.02, 1, 0)
  )

fnb_mdf <- fnb_mdf %>% select(time, fnb_mth, D)

##########################
# loop to predict fnb_qtr
##########################
fnb_mdf.train <- list()
fnb_mdf.test <- list()
fnb_ts <- list()
fnb_mdf.mod <- list()
fnb_mod.pred <- list()
fnb.fcst <- list()

# obs 146 is 2010 M2. Obs 278 is 2021 M2.
# we observe the first 2 months in a given quarter, and have to predict the 3rd month.
qtrs <- seq(146, 278, by=3)

for (i in qtrs){
  fnb_mdf.train[[i]] <- fnb_mdf %>% slice(1:i)
  
  fnb_mdf.test[[i]] <- fnb_mdf %>% slice((i+1))
  
  fnb_ts[[i]] <- fnb_mdf.train[[i]][,c(1:2)] %>% ts_long() %>% ts_ts()
  
  fnb_mdf.mod[[i]] <- auto.arima(fnb_ts[[i]],xreg=fnb_mdf.train[[i]][,3])
  
  fnb_mod.pred[[i]] <- forecast(fnb_mdf.mod[[i]], h=1, xreg = fnb_mdf.test[[i]][,3])
  
  fnb.fcst[[i]] <- c(fnb_mdf.train[[i]][(nrow(fnb_mdf.train[[i]])-1),2], 
                 fnb_mdf.train[[i]][(nrow(fnb_mdf.train[[i]])),2],
                 as.numeric(fnb_mod.pred[[i]]$mean))
  
  fnb.fcst[[i]] <- sum(fnb.fcst[[i]])/3
}

#do.call
fnb.fcst_df <- do.call(rbind.data.frame, fnb.fcst)
colnames(fnb.fcst_df) <- "pred"

plot(spidat_df$fnb[137:180],type="l")
lines(diff(fnb.fcst_df$pred),col="red")

#####################################################
# Now that we have the quarterly fnb predictions,
# we can replace the actual quarterly transformed fnb 
# with our predictions and see how it performs 
# when it comes to predicting services output
#####################################################

# Since the first diff of fnb.fcst_df is quite similar to the actual quarterly fnb,
# we just replace the spidat_df$fnb with the predicted fnb (which is in first diff btw)

myspi <- spidat_df

myspi[137:180,3] <- diff(fnb.fcst_df$pred)

myspi <- myspi %>%
  mutate(
    fnb.L1 = c( rep(NA,1), fnb[1:(nrow(myspi)-1)]),
    fnb.L4 = c( rep(NA,4), fnb[1:(nrow(myspi)-4)])
  )

# try using a loop to see the overall result (the output is SPI in 1st diff)
myspi.train <- list()
myspi.test <- list()
myspi.mod <- list()
myspi.pred <- list()

# obs 135 is 2009 Q4, obs 179 is 2020 Q4
# we want to predict obs 136: 2010 Q1 and obs 180: 2021 Q1
qtrs <- seq(135,179,by=1)

for (i in qtrs){
  
  # start from obs 78 (1995 Q3), which is the earliest row where all variables are non-NA 
  myspi.train[[i]] <- myspi %>% slice(78:i)
  
  myspi.test[[i]] <- myspi %>% slice(i+1)
  
  myspi.mod[[i]] <- lm(spi ~  
                   fnb + svcexp + fnb.L1 + spi.L4 + 
                   D  -1
                 ,data=myspi.train[[i]])
  
  myspi.pred[[i]] <- predict(myspi.mod[[i]],myspi.test[[i]])
}

# do.call
myspi.pred_df <- do.call(rbind.data.frame,myspi.pred)
colnames(myspi.pred_df) <- "pred"
myspi.pred_df <- na.omit(myspi.pred_df)

# plot
plot(myspi$spi[136:180],type="l")
lines(myspi.pred_df$pred,col="red")
regr.eval(myspi$spi[136:180], myspi.pred_df$pred)
# mae    rmse
# 0.0089 0.01177


#### ** Convert the 1st diff SPI prediction back to YoY % ** ####
svc_pred <- ngdp[2] %>% as.data.frame() %>%
  filter(Level2.level_2 %in% "Services Producing Industries") %>%
  select(Level2.quarter,Level2.value) %>%
  rename(
    date_str = Level2.quarter,
    svcs = Level2.value
  ) %>%
  mutate_at(vars(-date_str),as.numeric) %>%
  filter(date_str >= "1998-Q1")

# convert back to YoY % after our prediction
svc_pred <- svc_pred %>%
  mutate(
    svcs.L1 = c(rep(NA,1),svcs[1:(nrow(svc_pred)-1)]),
    svcs.L5 = c(rep(NA,5),svcs[1:(nrow(svc_pred)-5)]),
    chg = (svcs.L1 - svcs.L5)/svcs.L5
  )

svc_pred <- svc_pred %>% filter(date_str >= "2010-Q1")

svc_pred <- svc_pred %>%
  mutate(
    pred = myspi.pred_df$pred,
    predyoy = myspi.pred_df$pred + chg
  )

# Plot the predicted YoY % SPI against the actual SPI (2010 Q1 to 2021 Q1)
plot( mykeydata$SPI[137:181],type="l")
lines(svc_pred$predyoy,col="red")
regr.eval( 100*(mykeydata$SPI[137:181]), 100*(svc_pred$predyoy))
# mae    rmse
# 0.894  1.177


############################################
# 4. Predict GDP using predicted GPI and SPI
############################################

# Now with our predicted GPI and SPI, we can predict GDP

# `mykeydata` will be our training data

# `testdata` will be a subset of `mykeydata`
# where its 2010 Q1 to 2021 Q1 values are replaced with 
# our predicted GPI and SPI values

testdata <- mykeydata %>% filter(date >= "2010-Q1")

testdata <- testdata %>%
  mutate(
    GPI = mygpi.fcst_df$pred,
    SPI = svc_pred$predyoy
  )


# Run a loop to predict one step ahead recursively for GDP ~ GPI + SPI + D
mydata.list <- list()
mymod.list <- list()
myGDPhat.list <- list()

# Note: obs 136 is 2009 Q4, we want to predict obs 137: 2010 Q1
# Note: obs 180 is 2020 Q4, we want to predict obs 181: 2021 Q1
for (i in c(136:180)){
  mydata.list[[i]] <- mykeydata %>% slice(1:i)
  
  mymod.list[[i]] <- lm(GDP ~ GPI + SPI + D -1, data = mydata.list[[i]] )
  
  myGDPhat.list[[i]] <- predict(mymod.list[[i]], newdata=testdata[(i-135),])
}

# do.call
myGDPhat_df <- do.call(rbind.data.frame,myGDPhat.list)
colnames(myGDPhat_df) <- "GDPhat"

# Put the actual and pred into a comparison dataframe
mycompare_df <- data.frame(
  time=seq(as.Date("2010-03-01"),length=(181-137+1),by="quarter"),
  actual = mykeydata$GDP[137:181],
  pred = myGDPhat_df$GDPhat
)

plot(mycompare_df$actual,las=1,type="o", pch=16,main="myGDP pred",ylim=c(-0.15,0.2))
lines(mycompare_df$pred,col="blue",las=1, type="o", pch=16)

regr.eval(100*mycompare_df$actual, 100*mycompare_df$pred)


##############################
# Others: MAE/RMSE comparisons

# SPI eqn: fnb + svcexp + fnb.L1 + fnb.L4 + D
# mae    rmse
# 1.021  1.280

# SPI eqn: fnb + svcexp + fnb.L1 + D
# mae    rmse
# 0.991  1.272

# SPI eqn: fnb + svcexp + D
# mae    rmse
# 0.996  1.279

# SPI eqn: fnb + svcexp + retail + D
# mae    rmse
# 1.074  1.298

# SPI eqn: fnb + svcexp + spi.L4 + D
# mae    rmse
# 0.879  1.117

# After setting D = 1 for GPI eqn in 2021 Q1
# mae    rmse
# 0.835  1.100


##### Export to excel #####
library(writexl)
write_xlsx(mycompare_df, "gdp_pred_2010Q1to2021Q1.xlsx")
