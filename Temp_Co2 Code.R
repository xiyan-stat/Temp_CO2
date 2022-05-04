########################################################
#              Temperature & CO2                       #
########################################################

########################################################
## 0. Data Preparation

# loading packages
library(ncdf4)
library(lubridate)
library(forecast)
library(ggplot2)
library(fields)
library(quantreg)
library(abind)
library(maps)

## 0.1 Actual Daily Minimum, Maximum, and Average Land Surface Temperature (.nc file)
# read the file
Tmin = NULL
Tmax = NULL
Tmean = NULL
index = NULL


#Open and Combine all files
files= list.files('HadGHCND_TXTN_acts_1950-2014_15102015.nc/',pattern='*.nc',full.names=TRUE)

for(i in seq_along(files)) {
  nc <- nc_open(files[i])
  
  # Read the whole nc file and read the length of the varying dimension [lon, lat, time]
  tmax <- ncvar_get(nc, "tmax")   
  tmin <- ncvar_get(nc, "tmin") 
  lon <- ncvar_get(nc, "longitude")
  lat <- ncvar_get(nc, "latitude")
  nmax <- dim(tmax)
  nmin <- dim(tmin)
  
  tmin <- ncvar_get(nc,'tmin',start=c(1,1,1),count=c(nmin[1],nmin[2],nmin[3]))
  tmax <- ncvar_get(nc,'tmax',start=c(1,1,1),count=c(nmax[1],nmax[2],nmax[3]))
  tmean <- (tmin + tmax)/2
  t <- ncvar_get(nc = nc, "time") 
  # extract month, year
  m <- as.numeric(format(as.Date(t,"0000-01-01"), format="%m"))
  mm <- sprintf("%02d", as.numeric(m))
  year <- as.numeric(format(as.Date(t,"0000-01-01"), format="%Y"))
  yearm <- as.numeric(paste0(year,mm))
  
  Tmin <- abind(Tmin, tmin, along=3)
  Tmax <- abind(Tmax, tmax, along=3)
  Tmean <- abind(Tmean, tmean, along=3)
  index <- c(index, yearm)
}

save(Tmin, Tmax, Tmean, index, file="Temp_act_all.RData")  #Save the daily temperatues for use

## 0.2 Generate the monthly minimum, maximum, and average temperatues 
# Tmin_month = min(Tmin_daily)
# Tmax_month = max(Tmax_daily)
# Tmean_month = mean(Tmean_daily)

## find how many months
require(data.table)
a <- data.table(x=index)
m_list <- a[ , list( id = list(.I) ) , by = x ]

## generate list of monthly temperature for each location
Tmin_month <- list()
Tmax_month <- list()
Tmean_month <- list()

for (i in 1:nrow(m_list)){
  sublist <- min(m_list[i]$id[[1]]):max(m_list[i]$id[[1]])
  Tmin_month[[i]] <- apply(Tmin[,,sublist], 1:2, min, na.rm=TRUE)
  Tmax_month[[i]] <- apply(Tmax[,,sublist], 1:2, max, na.rm=TRUE)
  Tmean_month[[i]] <- apply(Tmean[,,sublist], 1:2, mean, na.rm=TRUE)
  
  Tmin_month[[i]][is.infinite(Tmin_month[[i]])] <- NA
  Tmax_month[[i]][is.infinite(Tmax_month[[i]])] <- NA
  Tmean_month[[i]][is.infinite(Tmean_month[[i]])] <- NA
}

## The monthly data used in the quantile regressions
# Tmin_month: list of monthly minimum temperature 96 longitude * 73 latitude * 781 months
# Tmax_month: list of monthly maximum temperature 96 longitude * 73 latitude * 781 months
# Tmean_month: list of monthly average temperature 96 longitude * 73 latitude * 781 months
# x: a vector of ordered longitude (by 3.75) from -180->176.25
# y: a vector of ordered latitude (by 2.5) from -90->90
# lon: a vector of original order of longitude
# lat: a vector of original order of latitude
# lon_id: a vector of index used to convert from lon to x
# lat_id: a vector of index used to convert from lat to y
# co2.month: monthly co2 concentration
save(Tmin_month,Tmax_month,Tmean_month,m_list,x,y,lon, lat,lon_id, lat_id,co2.month, file="Data.RData")

########################################################
## 1. Quantilde Regression Model 1: by month

load("~/Data.RData")

Ind <- function(m){
  index <- seq(m, 672, by=12)
  return(index)
}

## Three different quantile regression results for each temp series with co2
res_qr.min <- array(NA, dim=c(96,73, 12))
# 
res_qr.max <- array(NA, dim=c(96,73, 12))
# 
res_qr.mean <- array(NA, dim=c(96,73, 12))
#
se_qr.min <- array(NA, dim=c(96,73, 12))
#
se_qr.max <- array(NA, dim=c(96,73, 12))
#
se_qr.mean <- array(NA, dim=c(96,73, 12))
## Three different Ljung Test rejection # 
lb.min.3 <- array(NA, dim=c(96,73, 12))
lb.min.5 <- array(NA, dim=c(96,73, 12))
#
lb.max.3 <- array(NA, dim=c(96,73, 12))
lb.max.5 <- array(NA, dim=c(96,73, 12))
#
lb.mean.3 <- array(NA, dim=c(96,73, 12))
lb.mean.5 <- array(NA, dim=c(96,73, 12))


## Quantile regression
for (i in 1:length(x)){
  for(j in 1:length(y)){
    temp.m <- ts.temp_m(i,j)
    temp1 <- temp.m$ts.min
    temp2 <- temp.m$ts.max
    temp3 <- temp.m$ts.mean
    
    for(m in 1:12){
      # output: 12 month quantile regression for min, max, mean
      index <- Ind(m)
      
      temp1.m <- temp1[index]
      temp2.m <- temp2[index]
      temp3.m <- temp3[index]
      n <- length(index)
      co2.log_ave.m <- rep(NA, n)
      
      if(m %in% c(1,2,3)){
        for (k in 2:n){
          co2.log_ave.m[k] <- 1/4*(log(co2.month[(index[k]-3)]) + log(co2.month[(index[k]-2)]) + log(co2.month[(index[k]-1)]) + log(co2.month[index[k]]))
        }
      }else{
        for (k in 1:n){
          co2.log_ave.m[k] <- 1/4*(log(co2.month[(index[k]-3)]) + log(co2.month[(index[k]-2)]) + log(co2.month[(index[k]-1)]) + log(co2.month[index[k]]))
        }
      }
      
      
      # change tau
      ## min
      if (!anyNA(temp1.m)){
        # quantile by each month
        if(m %in%  c(1,2,3)){
          fit1 <- rq(temp1.m[2:n] ~ co2.log_ave.m[2:n], tau = 0.05) 
        }else{
          fit1 <- rq(temp1.m[1:n] ~ co2.log_ave.m[1:n], tau = 0.05) 
        }
        
        coef1 <- as.numeric(fit1$coefficients[2])
        res1 <- summary.rq(fit1, se="boot")
        se1 <- as.numeric(res1$coefficients[2,2])
        pvalue1 <- as.numeric(res1$coefficients[2,4])
        
        res_qr.min[i,j,m] <- coef1
        se_qr.min[i,j,m] <- se1
        # ljung Box Test
        #fit1.1 <- arima(fit1$residuals, order=c(1,0,0))
        lb.min.3[i,j,m] <- Box.test(fit1$residuals, lag = 3, type ="Ljung-Box")$p.value
        lb.min.5[i,j,m] <- Box.test(fit1$residuals, lag = 5, type ="Ljung-Box")$p.value
      }
      ## max
      if (!anyNA(temp2.m)){
        # quantile by each month
        if(m %in%  c(1,2,3)){
          fit2 <- rq(temp2.m[2:n] ~ co2.log_ave.m[2:n], tau = 0.95) 
        }else{
          fit2 <- rq(temp2.m[1:n] ~ co2.log_ave.m[1:n], tau = 0.95) 
        }
        # get the s.e.
        coef2 <- as.numeric(fit2$coefficients[2]) 
        res2 <- summary.rq(fit2, se="boot")
        se2 <- as.numeric(res2$coefficients[2,2])
        pvalue2 <- as.numeric(res2$coefficients[2,4])
        
        res_qr.max[i,j,m] <- coef2
        se_qr.max[i,j,m] <- se2
        # ljung Box Test
        #fit2.1 <- arima(fit2$residuals, order=c(1,0,0))
        lb.max.3[i,j,m] <- Box.test(fit2$residuals, lag = 3, type ="Ljung-Box")$p.value
        lb.max.5[i,j,m] <- Box.test(fit2$residuals, lag = 5, type ="Ljung-Box")$p.value
      }
      ## mean
      if (!anyNA(temp3.m)){
        # quantile by each month
        if(m %in%  c(1,2,3)){
          fit3 <- rq(temp3.m[2:n] ~ co2.log_ave.m[2:n], tau = 0.5) 
        }else{
          fit3 <- rq(temp3.m[1:n] ~ co2.log_ave.m[1:n], tau = 0.5) 
        } 
        coef3 <- as.numeric(fit3$coefficients[2])
        res3 <- summary.rq(fit3, se="boot")
        se3 <- as.numeric(res3$coefficients[2,2])
        pvalue3 <- as.numeric(res3$coefficients[2,4])
        
        res_qr.mean[i,j,m] <- coef3
        se_qr.mean[i,j,m] <- se3
        # ljung Box Test
        #fit3.1 <- arima(fit3$residuals, order=c(1,0,0))
        lb.mean.3[i,j,m] <- Box.test(fit3$residuals, lag = 3, type ="Ljung-Box")$p.value
        lb.mean.5[i,j,m] <- Box.test(fit3$residuals, lag = 5, type ="Ljung-Box")$p.value
      }
    }
  }
}

save(res_qr.min,res_qr.max,res_qr.mean,se_qr.min, se_qr.max, se_qr.mean,lb.min.3,lb.max.3,lb.mean.3,lb.min.5,lb.max.5,lb.mean.5,file="MothlyQR.RData")

###########################################################
## 2. Quantilde Regression Model 2: monthly dummy variables
# temp = m1 + m2 + ...+ m12 + m1*co2 + m2*co2 + ....+ m12*co2 + 0

## Three different quantile regression results for each temp series with co2
res_qr.min <- array(NA, dim=c(96,73,12))
# 
res_qr.max <- array(NA, dim=c(96,73,12))
# 
res_qr.mean <- array(NA, dim=c(96,73,12))
#
se_qr.min <- array(NA, dim=c(96,73, 12))
#
se_qr.max <- array(NA, dim=c(96,73, 12))
#
se_qr.mean <- array(NA, dim=c(96,73, 12))
## Three different Ljung Test rejection # 
lb.min.3 <- array(NA, dim=c(96,73, 1))
lb.min.5 <- array(NA, dim=c(96,73, 1))
#
lb.max.3 <- array(NA, dim=c(96,73, 1))
lb.max.5 <- array(NA, dim=c(96,73, 1))
#
lb.mean.3 <- array(NA, dim=c(96,73,1))
lb.mean.5 <- array(NA, dim=c(96,73, 1))

Dum <- function(m){
  res = rep(0,12)
  res[m] = 1
  return(res)
}


## Quantile regression
for (i in 1:length(x)){
  for(j in 1:length(y)){
    # temp
    temp.m <- ts.temp_m(i,j)
    temp1 <- temp.m$ts.min
    temp2 <- temp.m$ts.max
    temp3 <- temp.m$ts.mean
    
    n <- length(temp1)
    
    # co2
    co2.log_ave <- rep(NA, n)
    
    # calculation of past 3 months co2 
    for (k in 4:n){
      co2.log_ave[k] <- 1/4*(log(co2.month[(k-3)]) + log(co2.month[(k-2)]) + log(co2.month[(k-1)]) + log(co2.month[k]))
    }
    
    m1 <- rep(Dum(1), 56)[4:n]
    m2 <- rep(Dum(2), 56)[4:n]
    m3 <- rep(Dum(3), 56)[4:n]
    m4 <- rep(Dum(4), 56)[4:n]
    m5 <- rep(Dum(5), 56)[4:n]
    m6 <- rep(Dum(6), 56)[4:n]
    m7 <- rep(Dum(7), 56)[4:n]
    m8 <- rep(Dum(8), 56)[4:n]
    m9 <- rep(Dum(9), 56)[4:n]
    m10 <- rep(Dum(10), 56)[4:n]
    m11 <- rep(Dum(11), 56)[4:n]
    m12 <- rep(Dum(12), 56)[4:n]
    
    llogco2 <- co2.log_ave[4:n]
    
    if (!anyNA(temp1)){
      # quantile with month dummy (no intercept)
      fit1 <- rq(temp1[4:n] ~ m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12+m2*llogco2+m3*llogco2+
                   m4*llogco2+m5*llogco2+m6*llogco2+m7*llogco2+m8*llogco2+m9*llogco2+m10*llogco2+m11*llogco2 +m12*llogco2 , tau = 0.05)
      
      res1 <- summary.rq(fit1, se="boot")
      #pvalue1 <- as.numeric(res1$coefficients[2,4])
      
      res_qr.min[i,j,1] <- as.numeric(fit1$coefficients[13])
      se_qr.min[i,j,1] <- as.numeric(res1$coefficients[13,2])
      
      for(m in 2:12){
        res_qr.min[i,j,m] <- as.numeric(fit1$coefficients[(12+m)]) + as.numeric(fit1$coefficients[13])
        # check the coeffcient of two models
        se_qr.min[i,j,m] <- as.numeric(res1$coefficients[(12+m),2])
      }
      # ljung Box Test
      #fit1.1 <- arima(fit1$residuals, order=c(1,0,0))
      lb.min.3[i,j,1] <- Box.test(fit1$residuals, lag = 3, type ="Ljung-Box")$p.value
      lb.min.5[i,j,1] <- Box.test(fit1$residuals, lag = 5, type ="Ljung-Box")$p.value
      
    }
    ## max
    if (!anyNA(temp2)){
      # quantile with month dummy (no intercept)
      fit2 <- rq(temp2[4:n] ~ m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12+m2*llogco2+m3*llogco2+
                   m4*llogco2+m5*llogco2+m6*llogco2+m7*llogco2+m8*llogco2+m9*llogco2+m10*llogco2+m11*llogco2 +m12*llogco2 , tau = 0.95)
      
      res2 <- summary.rq(fit2, se="boot")
      #pvalue1 <- as.numeric(res1$coefficients[2,4])
      
      res_qr.max[i,j,1] <- as.numeric(fit2$coefficients[13])
      se_qr.max[i,j,1] <- as.numeric(res2$coefficients[13,2])
      
      for(m in 2:12){
        res_qr.max[i,j,m] <- as.numeric(fit2$coefficients[(12+m)]) + as.numeric(fit2$coefficients[13])
        se_qr.max[i,j,m] <- as.numeric(res2$coefficients[(12+m),2])
      }
      # ljung Box Test
      #fit1.1 <- arima(fit1$residuals, order=c(1,0,0))
      lb.max.3[i,j,1] <- Box.test(fit2$residuals, lag = 3, type ="Ljung-Box")$p.value
      lb.max.5[i,j,1] <- Box.test(fit2$residuals, lag = 5, type ="Ljung-Box")$p.value
    }
    ## mean
    if (!anyNA(temp3)){
      # quantile with month dummy (no intercept)
      # quantile with month dummy (no intercept)
      fit3 <- rq(temp3[4:n] ~ m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12+m2*llogco2+m3*llogco2+
                   m4*llogco2+m5*llogco2+m6*llogco2+m7*llogco2+m8*llogco2+m9*llogco2+m10*llogco2+m11*llogco2 +m12*llogco2 , tau = 0.5)
      
      res3 <- summary.rq(fit3, se="boot")
      #pvalue1 <- as.numeric(res1$coefficients[2,4])
      
      res_qr.mean[i,j,1] <- as.numeric(fit3$coefficients[13])
      se_qr.mean[i,j,1] <- as.numeric(res3$coefficients[13,2])
      
      for(m in 2:12){
        res_qr.mean[i,j,m] <- as.numeric(fit3$coefficients[(12+m)]) + as.numeric(fit3$coefficients[13])
        se_qr.mean[i,j,m] <- as.numeric(res3$coefficients[(12+m),2])
      }
      # ljung Box Test
      #fit1.1 <- arima(fit1$residuals, order=c(1,0,0))
      lb.mean.3[i,j,1] <- Box.test(fit3$residuals, lag = 3, type ="Ljung-Box")$p.value
      lb.mean.5[i,j,1] <- Box.test(fit3$residuals, lag = 5, type ="Ljung-Box")$p.value
    }
  }
}

save(res_qr.min,res_qr.max,res_qr.mean,se_qr.min, se_qr.max, se_qr.mean,lb.min.3,lb.max.3,lb.mean.3,lb.min.5,lb.max.5,lb.mean.5,file="MothlyQR_1.RData")

## Rename the results from model 2 (be different as model 1)
# res_qr.min1 <- res_qr.min
# res_qr.max1 <- res_qr.max
# res_qr.mean1 <- res_qr.mean
# se_qr.min1 <- se_qr.min
# se_qr.max1 <- se_qr.max
# se_qr.mean1 <- se_qr.mean
# lb.min.3_1 <- lb.min.3
# lb.max.3_1 <- lb.max.3
# lb.mean.3_1 <- lb.mean.3
# lb.min.5_1 <- lb.min.5
# lb.max.5_1 <- lb.max.5
# lb.mean.5_1 <- lb.mean.5

#############################################################
## 3. Generate the monthly time series for north america: 
# Take model 1 regression results as an example
## 3.1 plot only the north american regions
load("MonthlyQR.RData")
#load("MonthlyQR_1.RData")

# col functions

my.cols <- function(dat, nlevels, col = c("blue", "white", "red"), na.color="#D3D3D3"){
  n = nlevels
  z = c(dat)
  zlim = range(z,na.rm=TRUE)
  n1 = floor(n * (abs(zlim[1])) / (zlim[2] - zlim[1]))
  n2 = n - n1
  col1 = colorRampPalette(col[1:2])(n1)
  col2 = colorRampPalette(col[2:3])(n2)
  cols1 = c(col1, col2)
  
  zstep <- (zlim[2] - zlim[1]) / length(cols1); # step in the color palette
  newz.na <- zlim[2] + 2 * zstep # new z for NA
  
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na
  
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na
  
  cols <- c(cols1, na.color) # we construct the new color range by including: na.color
  
  tmp1 <- ifelse(is.na(dat), newz.na, dat)
  
  return(list(cols=cols, zlim=zlim, data = tmp1))
}

tmp <- res_qr.min
res <- my.cols(dat=tmp, nlevels=400, col = c("blue", "white", "red"), na.color="#D3D3D3")
dat <- res$data
lon_id_n <- lon_id[15:36]
lat_id_n <- lat_id[45:57]
image.plot(x[15:36], y[45:57],
           dat[lon_id_n,lat_id_n,1], las = 1, main = paste0("Subregion: North America"),  
           col = res$cols, zlim = res$zlim,
           xlab = "lon", ylab = "lat")
map("world", add = T)

library(animation)

saveLatex({
  for (i in 1:12){
    tmp <- res_qr.min
    res1 <- my.cols(dat=tmp, nlevels=400, col = c("blue", "white", "red"), na.color="#D3D3D3")
    dat <- res1$data
    lon_id_n <- lon_id[15:36]
    lat_id_n <- lat_id[45:57]
    image.plot(x[15:36], y[45:57],
               dat[lon_id_n,lat_id_n,i], las = 1, main = paste0("5th Quantile regression: Method 1, coef, Tmin, Month=", i),  
               col = res1$cols, zlim = res1$zlim,
               xlab = "lon", ylab = "lat")
    map("world", add = T)
  }
}, img.name = "QR_min_coef_us_1", ani.opts = "controls,width=0.95\\textwidth",
latex.filename = ifelse(interactive(), "QR_min_coef_us.tex", ""), interval = 0.5,
nmax = 12, ani.dev = "pdf", ani.type = "pdf", ani.width = 10,
ani.height = 8,documentclass = paste("\\documentclass{article}",
                                     "\\usepackage[papersize={10in,8in},margin=0.05in]{geometry}",
                                     sep = "\n"))

load("~/Desktop/Varying coefficient models/Daily gridded temperature/Code/MothlyQR_1.RData")
library(animation)

saveLatex({
  for (i in 1:12){
    tmp <- res_qr.min
    res1 <- my.cols(dat=tmp, nlevels=400, col = c("blue", "white", "red"), na.color="#D3D3D3")
    dat <- res1$data
    lon_id_n <- lon_id[15:36]
    lat_id_n <- lat_id[45:57]
    image.plot(x[15:36], y[45:57],
               dat[lon_id_n,lat_id_n,i], las = 1, main = paste0("5th Quantile regression: Method 2, coef, Tmin, Month=", i),  
               col = res1$cols, zlim = res1$zlim,
               xlab = "lon", ylab = "lat")
    map("world", add = T)
  }
}, img.name = "QR_min_coef_us_2", ani.opts = "controls,width=0.95\\textwidth",
latex.filename = ifelse(interactive(), "QR_min_coef_us2.tex", ""), interval = 0.5,
nmax = 12, ani.dev = "pdf", ani.type = "pdf", ani.width = 10,
ani.height = 8,documentclass = paste("\\documentclass{article}",
                                     "\\usepackage[papersize={10in,8in},margin=0.05in]{geometry}",
                                     sep = "\n"))

library(animation)

my.cols1 <- function(dat, nlevels, col = c("blue", "white", "red"), na.color="#D3D3D3"){
  n = nlevels
  z = c(dat)
  zlim = range(z,na.rm=TRUE)
  n1 = floor(n * (abs(zlim[1])<1) / (zlim[2] - zlim[1]))
  n2 = n - n1
  col1 = colorRampPalette(col[1:2])(n1)
  col2 = colorRampPalette(col[2:3])(n2)
  cols1 = c(col1, col2)
  
  zstep <- (zlim[2] - zlim[1]) / length(cols1); # step in the color palette
  newz.na <- zlim[2] + 2 * zstep # new z for NA
  
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na
  
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na
  
  cols <- c(cols1, na.color) # we construct the new color range by including: na.color
  
  tmp1 <- ifelse(is.na(dat), newz.na, dat)
  
  return(list(cols=cols, zlim=zlim, data = tmp1))
}

saveLatex({
  for (i in 1:12){
    tmp <- se_qr.min/se_qr.min_1  # se1/se2 se1: model 1, se2: model 2
    res1 <- my.cols1(dat=tmp, nlevels=400, col = c("blue", "white", "red"), na.color="#D3D3D3")
    dat <- res1$data
    lon_id_n <- lon_id[15:36]
    lat_id_n <- lat_id[45:57]
    image.plot(x[15:36], y[45:57],
               dat[lon_id_n,lat_id_n,i], las = 1, main = paste0("5th Quantile regression: se1/se2, Tmin, Month=", i),  
               col = res1$cols, zlim = res1$zlim,
               xlab = "lon", ylab = "lat")
    map("world", add = T)
  }
}, img.name = "QR_min_se1_se2", ani.opts = "controls,width=0.95\\textwidth",
latex.filename = ifelse(interactive(), "QR_min_se12.tex", ""), interval = 0.5,
nmax = 12, ani.dev = "pdf", ani.type = "pdf", ani.width = 10,
ani.height = 8,documentclass = paste("\\documentclass{article}",
                                     "\\usepackage[papersize={10in,8in},margin=0.05in]{geometry}",
                                     sep = "\n"))


##3.2 functions used

## 3.2.1. location transfer function
loc.trans <- function(x0, y0){
  i = NA
  j = NA
  lon0 = x[x0]  # i,x
  lat0 = y[y0] # j, y
  
  if(x0 <= 48){
    i = x0 + 48
  }else{ i = x0 - 48}
  
  j = which(lat == y[y0])
  
  return(list(i=i,j=j,lon0=lon0, lat0=lat0))
}

## 3.2.2. temperature data transfer function
ts.temp_m <- function(lon_ind, lat_ind){
  ts.min <- NULL
  ts.max <- NULL
  ts.mean <- NULL
  
  
  for (t in 109:780){
    tmp1 <- Tmin_month[[t]][lon_ind, lat_ind]
    tmp2 <- Tmax_month[[t]][lon_ind, lat_ind]
    tmp3 <- Tmean_month[[t]][lon_ind, lat_ind]
    
    ts.min <- c(ts.min, tmp1)
    ts.max <- c(ts.max, tmp2)
    ts.mean <- c(ts.mean, tmp3)
  }
  
  return(structure(list(ts.min=ts.min, ts.max=ts.max, ts.mean=ts.mean)))
}

### 3.3.3. Generate time series for each location (index) using Model I: quantile regression by month
ts.min <- function(index){
  int_min <- c()
  coef_min <- c()
  lb_min <- c()
  ub_min <- c()
  se_min <- c()
  
  index = as.numeric(index)
  ind = loc.trans(index[1], index[2])
  # model 1
  temp.m <- ts.temp_m(ind$i,ind$j)
  temp1 <- temp.m$ts.min
  
  for(m in 1:12){
    # output: 12 month quantile regression for min, max, mean
    index <- Ind(m)
    
    temp1.m <- temp1[index]
    temp2.m <- temp2[index]
    temp3.m <- temp3[index]
    n <- length(index)
    co2.log_ave.m <- rep(NA, n)
    
    if(m %in% c(1,2,3)){
      for (k in 2:n){
        co2.log_ave.m[k] <- 1/4*(log(co2.month[(index[k]-3)]) + log(co2.month[(index[k]-2)]) + log(co2.month[(index[k]-1)]) + log(co2.month[index[k]]))
      }
    }else{
      for (k in 1:n){
        co2.log_ave.m[k] <- 1/4*(log(co2.month[(index[k]-3)]) + log(co2.month[(index[k]-2)]) + log(co2.month[(index[k]-1)]) + log(co2.month[index[k]]))
      }
    }
    
    
    # change tau
    ## min
    if (!anyNA(temp1.m)){
      # quantile by each month
      if(m %in%  c(1,2,3)){
        fit1 <- rq(temp1.m[2:n] ~ co2.log_ave.m[2:n], tau = 0.05) 
        QR.b <- boot.rq(cbind(1,co2.log_ave.m[2:n]),temp1.m[2:n],tau=0.05, R=10000)
      }else{
        fit1 <- rq(temp1.m[1:n] ~ co2.log_ave.m[1:n], tau = 0.05) 
        QR.b <- boot.rq(cbind(1,co2.log_ave.m[1:n]),temp1.m[1:n],tau=0.05, R=10000)
      }
      
      coef0 <- as.numeric(fit1$coefficients[1])
      coef1 <- as.numeric(fit1$coefficients[2])
      # bootstrap se
      res1 <- summary.rq(fit1, se="boot")
      se1 <- as.numeric(res1$coefficients[2,2])
      pvalue1 <- as.numeric(res1$coefficients[2,4])
      # confidence interval 95%
      conf1 <- t(apply(QR.b$B, 2, quantile, c(0.025,0.975)))
      lb1 <- conf1[2,1]
      ub1 <- conf1[2,2]
      
      int_min[m] <- coef0
      coef_min[m] <- coef1
      se_min[m] <- se1
      lb_min[m] <- lb1
      ub_min[m] <- ub1
    }
  }
  
  return(list(int_min=int_min, coef_min=coef_min, se_min=se_min, lb_min=lb_min, ub_min=ub_min))
}

### 3.3.4. Generate time series for each location (index) using Model II: quantile regression with monthly dummy variables
ts.min1 <- function(index){
  int_min <- c()
  coef_min <- c()
  lb_min <- c()
  ub_min <- c()
  se_min <- c()
  
  index = as.numeric(index)
  ind = loc.trans(index[1], index[2])
  # model 1
  temp.m <- ts.temp_m(ind$i,ind$j)
  temp1 <- temp.m$ts.min
  
  n <- length(temp1)
  
  # co2
  co2.log_ave <- rep(NA, n)
  
  # calculation of past 3 months co2 
  for (k in 4:n){
    co2.log_ave[k] <- 1/4*(log(co2.month[(k-3)]) + log(co2.month[(k-2)]) + log(co2.month[(k-1)]) + log(co2.month[k]))
  }
  
  m1 <- rep(Dum(1), 56)[4:n]
  m2 <- rep(Dum(2), 56)[4:n]
  m3 <- rep(Dum(3), 56)[4:n]
  m4 <- rep(Dum(4), 56)[4:n]
  m5 <- rep(Dum(5), 56)[4:n]
  m6 <- rep(Dum(6), 56)[4:n]
  m7 <- rep(Dum(7), 56)[4:n]
  m8 <- rep(Dum(8), 56)[4:n]
  m9 <- rep(Dum(9), 56)[4:n]
  m10 <- rep(Dum(10), 56)[4:n]
  m11 <- rep(Dum(11), 56)[4:n]
  m12 <- rep(Dum(12), 56)[4:n]
  
  llogco2 <- co2.log_ave[4:n]
  
  if (!anyNA(temp1)){
    # quantile with month dummy (no intercept)
    fit1 <- rq(temp1[4:n] ~ m2+m3+m4+m5+m6+m7+m8+m9+m10+m11+m12+m2*llogco2+m3*llogco2+
                 m4*llogco2+m5*llogco2+m6*llogco2+m7*llogco2+m8*llogco2+m9*llogco2+m10*llogco2+m11*llogco2 +m12*llogco2 , tau = 0.05)
    QR.b <- boot.rq(cbind(1,co2.log_ave.m[2:n]),temp1.m[2:n],tau=0.05, R=10000)
    
    res1[1] <- summary.rq(fit1, se="boot")
    #pvalue1 <- as.numeric(res1$coefficients[2,4])
    
    coef1 <- as.numeric(fit1$coefficients[13])
    se1 <- as.numeric(res1$coefficients[13,2])
    
    for(m in 2:12){
      res_qr.min[i,j,m] <- as.numeric(fit1$coefficients[(12+m)]) + as.numeric(fit1$coefficients[13])
      # check the coeffcient of two models
      se_qr.min[i,j,m] <- as.numeric(res1$coefficients[(12+m),2])
    }
    # ljung Box Test
    #fit1.1 <- arima(fit1$residuals, order=c(1,0,0))
    lb.min.3[i,j,1] <- Box.test(fit1$residuals, lag = 3, type ="Ljung-Box")$p.value
    lb.min.5[i,j,1] <- Box.test(fit1$residuals, lag = 5, type ="Ljung-Box")$p.value
    
  }
  
  for(m in 1:12){
    # output: 12 month quantile regression for min, max, mean
    index <- Ind(m)
    
    temp1.m <- temp1[index]
    temp2.m <- temp2[index]
    temp3.m <- temp3[index]
    n <- length(index)
    co2.log_ave.m <- rep(NA, n)
    
    if(m %in% c(1,2,3)){
      for (k in 2:n){
        co2.log_ave.m[k] <- 1/4*(log(co2.month[(index[k]-3)]) + log(co2.month[(index[k]-2)]) + log(co2.month[(index[k]-1)]) + log(co2.month[index[k]]))
      }
    }else{
      for (k in 1:n){
        co2.log_ave.m[k] <- 1/4*(log(co2.month[(index[k]-3)]) + log(co2.month[(index[k]-2)]) + log(co2.month[(index[k]-1)]) + log(co2.month[index[k]]))
      }
    }
    
    
    # change tau
    ## min
    if (!anyNA(temp1.m)){
      # quantile by each month
      if(m %in%  c(1,2,3)){
        fit1 <- rq(temp1.m[2:n] ~ co2.log_ave.m[2:n], tau = 0.05) 
        QR.b <- boot.rq(cbind(1,co2.log_ave.m[2:n]),temp1.m[2:n],tau=0.05, R=10000)
      }else{
        fit1 <- rq(temp1.m[1:n] ~ co2.log_ave.m[1:n], tau = 0.05) 
        QR.b <- boot.rq(cbind(1,co2.log_ave.m[1:n]),temp1.m[1:n],tau=0.05, R=10000)
      }
      
      coef0 <- as.numeric(fit1$coefficients[1])
      coef1 <- as.numeric(fit1$coefficients[2])
      # bootstrap se
      res1 <- summary.rq(fit1, se="boot")
      se1 <- as.numeric(res1$coefficients[2,2])
      pvalue1 <- as.numeric(res1$coefficients[2,4])
      # confidence interval 95%
      conf1 <- t(apply(QR.b$B, 2, quantile, c(0.025,0.975)))
      lb1 <- conf1[2,1]
      ub1 <- conf1[2,2]
      
      int_min[m] <- coef0
      coef_min[m] <- coef1
      se_min[m] <- se1
      lb_min[m] <- lb1
      ub_min[m] <- ub1
    }
  }
  
  return(list(int_min=int_min, coef_min=coef_min, se_min=se_min, lb_min=lb_min, ub_min=ub_min))
}

## 3.3 Generate the time series for north american grids
x0 = 15:36
y0 = 45:57
nx <- length(x0)
ny <- length(y0)

locations = expand.grid(x0,y0)
nloc = nrow(locations)

ts_min_results <- list()

for (k in 1:nloc){
  results <- ts.min(locations[k,])
  ts_min_results[[k]] <- results
}

## find NA locations
ind <- c()

for (i in 1:nloc){
  if (!is.null(ts_min_results[[i]]$coef_min) & length(ts_min_results[[i]]$coef_min)==12){
    ind <- c(ind, i)
  }
}

lon_n <- lon_id[locations[ind,1]]  # the original index
lat_n <- lat_id[locations[ind,2]]

ind1 <- ind[-5] ## 121

##############################################################
## 4.  Spatial clustering
library(fda.usc)

mlearn <- X  #quantile regression coefficient matrix generated above

# k = 3
out.fd1 = kmeans.fd(mlearn,ncl=3,draw=TRUE)
cluster1 <- as.numeric(out.fd1$cluster)

## plot the cluster results

## cluster group = 3 of functional data: 121 grids
cluster_m1 <- matrix(NA, ncol=73, nrow=96)

for(i in 1:length(ind1)){
  xx <- locations[ind1[i], 1]  #the adjusted index
  yy <- locations[ind1[i], 2]
  
  cluster_m1[xx,yy] <- cluster1[i]
}

#####################################################################
## Question: How to create the categorical color scale instead of continuous scalings?
my.cols_cate3 <- function(dat, nlevels, col = c("blue","white","red"), na.color="#D3D3D3"){
  n = nlevels
  z = c(dat)
  z1 = z[!is.na(z)]
  zlim = range(z,na.rm=TRUE)
  #n1 = floor(n * (length(z1[z1<0.05])) / length(z1))
  n1 = 1
  n2 = 1
  n3 = 1
  col1 = colorRampPalette(col[1])(n1)
  col2 = colorRampPalette(col[2])(n2)
  col3 = colorRampPalette(col[3])(n3)
  cols1 = c(col1, col2,col3)
  
  zstep <- (zlim[2] - zlim[1]) / length(cols1); # step in the color palette
  newz.na <- zlim[2] + 2 * zstep # new z for NA
  
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na
  
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na
  
  cols <- c(cols1, na.color) # we construct the new color range by including: na.color
  
  tmp1 <- ifelse(is.na(dat), newz.na, dat)
  
  return(list(cols=cols, zlim=zlim, data = tmp1))
}

res1 <- my.cols_cate3(dat=cluster_m1, nlevels=3, c("blue","white","red"), na.color="#D3D3D3")
#res1 <- my.cols(dat=cluster_m1, nlevels=400, c("blue","white","red"), na.color="#D3D3D3")
dat <- res1$data
dat1 <- dat[15:36, 45:57]
#par(mfrow = c(1, 1))
image.plot(x[15:36], y[45:57],
           dat1,las = 1, main = "Spatial Clustering, k=3",  
           col = res1$cols, zlim = res1$zlim,
           xlab = "lon", ylab = "lat")
map("world", add = T)




###################################################################
## 5. Pick one location to compare the linear quantile regression and nonparametric qr methods

i = 16
j = 53
loc1 = loc.trans(i, j)  # functions: see part 3.2
m = 6

temp.m <- ts.temp_m(loc1$i,loc1$j)
temp1 <- temp.m$ts.min

index <- Ind(m)

temp1.m <- temp1[index]

n <- length(index)
co2.log_ave.m <- rep(NA, n)

for (k in 1:n){
  co2.log_ave.m[k] <- 1/4*(log(co2.month[(index[k]-3)]) + log(co2.month[(index[k]-2)]) + log(co2.month[(index[k]-1)]) + log(co2.month[index[k]]))
}

## 5.1 linear qr method: model 1 on June
# fit at other quantile levels and plot
taus = c(0.05, 0.25, 0.5, 0.75, 0.95)
f1 <- rq(temp1.m[1:n] ~ co2.log_ave.m[1:n], tau = taus) 
f1 
for (i in 1: length(taus)){
  abline(coef(f1)[,i], col=i)
}

# plots of confidence band with different taus
taus = 1:49/50
fm = rq(temp1.m ~ co2.log_ave.m, taus)
sfm = summary(fm, se="boot", alpha=0.05)
plot(sfm, mfrow=c(1,2))


## 5.2 Nonparametric quantile regressions
par(mfrow=c(1,2))
# time
taus = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)
years = 1959:2014
plot(temp1.m ~ years,col="grey", main="Locally polynomial")
for(i in 1:length(taus)){
  fit = lprq(years,temp1.m, h=20, tau=taus[i]) # select the optimal h
  lines(fit$fv~fit$xx, col=i)
}

# compare with linear quantile regression with co2 (model 1 by month)
plot(temp1.m ~ years,col="grey", main="Linear QR")
for(i in 1:length(taus)){
  fit = rq(temp1.m ~ co2.log_ave.m, tau = taus[i])
  #abline(fit, col=i)
  lines(fit$fitted.values~years, col=i)
}

## co2
# plot(temp1.m ~ co2.log_ave.m,col="grey", main="B-spline (5 knots)")
# for(i in 1:length(taus)){
#   fit = rq(temp1.m~bs(co2.log_ave.m, knots=5), tau=taus[i])
#   lines(fit$fitted~co2.log_ave.m, col=i)
# }

## 5.3 Appendix: other comparisons---linear with nonparametric

## 0. linear quantile regression

plot(temp1.m ~ co2.log_ave.m,col="grey", main="Linear QR")
for(i in 1:length(taus)){
  fit = rq(temp1.m ~ co2.log_ave.m, tau = taus[i])
  abline(fit, col=i)
}


## 1.locally polynomial method
## co2
plot(temp1.m ~ co2.log_ave.m,col="grey", main="Locally polynomial")
for(i in 1:length(taus)){
  fit = lprq(co2.log_ave.m,temp1.m, h=0.08, tau=taus[i])
  lines(fit$fv~fit$xx, col=i)
}

## 2. B-spline method
plot(temp1.m ~ co2.log_ave.m,col="grey", main="B-spline (10 knots)")
for(i in 1:length(taus)){
  fit = rq(temp1.m~bs(co2.log_ave.m, knots=10), tau=taus[i])
  lines(fit$fitted~co2.log_ave.m, col=i)
}





#############################################################################
library(animation)
for (i in 1:12){
  QRFit_by_mon[[i]] <- apply(tmx_mon[,, i,], 1:2, function(x) rq(x ~ log_CO2_4mon[seq(0, 456, 12) + i], tau = 0.5)$coefficients)
}

rg <- range(unlist(lapply(QRFit_by_mon, range)))

saveLatex({
  for (i in 1:12){
    maps::map("world", xlim = range(lon), ylim = range(lat))
    image.plot(lon, lat, QRFit_by_mon[[i]][2,,], add = T, zlim = rg)
    maps::map("world", add = T)
    maps::map("state", add = T)
    title(month.abb[i])
  }
}, img.name = "QR_max_model1", ani.opts = "controls,width=0.95\\textwidth",
latex.filename = ifelse(interactive(), "QR_max_model1_us.tex", ""), interval = 0.5,
nmax = 12, ani.dev = "pdf", ani.type = "pdf", ani.width = 10,
ani.height = 8,documentclass = paste("\\documentclass{article}",
                                     "\\usepackage[papersize={10in,8in},margin=0.05in]{geometry}",
                                     sep = "\n"))
