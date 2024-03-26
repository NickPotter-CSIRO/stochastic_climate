# NB basepath variable defined outside this script. 

# Helper functions for the ARMA modelling:
xtract_model <- function(arma.model){
  nar = arma.model$arma[1]
  nma = arma.model$arma[2]
  if ((nar) == 0){
    return(list(ar=c(),
                ma=arma.model$coef[1:nma]))
    
  } else{
    if ((nma) == 0){
      return(list(ar=arma.model$coef[1:nar],
                  ma = c()))
    } else{
      return(list(ar=arma.model$coef[1:nar],
                  ma=arma.model$coef[(nar+1):(nar+nma)]))
    }
    
  }
}

which.cluster_singleyear <- function(x_, centres){
  # compute squared euclidean distance from a sample to each cluster center
  x_ = as.numeric(x_)
  #tmp = apply(centers, 1,
  #            function(v) sum((x_-v)^2))
  
  tmp = apply(sweep(centres,2,x_,FUN="-")**2,1,sum)
  #hist(tmp)
  return(which.min(tmp))
}
which.cluster <- function(x, centers) {
  # compute squared euclidean distance from each sample to each cluster center
  tmp <- sapply(seq_len(nrow(x)),
                function(i) apply(centers, 1,
                                  function(v) sum((x[i, ]-v)^2)))
  max.col(-t(tmp))  # find index of min distance
}  

get_MoF_fragyrs <- function(asim, compare, SIM_START_YEAR){
  #C = NULL
  ##browser()
  nC = dim(asim)[2]
  frags = data.frame(Year = numeric(), fragyr = numeric())
  
  # scale simulation
  if (nC == 1){
    annual_sd = sd(compare[,2])
    asim.scaled = asim / annual_sd
  } else{
    annual_sds = apply(compare[,2:(nC+1)], 2, sd)
    asim.scaled = sweep(asim, 2, annual_sds, FUN="/")
  }
  # scale comparison values
  compare.scaled = cbind(data.frame(Year=compare[,1]),
                         sweep(compare[,-1], 2, 
                               apply(compare[,-1],2,sd),
                               FUN="/"))
  for (colname in names(compare)[2:(1+2*nC)]){
    frags[colname] = numeric()
  }
  nsim = 0
  sim_year = SIM_START_YEAR + nsim
  if (nC==1){
    fragix0 = which.min((asim.scaled[1,] - compare.scaled[,2])**2)
  } else{
    fragix0 = which.cluster_singleyear(asim.scaled[1,],
                                       compare.scaled[,2:(nC+1)])
  }
  #  stop('here')
  fragyr = compare$Year[fragix0]
  cc = compare[fragix0,-1]
  csc=compare.scaled[fragix0,-1]
  
  #print(data.frame(Year=sim_year,
  #                 fragyr=fragyr))
  newrow = cbind(data.frame(Year=sim_year,
                            fragyr=fragyr),
                 csc)
  frags = rbind(frags, newrow)
  for (nsim in 1:(dim(asim)[1]-1)){
    sim_year = sim_year + 1
    fragix = which.cluster_singleyear(cbind(data.frame(t(asim.scaled[nsim+1,])),
                                            frags[nsim,(nC+3):(2*nC+2)]),
                                      compare.scaled[,-1])
    #matplot(t(compare.scaled[,2:(2*nC+1)]))
    #lines(as.numeric(cbind(data.frame(t(asim.scaled[nsim+1,])), frags[nsim, (nC+3):(2*nC+2)])))
    #lines(as.numeric(compare.scaled[fragix,2:(2*nC+1)]), col='red')
    fragyr = compare$Year[fragix]
    cc = compare[fragix,-1]
    csc = compare.scaled[fragix,-1]
    newrow = cbind(data.frame(Year=sim_year,
                              fragyr=fragyr),
                   csc)
    frags = rbind(frags, newrow)  
    #print(frags)
  }
  #frags
  # Multiply the standard deviation back in
  frags_unscaled = cbind(frags[,1:2],
                         sweep(frags[,-c(1:2)], 2, 
                               apply(compare[,-1],2,sd),
                               FUN="*"))
  return(frags_unscaled)  
}


MoF_annual <- function(asim, compare, SIM_START_YEAR, nClusters=NA){
  nC = dim(asim)[2]
  #cat('Number of columns: ',nC,'\n')
  if (is.na(nClusters)){
    # no clusters - regular method of fragments
    unclustered_res = get_MoF_fragyrs(asim, compare, SIM_START_YEAR)
    return(list(results=unclustered_res,
                clusters=NULL))
  }
  else{
    nCluster = NA
    comp = NA
    comp.scaled = NA
    ##browser()
    compare.scaled = cbind(data.frame(Year=compare[,1]),
                           sweep(compare[,-1], 2, 
                                 apply(compare[,-1],2,sd),
                                 FUN="/"))
    clusters = kmeans(compare.scaled[,2:(2*nC+1)],
                      nClusters)
    cluster_ixes = lapply(1:nClusters, function(x) which(clusters$cluster==x))
    cluster_comp.scaled = cbind(data.frame(Year=1:nClusters),
                                clusters$centers)
    cluster_comp = cbind(data.frame(Year=cluster_comp.scaled$Year), 
                         sweep(cluster_comp.scaled[,-1], 2, 
                               apply(compare[,2:(2*nC+1)], 2, sd), 
                               FUN="*"))
    clustered_res = get_MoF_fragyrs(asim, cluster_comp, SIM_START_YEAR)
    
    # Take random samples within clusters:
    random_fragixes = rep.int(NA, dim(clustered_res)[1])
    for (i in 1:nClusters){
      ixes = clustered_res$fragyr==i
      random_fragixes[ixes] = sample(cluster_ixes[[i]], 
                                     size = sum(ixes),
                                     replace = T)
      
    }
    random_fragyrs = compare$Year[random_fragixes]
    clustered_res$fragyr = random_fragyrs
    return(list(results=clustered_res,
                clusters=clusters))
    
  }
}

### Date functions

ndays.month <- function(month, year){
  if (month %in% c(1,3,5,7,8,10,12)){
    return(31)
  } else{
    if (month == 2){
      if (is.leapyear(year)){
        return (29)
      } else{
        return (28)
      }
    } else{
      return (30)
    }
  }
  
  stop('Exception in ndays.month function')
}

is.leapyear <- function(year){
  # is year leap-year or not?
  if ((year %% 4) == 0) {
    if ((year %% 100) == 0) {
      if ((year %% 400) == 0) {
        return (TRUE)
      } else {
        return (FALSE)
      }
    } else {
      return (TRUE)
    }
  }
  return (FALSE)
}


# Date functions:
calendar2water <- function(cyear, cmonth, start_month = 4){
  temp = data.frame(year = cyear,
                    month = cmonth)
  if (start_month == 1){
    temp$wyear = temp$year
    temp$wmonth = temp$month - 1
  } else{
    if (start_month %in% seq(2,12)){
      temp$wyear = temp$year
      temp$wyear[temp$month %in% seq(1, start_month-1)] = 
        temp$wyear[temp$month %in% seq(1, start_month-1)] - 1 
      temp$wmonth = (temp$month - start_month) %% 12
      
    } else{
      stop('calendar2water(): start_month must be in [1, 12]')
    }
  }
  return(temp[,c('wyear', 'wmonth')])
}
water2calendar <- function(wyr, wmth, start_month = 4){
  temp = data.frame(wyear=wyr, wmonth=wmth)
  if (start_month == 1){
    temp$cmonth = temp$wmonth + 1
    temp$cyear = temp$wyear
  } else{
    if (start_month %in% seq(2,12)){
      temp$cmonth = (temp$wmonth + (start_month - 1)) %% 12 + 1
      temp$cyear = temp$wyear
      temp$cyear[temp$wmonth >= (12 - start_month + 1)] = 
        temp$cyear[temp$wmonth >= (12 - start_month + 1)] + 1
    } else{
      stop('water2calendar(): start_month must be in [1, 12]')
    }
  }
  return (temp[,c('cyear', 'cmonth')])
}

# Main function to disaggregate from annual to daily
# NB this is wrapped up in Simulation() methods now


disaggregate <- function(annual_frags, asim, 
                         aacd=all_annual_climate_data, 
                         adcd=all_daily_climate_data){
  # annual_frags contains: sim year in 1st column (typically a seq-type variable),
  #                      , frag year in 2nd column
  # asim is a (NSIMYEARS x nC) matrix of annual results
  # aacd is a copy of all_annual_climate_data
  # adcd is a copy of all_daily_climate_data
  
  nC = dim(asim)[2] # number of catchments
  
  # merge simulated annual data and observed annual rainfall from frag years
  df = data.frame(annual_frags[,1:2], 
                  asim) %>%  # Just in case column names have "X" preprended
                    dplyr::rename_with(function(xx) gsub(x=xx,pattern='X',replacement=''))
  if (all(grepl(pattern='X', colnames(asim)))){
    if (any(!grepl(pattern='X', colnames(asim)))){
      stop("Some columns contain X, others don't")
    }
    X_present = T
    
  } else{
    X_present = F
  }
  #names(df)[-c(1,2)] = gsub(names(df)[-c(1,2)], pattern="X",replacement="")
  df2=merge(df, aacd %>% dplyr::select(matches('Year') |
                                         starts_with('Rain')),
            by.x='fragyr', by.y='Year')
  # Calculate scaling factors (simulated annual divided by observed frag year rainfall)
  df2.m = melt(df2, id.vars=c('fragyr','Year'))
  df2.m$Obs = paste0("Obs_",grepl(df2.m$variable, pattern='Rain.'))
  df2.m$variable = gsub(df2.m$variable, pattern='Rain.X',replacement='')
  df2.m$variable = gsub(df2.m$variable, pattern='Rain.',replacement='')
  
  df2.c = dcast(df2.m, Year + fragyr + variable ~ Obs) %>% 
    dplyr::rename(c('Sim.Year'='Year'))
  df2.c$scaling_factor = df2.c$Obs_FALSE / df2.c$Obs_TRUE
  
  # Take a copy of daily data from frag years (for analogue simulation)
  daily.rain = adcd %>% dplyr::select(matches('date') |
                                        starts_with('Rain')) %>% 
    mutate(Year=as.numeric(substr(date,1,4)), 
           Month=as.numeric(substr(date,6,7)), 
           Day=as.numeric(substr(date, 9,10)))
  daily.rain.sim = merge(annual_frags[,1:2] %>% dplyr::rename(c('Sim.Year'='Year')), 
                         daily.rain, by.x='fragyr',by.y='Year')
  daily.rain.sim.m = melt(daily.rain.sim, id.vars=c('fragyr','Sim.Year','date','Month','Day'))
  daily.rain.sim.m$variable = gsub(daily.rain.sim.m$variable, 
                                   pattern = 'Rain.X',
                                   replacement='')
  daily.rain.sim.m$variable = gsub(daily.rain.sim.m$variable, 
                                   pattern = 'Rain.',
                                   replacement='')
  #print(daily.rain.sim.m %>% head)
  #browser() 
  # merge with scaling factor and scale:
  mg2 = merge(daily.rain.sim.m, df2.c[,c('Sim.Year','variable','scaling_factor')],
              by=c('Sim.Year', 'variable'))
  mg2$scaled = mg2$value * mg2$scaling_factor

  
  daily.sim = dcast(mg2, Sim.Year + Month + Day ~ variable, value.var='scaled')
  # Note that these columns are not in the same order as those from asim
  # reorder them...
  if (X_present){ # put 'X' back in front of column names
    names(daily.sim)[4:dim(daily.sim)[2]] = paste0('X',names(daily.sim)[4:dim(daily.sim)[2]])
  }
  daily.sim = daily.sim %>% dplyr::relocate(names(daily.sim)[1:3], 
                                            colnames(asim))

  
  # Fix up Feb 28s and 29s
  nd = table(subset(daily.sim, daily.sim$Month == 2)$Sim.Year)
  nd2 = data.frame(nd)
  names(nd2) = c('Sim.Year','ndays.sim')
  nd2$Sim.Year = as.numeric(as.character(nd2$Sim.Year))
  
  nd2$ndays.actual = sapply(nd2[,1],FUN=is.leapyear)*1+28
  #browser()
  # Just add a zero rainfall day on 29th Feb that shouldn't be there
  toadd = nd2[which(nd2$ndays.sim<nd2$ndays.actual),'Sim.Year']
  feb29s = data.frame(matrix(rep.int(0,length(toadd)*nC), nrow=length(toadd)))
  names(feb29s) = names(daily.sim)[-c(1:3)]
  feb29s.df = cbind(data.frame(Sim.Year=toadd, Month=2, Day=29), feb29s)
  daily.sim = rbind(daily.sim, feb29s.df)
  
  # Delete extra Feb 29ths 
  toremove = nd2[which(nd2$ndays.sim>nd2$ndays.actual),'Sim.Year']
  dropix = which((daily.sim$Month == 2) & (daily.sim$Day==29) & (daily.sim$Sim.Year %in% toremove))
  # Copy rainfall from extra 29ths onto the day before to preserve monthly/annual totals
  replaceix = dropix - 1
  daily.sim[replaceix,-c(1:3)] = daily.sim[replaceix,-c(1:3)] + daily.sim[dropix,-c(1:3)]
  # Drop 29ths
  daily.sim = daily.sim[-dropix, ]
  
  # resort by date
  daily.sim = daily.sim[order(daily.sim$Sim.Year,daily.sim$Month,daily.sim$Day),]
  
  return(daily.sim)
}



#################################
#################################
# Classes to hold rainfall and PET observations and simulations
#################################
#################################


library(R6)
Rainfall <- R6Class("Rainfall",
                    public = list(
                      list_of_data = NULL,
                      daily_data = NULL,
                      CIDs = NULL,
                      annual_values = NULL,
                      raw_vcvmat = NULL,
                      tr_vcvmat = NULL,
                      monthly_values = NULL,
                      compare = NULL,
                      ENSO_rainfall_model_parameters = NULL,
                      ENSO_model_fitted_dates = NULL,
                      initialize = function(list_of_data) {
                        self$list_of_data <- list_of_data
                        self$CIDs <- names(list_of_data)
                        self$daily_data = data.frame(YMD = character(),
                                                     Rain = numeric())
                        for (cid in self$CIDs){
                          self$daily_data = merge(x=self$daily_data,
                                                  y=self$list_of_data[[cid]],
                                                  by='YMD',all=T,
                                                  suffixes = c('',paste0('.',cid)))
                        }
                        self$daily_data$Rain = NULL
                        self$daily_data$Year = as.numeric(substr(self$daily_data$YMD, 1, 4))
                        self$daily_data$Month = as.numeric(substr(self$daily_data$YMD, 6, 7))
                        self$daily_data$Day = as.numeric(substr(self$daily_data$YMD, 9, 10))
                      },
                      
                      monthly_aggregation = function(){#whole.months.only = TRUE){
                        nC = length(self$CIDs)
                        self$monthly_values = aggregate(self$daily_data[,2:(1+nC)],
                                                        by = list(Year=self$daily_data$Year,
                                                                  Month=self$daily_data$Month),
                                                        FUN=sum)
                        #if (whole.months.only){
                        #  month_counts = aggregate(self$daily_data[,2:(1+nC)],
                        #                           by = list(Year=self$daily_data$Year,
                        #                                     Month=self$daily_data$Month),
                        #                           FUN=count)
                        #  month_days = mapply(ndays.month,
                        #                      self$monthly_values$Month,
                        #                      self$monthly_values$Year)
                        #  self$monthly_values = self$monthly_values[which(month_counts == month_days)]
                        #}
                      },
                      annual_aggregation = function(){
                        nC = length(self$CIDs)
                        self$annual_values = aggregate(self$daily_data[,2:(1+nC)],
                                                       by = list(Year=self$daily_data$Year),
                                                       FUN=sum)
                      },
                      make_compare = function(){
                        if (is.null(self$monthly_values)){
                          self$monthly_aggregation()
                        }
                        if (is.null(self$annual_values)){
                          self$annual_aggregation()
                        }
                        compare = merge(self$annual_values,
                                        self$monthly_values %>% 
                                          subset(Month==12) %>% # filter to December
                                          dplyr::select(-Month),       # drop Month column
                                        by='Year',
                                        suffixes = c('.Annual','.Dec'))
                        self$compare = compare[complete.cases(compare),]
                        ### Include complete years only...
                        ndays_year = tapply(rainfall_obs$daily_data$Year, rainfall_obs$daily_data$Year, FUN=length)
                        ndays_year = data.frame(yr=as.numeric(names(ndays_year)), 
                                                ndays_data = ndays_year) 
                        ndays_year$ndays_actual = 365 + 1 * sapply(ndays_year$yr, FUN=is.leapyear)
                        ndays_year$include = ndays_year$ndays_data == ndays_year$ndays_actual
                        incl_years = ndays_year[which(ndays_year$include),'yr']
                        self$compare = self$compare[which(self$compare$Year %in% incl_years),]
                        ###
                      },
####################### Fit ENSO-rainfall model #######################
                      fit_ENSO_model = function(ENSO_Obj_Instance,
                                                start_year=NULL,end_year=NULL){
                        library(MASS) # boxcox
                        library(forecast) # BoxCox
                        r = 1
                        ENSO_rainfall_model_parameters_ = data.frame(CID = numeric(),
                                                                     raw_intercept = numeric(),
                                                                     raw_slope = numeric(),
                                                                     raw_s2 = numeric(),
                                                                     Lambda = numeric(),
                                                                     transformed_intercept = numeric(),
                                                                     transformed_slope = numeric(),
                                                                     transformed_s2 = numeric())
                        if (is.null(self$annual_values) | is.null(self$annual_values)){
                          stop('Error in fit_ENSO_model(): annual_values missing, need to call annual_aggregation()?')
                        }
                        mgd = merge(self$annual_values, 
                                    ENSO_Obj_Instance$annual_data,
                                    by = 'Year')
                        if (!is.null(start_year) & !is.null(end_year)){
                          mgd = mgd %>% dplyr::filter(Year >= start_year &
                                                        Year <= end_year)
                        }
                        self$ENSO_model_fitted_dates = c(min(mgd$Year),
                                                         max(mgd$Year))
                        for (cid in self$CIDs){
                          col = paste0('Rain.', cid)
                          # First fit model to untransformed rainfall:
                          Fo = formula(paste(paste0(col, ' ~ ENSO_annual')))
                          
                          linear_fit = lm(Fo, data=mgd)
                          ENSO_rainfall_model_parameters_[r,1] = cid
                          ENSO_rainfall_model_parameters_[r,2:3] = linear_fit$coefficients
                          ENSO_rainfall_model_parameters_[r,4] = anova(linear_fit)$`Mean Sq`[2]
                          
                          # ML estimate of box-cox parameter:
                          tmp = boxcox(Fo, data=mgd,
                                       lambda=seq(-2,2,1/100),plotit=F)
                          Lam = tmp$x[which.max(tmp$y)]
                          ENSO_rainfall_model_parameters_[r,5] = Lam
                          mgd[,paste0('bc',col)] = BoxCox(mgd[,col],
                                                          Lam)
                          
                          # fit model to BC-transformed rainfall:
                          bcFo = formula(paste0('bc',col,'~ENSO_annual'))
                          transformed_linear_fit = lm(bcFo, data=mgd)

                          ENSO_rainfall_model_parameters_[r,6:7] = transformed_linear_fit$coefficients
                          ENSO_rainfall_model_parameters_[r,8] = anova(transformed_linear_fit)$`Mean Sq`[2]
                          
                          
                          r = r + 1
                        }
                        self$ENSO_rainfall_model_parameters = ENSO_rainfall_model_parameters_
                      
                        # Now calculate variance-covariance matrices
                        
                        self$raw_vcvmat = cov(mgd[,paste0('Rain.', self$CIDs)])
                        self$tr_vcvmat = cov(mgd[,paste0('bcRain.', self$CIDs)])
                        
                        
                      },
                      simulate = function(ENSO_ts, nR, transformed=F){# pass in a vector only, realisation 
                        # of an ENSO replicate (or observation)
##################                        # No argument? Uses observed IPO timeseries
##################                        # that the model was fitted on.
                        
                        if (is.null(self$ENSO_rainfall_model_parameters)){
                          stop('Error in simulate(): need to call fit_ENSO_model first')
                        }
                        if (transformed){
                          covmat = self$tr_vcvmat
                          
                        } else{
                          covmat = self$raw_vcvmat
                        }
                        # Previously, random errors were added in after inverse box-cox transform
                        # # Previously, random errors were added in after inverse box-cox transform
                        # gen_replicate <- function(sim_soi, bcvcvmat, lr_params = ENSO_rain_params,
                        #                           NR=1, CIDS_=CIDs){
                        #   CIDs = NA
                        nC = length(self$CIDs)
                        nY = length(ENSO_ts)
                        simresids = array(NA,
                                          c(nY, nC, nR))
                        # Generate random (correlated) errors
                        for (i in 1:nR){
                          simresids[,,i] = mvrnorm(n=nY,
                                                   mu=rep.int(0,nC),
                                                   Sigma=covmat)#self$tr_vcvmat)
                        }
                        dimnames(simresids)[[2]] = self$CIDs
                        #   
                        #bcsims = simresids # copy
                        #sims = bcsims #copy
                        #   
                        bcsims = NULL
                        ###browser()
                        for (cid in self$CIDs){
                          bcmean_rainf = self$predict_rain_from_ENSO(ENSO_ts_=ENSO_ts, 
                                                                        cid_=cid, 
                                                                        transformed=transformed)
                          bcsims = cbind(bcsims, array(bcmean_rainf, c(nY, 1)))
                        }
                        
                        # broadcast bcsims out along third dimension:
                        dim(bcsims) = c(dim(bcsims), 1)
                        if (nR > 1){
                          bcsims = bcsims[,,rep(1,nR)]
                        }
                        
                        dimnames(bcsims)[[2]] = self$CIDs
                        
                        sims = bcsims + simresids
                        for (xcid in dimnames(bcsims)[[2]]){
                          ix = which(as.character(self$ENSO_rainfall_model_parameters$CID) == as.character(xcid))
                          if (transformed){
                            sims[,xcid,] = InvBoxCox(sims[,xcid,],
                                                     #biasadj=TRUE,
                                                     #fvar=self$ENSO_rainfall_model_parameters[ix,'raw_s2'],
                                                     lambda = self$ENSO_rainfall_model_parameters[ix,'Lambda'])
                          }
                        }
                        return(sims)
                      },
                      predict_rain_from_ENSO = function(ENSO_ts_, cid_, transformed=F){
                        ENSO_ts = NULL
                        cid = NULL
                        if (transformed){
                          col_prefix = 'transformed_'
                        } else{
                          col_prefix = 'raw_'
                        }
                        #cid_n = as.numeric(cid_)
                        
                        ix = which(as.character(self$ENSO_rainfall_model_parameters[,'CID']) == as.character(cid_))
                        # print(as.character(self$ENSO_rainfall_model_parameters[,'CID']))
                        # print(as.character(cid_))
                        # print(ix)
                        # print(self$ENSO_rainfall_model_parameters[ix,])
                        return(ENSO_ts_ * self$ENSO_rainfall_model_parameters[ix, paste0(col_prefix, 'slope')] +
                                 self$ENSO_rainfall_model_parameters[ix, paste0(col_prefix, 'intercept')])
                      }
                    ), # End public methods
                    private = list(
                      
                    ) # End private methods
)

PETData <- R6Class("PETData",
                   public = list(
                     list_of_data = NULL,
                     daily_data = NULL,
                     CIDs = NULL,
                     cosine_parameters = NULL,
                     residual_model = NULL,
                     #annual_values = NULL,
                     #monthly_values = NULL,
                     #compare = NULL,
                     initialize = function(list_of_data) {
                       self$list_of_data <- list_of_data
                       self$CIDs <- names(list_of_data)
                       self$daily_data = data.frame(YMD = character(),
                                                    PET = numeric())
                       for (cid in self$CIDs){
                         self$daily_data = merge(x=self$daily_data,
                                                 y=self$list_of_data[[cid]],
                                                 by='YMD',all=T,
                                                 suffixes = c('',paste0('.',cid)))
                       }
                       self$daily_data$PET = NULL
                       self$daily_data$Year = as.numeric(substr(self$daily_data$YMD, 1, 4))
                       self$daily_data$Month = as.numeric(substr(self$daily_data$YMD, 6, 7))
                       self$daily_data$Day = as.numeric(substr(self$daily_data$YMD, 9, 10))
                       self$daily_data$jday_offset0 = as.numeric(as.Date(self$daily_data$YMD) - 
                                                                   as.Date(paste0(self$daily_data$Year, '-01-01')))
                     },
                     fit_cosine = function(start_year, end_year, envelope=T){
                       ix = which(self$daily_data$Year >= start_year & 
                                    self$daily_data$Year <= end_year)
                       self$cosine_parameters = list(asdf=1)
                       
                       
                       
                       PET_model_error = function(pars, CID, envelope){
                         if (envelope){
                            return (sum((fy_envelope - self$PETf(fx_envelope, pars))**2))
                         
                         } else{
                           predictions = self$PETf(j, pars)
                           return(sum((predictions - d2)**2,
                                      na.rm=T))
                         }
                       }


                       for (cid_ in self$CIDs){    
                         j = self$daily_data[ix, 'jday_offset0']  
                         d2 = self$daily_data[ix,
                                              paste0('PET.',
                                                     cid_)]
                         #browser()
                         # First calculate 95% envelope of data:
                         if (envelope){
                           f_envelope = tapply(X=d2,
                                               INDEX=self$daily_data$jday_offset0[ix],
                                               FUN=quantile, probs=c(0.95))
                           fx_envelope = as.numeric(names(f_envelope))
                           fy_envelope = as.numeric(f_envelope)
                         }
                         init_params = c(mean(d2, na.rm=T), # param 1 is mean of cosine curve
                                       mean(d2, na.rm=T), # param2 is amplitude
                                       # param3 is occurrence of maximum value (jday)
                                       self$daily_data$jday_offset0[ix][which(d2 == max(d2,na.rm=T))[1]])
                       
                         # 
                         res = optim(par=init_params,
                                     fn = PET_model_error,
                                     CID=cid_,
                                     envelope=envelope)
                         self$cosine_parameters[[as.character(cid_)]] = res$par
                         
                         self$cosine_parameters[['asdf']] = NULL
                       }
                       
                     },
                     fit_model = function(Rainfall_Observations_Instance, start_year, end_year,
                                          span=0.82){
                       if (is.null(self$cosine_parameters)){
                         stop('Error in fit_model(): need to call fit_cosine() first')
                       }
                       ix = which(self$daily_data$Year >= start_year & 
                                    self$daily_data$Year <= end_year)
                       pix = which(Rainfall_Observations_Instance$daily_data$Year >= start_year & 
                                     Rainfall_Observations_Instance$daily_data$Year <= end_year)
                       
                       resids = self$daily_data[ix, c('YMD','jday_offset0')]
                       for (cid in self$CIDs){
                         resids[,paste0('PET_resid.', cid)] = self$daily_data[ix, paste0('PET.', cid)]/ # - 
                           self$PETf(self$daily_data[ix, 'jday_offset0'], self$cosine_parameters[[as.character(cid)]])
                       }
                       
                       self$residual_model = list(asfd=1)
                       for (cid in self$CIDs){
                         cat('---',cid,'\n')
                         
                         cat('---',Rainfall_Observations_Instance$daily_data %>% names,'\n')
                         rain_obs_cid = Rainfall_Observations_Instance$daily_data[pix, paste0('Rain.', cid)]
                         #plot(resids[,paste0('PET_resid.',cid)] ~ rain_obs_cid)
                         m1 = loess(resids[,paste0('PET_resid.',cid)] ~ rain_obs_cid,
                                    span=span)
                         pr = list()
                         # Take 99 quantiles of non-zero rainfall and guards for zero and very large rainfall.
                         pr$x = c(0,
                                  as.numeric(quantile(Filter(function(x) x > 0, rain_obs_cid), 
                                                      probs = seq(from=0,to=0.99,length.out=99))))
                         #10**24)
                         pr$x[1] = -10**-10 # This allows zero rainfall days in the call to "approx" when simulating
                         pr$y = predict(m1, newdata=pr$x)
                         pr$y[1] = 1
                         pr$x[length(pr$x) + 1] = 10**24 # Extrapolation: anything very large just uses the last value
                         pr$y[length(pr$y) + 1] = pr$y[length(pr$y)]
                         lines(pr$x, pr$y, col='red')
                         self$residual_model[[as.character(cid)]] = pr
                         #print(m1)
                       }
                       self$residual_model[['asfd']] = NULL
                       #return(m1)
                     },
                     PETf = function(x, pars) { 
                       
                       pars[1] + pars[2] * cos((x - pars[3])/366*2*pi)
                     },
                     
                     
                     model_PET = function(Rainfall_Observations_Instance){
                       if (is.null(self$residual_model)){
                         print('Error in model_PET(): need to call fit_model() first')
                       }
                       rain_cols = Filter(f = function(zz) grepl(pattern='Rain.', fixed=T,
                                                                 x = zz),
                                          x=names(Rainfall_Observations_Instance$daily_data))
                       results = Rainfall_Observations_Instance$daily_data[,c('YMD', 'Year', 'Month', 'Day')]
                       results$jday_offset0 = as.numeric(as.Date(Rainfall_Observations_Instance$daily_data$YMD) - 
                                                           as.Date(paste0(Rainfall_Observations_Instance$daily_data$Year, '-01-01')))
                       for (col in rain_cols){
                         print('----------')
                         print(col)
                         cid_str = gsub(col, pattern='Rain.', replacement='', fixed=T)
                         print(cid_str)
                         print('----------')
                         fitted_cos = self$PETf(results$jday_offset0, self$cosine_parameters[[cid_str]]) 
                         print(self$residual_model)
                         print(cid_str)
                         print(self$residual_model[[cid_str]])
                         print('----------')
                         rainfall_residual = approx(x=self$residual_model[[cid_str]]$x,
                                                    y=self$residual_model[[cid_str]]$y,
                                                    xout = Rainfall_Observations_Instance$daily_data[,col])$y
                         # Extrapolation: just use the largest value in the table:
                         rainfall_residual[which(Rainfall_Observations_Instance$daily_data[,col] > 
                                                   max(self$residual_model[[cid_str]]$x))] = 
                           max(self$residual_model[[cid_str]]$y)
                         results[,gsub(col, pattern='Rain.', replacement='modelledPET.', fixed=T)] = 
                           fitted_cos * rainfall_residual
                         
                         
                       }
                       return(results)
                     },
                     generate_PET = function(simulated_rainfall){
                       # Generates PET for rainfall timeseries 
                       # pass through simulated_rainfall, an output from 
                       # Simulation$disagg.
                       # This will be a data.frame with 
                       # Sim.Year, Month, Day and CIDs as columns
                       datestr = sprintf('%d-%02d-%02d',
                                         simulated_rainfall$Sim.Year,
                                         simulated_rainfall$Month,
                                         simulated_rainfall$Day)
                       jday = private$datestr2jday(datestr)
                       pet_sim = simulated_rainfall
                       pet_sim[,self$CIDs] = NA
                       
                       for (cid in self$CIDs){
                       #for (ic in 4:dim(pet_sim)[2]){
                         ccid = as.character(cid)
                         ic = which(names(pet_sim)==ccid)
                         
                         #cid = names(pet_sim)[ic]
                         #print(ccid)
                         #print(self$residual_mode[[ccid]])
                         
                         pet_sim[,ic] = self$PETf(jday, self$cosine_parameters[[cid]])
                         PET_rainfall_factors = approx(x=self$residual_model[[ccid]]$x, 
                                                       y=self$residual_model[[ccid]]$y, 
                                                       xout = simulated_rainfall[,ic])$y
                         pet_sim[,ic] = pet_sim[,ic] * PET_rainfall_factors
                         
                       }
                       return(pet_sim)
                     }
                     
                   ), # end public functions
                   private = list(
                     datestr2jday = function(datestr){
                       #return(as.Date(datestr))
                       #return()
                       return(as.numeric(as.Date(datestr) - 
                                           as.Date(paste0(substr(datestr,1,4), '-01-01'))))
                       
                     }
                     
                   ) # end private functions
)


Simulation <- R6Class("Simulation",
                      public = list(
                        CIDs = NULL,
                        annual_values = NULL,
                        fragyrs = NULL,
                        fragvals = NULL,
                        clusters = NULL,
                        SIM_START_YEAR = NULL,
                        daily_sim = NULL,
                        initialize = function(CIDs = NA, annual_values = NA, SIM_START_YEAR = 2030) {
                          self$CIDs <- CIDs
                          self$annual_values <- annual_values
                          self$SIM_START_YEAR <- SIM_START_YEAR
                        },
                        get_frags = function(Observations_Instance, SIM_START_YEAR=self$SIM_START_YEAR, 
                                             n_clusters = NA, restricted_CIDs = NULL){
                          if (is.null(Observations_Instance$compare)){
                            stop('Error in Simulation: need to call make_compare() on Observations_Instance')
                            
                          }
                          if (is.null(restricted_CIDs)){ # calculate on all CIDs
                            tmp = MoF_annual(asim=self$annual_values, 
                                             compare=Observations_Instance$compare, 
                                             SIM_START_YEAR = SIM_START_YEAR,
                                             nClusters=n_clusters)
                            return(tmp)
                          } else{ # calculate only using CIDs in restricted_CIDs
                            restricted_sim = self$annual_values[,colnames(self$annual_values) %in% 
                                                                  as.character(restricted_CIDs)]
                            if (length(restricted_CIDs) == 1){
                              restricted_sim = matrix(restricted_sim, ncol=1)
                            }
                            tmp = MoF_annual(asim=restricted_sim,
                                             compare=Observations_Instance$compare %>% 
                                                       dplyr::select(matches(paste0('Year|',
                                                                                    paste(restricted_CIDs,
                                                                                          collapse="|")))),
                                             SIM_START_YEAR=SIM_START_YEAR,
                                             nClusters=n_clusters)
                            return(tmp)
                          }
                        },
                        disagg_annual = function(Observations_Instance, SIM_START_YEAR=self$SIM_START_YEAR, 
                                                 n_clusters = NA, restricted_CIDs = NULL){
                          frags=self$get_frags(Observations_Instance, 
                                               SIM_START_YEAR,
                                               n_clusters, 
                                               restricted_CIDs)
                          self$fragyrs = frags$results$fragyr
                          # Save a copy of the frag results merged with the observational analogue years for
                          # future reference
                          self$fragvals = merge(frags$results %>% 
                                                  mutate(Sim.Year=Year) %>% 
                                                  dplyr::select(-Year) %>% 
                                                  mutate(fragyr=fragyr), 
                                                Observations_Instance$annual_values,by.x='fragyr',by.y='Year') %>% 
                            dplyr::arrange(Sim.Year)
                          self$clusters = frags$clusters
                          
                        },
                        disagg = function(Observations_Instance){
                          
                          if (is.null(self$fragyrs)){
                            stop("Error in Simulation: disagg_annual() needs to be called before disagg().\n  Did you call get_frags() (doesn't set attributes)")
                          }
                          nY = length(self$fragyrs) # number of years of simulation
                          sim_years = seq(from = self$SIM_START_YEAR,
                                          by = 1,
                                          length.out = nY)
                          
                          return(disaggregate(data.frame(Year=sim_years,
                                                         fragyr=self$fragyrs),
                                              self$annual_values,
                                              aacd=Observations_Instance$annual_values,
                                              adcd=Observations_Instance$daily_data %>% 
                                                mutate(date=YMD)))
                        }
                      )
)

###### 
# Climate objects

IPO_Obs <- R6Class("IPO_Obs",
                   public = list(
                     annual_data = NULL,
                     model_nointercept = NULL,
                     model_intercept = NULL,
                     year_filter = NULL,
                     initialize = function(fpath=paste0(basepath, #="//fs1-cbr.nexus.csiro.au/{lw-hydroclimate}/work/
                                                        "pot063/SEA_stoch/new/tpi.timeseries.ersstv5.data"),
                                           start_year = 1900,
                                           end_year = 2030,
                                           start_month = 1){
                       IPOwide = read.table(fpath,
                                        skip=1)
                       names(IPOwide) = c('YEAR', 1:12)
                       IPOlong = IPOwide %>% melt(id.vars=1)
                       cd=calendar2water(IPOlong$YEAR, 
                                         as.numeric(as.character(IPOlong$variable)),
                                         start_month=start_month)
                       IPOlong$wyear = cd$wyear
                       IPOlong$wmonth = cd$wmonth

                       summed = tapply(X=IPOlong$value, INDEX=IPOlong$wyear,FUN=sum)
                       annual = data.frame(Year=as.numeric(row.names(summed)),
                                           IPO_annual=as.vector(summed))
                       
                       self$year_filter = c(start_year, end_year)
                       self$annual_data = annual %>% 
                         dplyr::filter(Year >= start_year & Year <= end_year)
                     },
                     fit = function(){
                       ipomod = auto.arima(as.vector(self$annual_data$IPO_annual - mean(self$annual_data$IPO_annual)),
                                           max.order=20)
                       self$model_nointercept = ipomod
                       self$model_intercept = mean(self$annual_data$IPO_annual)
                     },
                     fsim = function(n){
                       if (is.null(self$model_nointercept)){
                         stop('Error in fsim(): need to fit model first')
                       }
                       zeromean_simulated = arima.sim(model=xtract_model(self$model_nointercept), 
                                               n=n, 
                                               sd=sqrt(self$model_nointercept$sigma2))

                       simulated = zeromean_simulated + self$model_intercept
                         
                       return(simulated)  
                     },
                     plot_sims = function(NREPS=50){
                       if (is.null(self$model)){
                         stop('Error in plot_sims(): need to fit model first')
                       }
                         plot(ecdf(self$annual_data$IPO_annual),
                              main=paste(NREPS,'replicates of IPO from ARMA model'))
                         nipo = length(self$annual_data$IPO_annual)
                         ppos = (1:nipo-0.5)/nipo
                         for (i in 1:NREPS){
                           sim = self$fsim(n=nipo)
                           lines(sort(sim),ppos, col='red')
                         }
                         plot(ecdf(self$annual_data$IPO_annual),add=T)
                         
                       #}
                       
                     
                    
                     
                       
                     }
                   )
)

ENSO_Obs <- R6Class("ENSO_Obs",
                    public = list(
                      annual_data = NULL,
                      year_filter = NULL,
                      model_nointercept = NULL,
                      model_intercept = NULL,
                      obs_ipo_ts = NULL,
                      fitted_dates = NULL, # Window of IPO/ENSO overlap
                      initialize = function(fpath=paste0(basepath, #="//fs1-cbr.nexus.csiro.au/{lw-hydroclimate}/work/
                                                         "pot063/SEA_stoch/new/soi.txt"),
                                            start_year = 1900,
                                            end_year = 2030,
                                            start_month = NULL){
                        ensowide=(read.table(fpath,skip=3,nrows=72,header=T))
                        names(ensowide) = c('YEAR', 1:12)
                        #browser()
                        if (!is.null(start_month)){
                          ensolong = ensowide %>% melt(id.vars=1)
                          cd=calendar2water(ensolong$YEAR, 
                                            as.numeric(as.character(ensolong$variable)),
                                            start_month=start_month)
                          ensolong$wyear = cd$wyear
                          ensolong$wmonth = cd$wmonth
                          
                          summed = tapply(X=ensolong$value, INDEX=ensolong$wyear,FUN=sum)
                          annual = data.frame(Year=as.numeric(row.names(summed)),
                                              ENSO_annual=as.vector(summed))
                        } else{
                          annual = data.frame(Year = ensowide$YEAR,
                                              ENSO_annual = apply(ensowide[,-1], 
                                                                  FUN=sum,
                                                                  1))
                        }
                        year_filter = c(start_year, end_year)
                        self$annual_data = annual %>% 
                          dplyr::filter(Year >= start_year & Year <= end_year)
                        
                      },
                      fit_model = function(IPO_Obj_Instance){
                        mgd = merge(self$annual_data,
                                    IPO_Obj_Instance$annual_data,
                                    by = 'Year')
                        self$fitted_dates = c(min(mgd$Year),
                                              max(mgd$Year))
                        
                        # self$model = auto.arima(x=mgd[,'ENSO_annual'], 
                        #                         xreg=mgd[,'IPO_annual'])
                        self$model_intercept = mean(mgd[,'ENSO_annual'])
                        self$model_nointercept = auto.arima(x=mgd[,'ENSO_annual'] - mean(mgd[,'ENSO_annual']), 
                                                            xreg=mgd[,'IPO_annual'])
                        self$obs_ipo_ts = mgd[,'IPO_annual']
                      },
                      sim = function(ipo_ts=NULL){ # pass in a vector only, realisation 
                                              # of an IPO replicate (or observation)
                                      # No argument? Uses observed IPO timeseries
                                      # that the model was fitted on.
                        if (is.null(self$model_nointercept)){
                          stop('Error in sim(): need to fit model first')
                        }
                        if (is.null(ipo_ts)){
                          ipo_ts = self$obs_ipo_ts
                        } 
                        # Fitted values (predictions), no intercept:
                        mgd = data.frame(IPO_annual=ipo_ts)
                        preds = as.vector(predict(self$model_nointercept, 
                                                  newxreg=mgd)$pred)
                        # add intercept:
                        preds = preds + self$model_intercept
                        # add in noise:
                        sd = sqrt(self$model_nointercept$sigma2)
                        sims = preds + rnorm(length(preds), 0, sd)
                        return(sims)
                      }
                    )
)
