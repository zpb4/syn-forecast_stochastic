#Forecast Verification Functions
library(lubridate)

wy_fun<-function(date_vec){
  wy_vec <- date_vec$year
  wy_vec[date_vec$mo%in%c(9,10,11)] <- wy_vec[date_vec$mo%in%c(9,10,11)]+1
  date_vec_wy <- date_vec
  date_vec_wy$year <- wy_vec
  return(date_vec_wy)
}

fwd_forecast_rearrange<-function(forecast){
  forecast_out<-array(0,dim(forecast))
  for(i in 1:dim(forecast)[3]){
    forecast_out[,(i+1):dim(forecast)[2],i]<-forecast[,1:(dim(forecast)[2]-i),i]
  }
  return(forecast_out)
}

ens_rank<-function(ens_forc,obs){
  f_sort<-sort(ens_forc)
  dif<-f_sort - obs
  if(all(dif<0)==T){rnk<-length(dif)+1} else {y<-min(dif[dif>=0]);rnk<-which(dif==y)}
  if(length(rnk)>1){rnk<-sample(rnk,1)}
  return(rnk)
}

roll_sum<-function(x){out<-c();out[1]<-x[1];
for(i in 1:(length(x)-1)){out[i+1]<-out[i]+x[i+1]};return(out)}

chi2<-function(freqs){
  m<-length(freqs) #already 'm + 1' in current formulation
  n<-sum(freqs)
  chi2<-((m)/n)*sum((freqs - (n / m))^2)
  prob<-pgamma(chi2,shape=((m-1)/2),scale=2,lower.tail=F)
  return(list(chi2,prob))
}

chi2_comp<-function(freq_ref,freq_tst){
  freq_ref[freq_ref<1]<-1
  m<-length(freq_ref) 
  n<-sum(freq_ref)
  chi2<-sum(((freq_tst - freq_ref)^2)/freq_ref)
  prob<-pgamma(chi2,shape=((m-1)/2),scale=2,lower.tail=F)
  return(list(chi2,prob))
}

RI<-function(freqs){
  m<-length(freqs) #already 'm + 1' in current formulation
  n<-sum(freqs)
  ri<-(1/n)*sum(abs(freqs - (n / (m+1))))
  return(ri)
}

Ent1<-function(freqs){
  freqs<-freqs[freqs>0]
  m<-length(freqs) #already 'm + 1' in current formulation
  n<-sum(freqs)
  ent<-((-1)/log(m+1))*sum((freqs/n)*log(freqs/n))
  return(ent)
}

Ent<-function(freqs){
  freqs[freqs<1]<-1
  m<-length(freqs) #already 'm + 1' in current formulation
  n<-sum(freqs)
  ent<-((-1)/log(m+1))*sum((freqs/n)*log(freqs/n))
  return(ent)
}

mn_ens<-function(ens){
  x_bar_t<-(1/length(ens))*sum(ens)
  s2<-(1/length(ens))*sum((ens - x_bar_t)^2)
  return(x_bar_t)
}

st2_ens<-function(ens){
  x_bar_t<-(1/length(ens))*sum(ens)
  s2<-(1/length(ens))*sum((ens - x_bar_t)^2)
  return(s2)
}

bin_spread<-function(ens_vars,ens_num){
  sprd<-(ens_num/(ens_num-1))*(1/length(ens_vars))*sum(ens_vars)
  return(sqrt(sprd))
}

bin_mse<-function(ens_mns,obs,ens_num){
  mse<-(ens_num/(ens_num+1))*(1/length(ens_mns))*sum((ens_mns - obs)^2)
  return(sqrt(mse))
}

eCRPS<-function(ens,obs){
  m<-length(ens)
  t1<-(1/m)*sum(abs(ens-obs))
  #t2<-array(0,c(m,m))
  #for(i in 1:(m-1)){
    #for(j in (i+1):m){
      #t2[i,j]<-abs(ens[i]-ens[j])
    #}
  #}
  
  mat1 = matrix(rep(ens,each=length(ens)),ncol=length(ens),byrow=F)
  mat2 = matrix(rep(ens,length(ens)),ncol=length(ens),byrow=F)
  t2 = abs(mat1[lower.tri(mat1)]-mat2[lower.tri(mat2)])
  
  t2_res<-(1 / (m*(m-1)))*sum(t2)
  return(t1 - t2_res)
}

eDSS<-function(var_ens,mn_ens,obs){
  ss<-log(var_ens) + (((obs - mn_ens)^2)/var_ens^2)
  return(ss)
}

declust_evts_extract <- function(Q,n_evts,sep,max_lds){
  rnk_data = rank(Q,ties.method = 'first')
  srt_rnks_idx = order(rnk_data,decreasing = T)[1:(n_evts*10)]

  vec = c()
  for(i in 1:length(srt_rnks_idx)){
    evt <- which(rnk_data==max(rnk_data[max(srt_rnks_idx[i]-sep,1):min(srt_rnks_idx[i]+sep,length(rnk_data))]))
    vec[i] <- evt
    if(i>1 & evt%in%(vec[1:(i-1)]+matrix(rep(-sep:sep,length(vec[1:(i-1)])),nrow=length(vec[1:(i-1)]),byrow=T))){vec[i]<-NA}}
  declust_evts=unique(vec)
  sel_idx=order(Q[declust_evts],decreasing = T)
  evt_idx = declust_evts[sel_idx][0:n_evts]
  
  rmv_idx <- which(evt_idx<=leads)
  if(length(rmv_idx)>0){
    pool <- declust_evts[sel_idx][(n_evts+1):(length(declust_evts[sel_idx]))]
    pool <- pool[pool>leads]
    evt_idx <- evt_idx[-c(rmv_idx)]
    evt_idx <- c(evt_idx,pool[1:length(rmv_idx)])}
  
  return(evt_idx)}

declust_evts_extract_xhc <- function(Q,n_evts,sep,max_lds,idx_86,idx_hefs,idx_all,has_86){
  ixx_all = as.character(idx_all)
  rnk_data = rank(Q,ties.method = 'first')
  srt_rnks_idx = order(rnk_data,decreasing = T)[1:(n_evts*10)]
  
  vec = c()
  for(i in 1:length(srt_rnks_idx)){
    evt <- which(rnk_data==max(rnk_data[max(srt_rnks_idx[i]-sep,1):min(srt_rnks_idx[i]+sep,length(rnk_data))]))
    vec[i] <- evt
    if(i>1 & evt%in%(vec[1:(i-1)]+matrix(rep(-sep:sep,length(vec[1:(i-1)])),nrow=length(vec[1:(i-1)]),byrow=T))){vec[i]<-NA}}
  declust_evts=unique(vec)
  sel_idx=order(Q[declust_evts],decreasing = T)
  evt_idx = declust_evts[sel_idx][0:n_evts]
  obs_dates = ixx_all[evt_idx]
  
  rmv_idx <- which(evt_idx<=leads)
  if(has_86==T){
    idx_close_86 <- seq(as.Date(tail(idx_86,1)+(60*60*24)),as.Date(tail(idx_86,1)+(60*60*24*max_lds)),by='day')
    if(any(obs_dates%in%idx_close_86==T)){
      rmv_idx <- c(rmv_idx,(1:length(obs_dates))[obs_dates%in%idx_close_86])}
  }
  #remove any dates that are within the lead window of the end of the HEFS period to ensure no cutoffs
  idx_close_hefsend <- seq(as.Date(tail(idx_hefs,1)+(60*60*24)),as.Date(tail(idx_hefs,1)+(60*60*24*max_lds)),by='day')
  if(any(obs_dates%in%idx_close_hefsend==T)){
    rmv_idx <- c(rmv_idx,(1:length(obs_dates))[obs_dates%in%idx_close_hefsend])}
  if(length(rmv_idx)>0){
    pool <- declust_evts[sel_idx][(n_evts+1):(length(declust_evts[sel_idx]))]
    pool <- pool[pool>leads]
    evt_idx <- evt_idx[-c(rmv_idx)]
    evt_idx <- c(evt_idx,pool[1:length(rmv_idx)])}
  
  return(evt_idx)}

wy_fun<-function(date_vec){
  wy_vec <- date_vec$year
  wy_vec[date_vec$mo%in%c(9,10,11)] <- wy_vec[date_vec$mo%in%c(9,10,11)]+1
  date_vec_wy <- date_vec
  date_vec_wy$year <- wy_vec
  return(date_vec_wy)
}


climo_forecast <- function(dtg_vec,hefs_array,obs_forward_array){
  yrs <- unique(dtg_vec$year)
  if(length(which(dtg_vec$year==yrs[1]))<365){
    yrs <- yrs[-c(1)]
  }
  if(length(which(dtg_vec$year==yrs[length(yrs)]))<365){
    yrs <- yrs[-c(length(yrs))]
  }
  yrs_idx <- unique(dtg_vec$year)
  climo_farray <- array(NA,dim(hefs_array))
  doy <- yday(dtg_vec)

  for(i in 1:length(yrs_idx)){
    set.seed(1)
    n_ens <- dim(hefs_array)[1]
    if(n_ens<=(length(yrs)-1)){
      climo_yrs <- sample(yrs[!yrs%in%yrs_idx[i]],n_ens,replace = F)}
    if(n_ens>(length(yrs)-1)){
      climo_yrs <- sample(yrs[!yrs%in%yrs_idx[i]],n_ens,replace = T)}
    yr_len <- length(which(dtg_vec$year==yrs_idx[i]))
    climo_inp <- climo_farray[,1:yr_len,]
    for(j in 1:yr_len){
      jday <- min(j,365) #for leap years, just repeat the end of year climatology forecast because every year does not have a day #366
      #for the particular case where the beginning year is also a leap year, don't use jday index 366
      if(i==1 & doy[jday]==366){
        jday = jday-1
      }
      dy <- doy[dtg_vec$year%in%yrs_idx[i]][jday]
      doy_idx <- rep(NA,length(dtg_vec))
      doy_idx[dtg_vec$year%in%climo_yrs] <- doy[dtg_vec$year%in%climo_yrs]
      d_idx <- c()
      for(k in 1:length(climo_yrs)){
        d_idx[k] = which(doy_idx==dy & dtg_vec$year==climo_yrs[k])
      }
      climo_inp[,j,] <- obs_forward_array[d_idx,]
    }
    climo_farray[,dtg_vec$year%in%yrs_idx[i],] <- climo_inp
  }
  return(climo_farray)
}

#dtg_vec = ixx_gen_swm
#hefs_array = syn_hefs_forward_swm[1,1,,,]
#obs_forward_array = obs_forward_all_leads_swm[1,1,,2:(leads+1)]

climo_forecast_swm <- function(dtg_vec,hefs_array,obs_forward_array){
  yrs <- unique(dtg_vec$year)
  yrs_idx <- yrs[-c(1,length(yrs))]
  dtg_vec_samp <- dtg_vec[dtg_vec$year%in%yrs_idx]
  
  climo_farray <- array(NA,dim(hefs_array))
  doy <- yday(dtg_vec)
  doy_samp <- yday(dtg_vec_samp)
  
  obs_fwd_samp <- obs_forward_array[dtg_vec$year%in%yrs_idx,]
  obs_fwd_noleap <- obs_fwd_samp[!doy_samp%in%366,]
  dtg_vec_noleap <- dtg_vec_samp[!doy_samp%in%366]
  
  for(i in 1:length(yrs_idx)){
    set.seed(1)
    yr_len <- length(which(dtg_vec_samp$year==yrs_idx[i]))
    climo_inp <- climo_farray[,1:yr_len,]
    if(n_ens<=(length(yrs)-1)){
      climo_yrs <- sample(yrs_idx[!yrs_idx%in%yrs_idx[i]],n_ens,replace = F)
      inp_samp <- obs_fwd_noleap[dtg_vec_noleap$year%in%climo_yrs,]}
    if(n_ens>(length(yrs)-1)){
      climo_yrs <- sample(yrs_idx[!yrs_idx%in%yrs_idx[i]],n_ens,replace = T)
      inp_samp <- array(NA,c(365*n_ens,leads))
      for(k in 1:length(climo_yrs)){
        inp_samp[((k-1)*365+1):(k*365),] <- obs_fwd_noleap[dtg_vec_noleap$year%in%climo_yrs[k],]
      }}
    
    #no leap year sampling produces 365 day samples only
    if(yr_len==365){
      for(j in 1:365){
        climo_inp[,j,] <- inp_samp[seq(j,365*n_ens,365),]}}
    #if leap year, recycle the day 365 sample
    if(yr_len==366){
      for(j in 1:365){
        climo_inp[,j,] <- inp_samp[seq(j,365*n_ens,365),]}
      climo_inp[,366,] <- climo_inp[,365,]}
    inp_idx <- dtg_vec$year%in%yrs_idx[i]
    climo_farray[,inp_idx,] <- climo_inp
  }
  #add data back to first year of dataset by pulling from second  year
  init_yr_inp <- climo_farray[,dtg_vec$year%in%yrs_idx[1],][,doy[dtg_vec$year%in%yrs[1]],]
  last_yr_inp <- climo_farray[,dtg_vec$year%in%yrs_idx[length(yrs_idx)],][,doy[dtg_vec$year%in%yrs[length(yrs)]],]
  climo_farray[,dtg_vec$year%in%yrs[1],] <- init_yr_inp
  climo_farray[,dtg_vec$year%in%yrs[length(yrs)],] <- last_yr_inp

  return(climo_farray)
}

skill_collapse_cumul <- function(tgt_forc,base_forc,mod_vec){
  lds = dim(base_forc)[2]
  cumul_base <- t(apply(base_forc,1,roll_sum))
  cumul_tgt <- t(apply(tgt_forc,1,roll_sum))
  vec_base_c <- sort(cumul_base[,lds],index.return=T)
  vec_tgt_c <- sort(cumul_tgt[,lds],index.return=T)
  cumul_out <- array(NA,dim(cumul_base))
  cumul_out[vec_base_c$ix,] <- cumul_base[vec_base_c$ix,] + (cumul_tgt[vec_tgt_c$ix,] - cumul_base[vec_base_c$ix,]) * matrix(rep(mod_vec,n_ens),ncol=lds,byrow=T)
  base_mod <- cbind(cumul_out[,1],t(apply(cumul_out,1,diff)))
  return(base_mod)
}

##################################END###################################################