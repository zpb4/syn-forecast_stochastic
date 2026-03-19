hyb_loess_fit<-function(obs,pred){
  ls_fit<-loess(pred~obs,span=0.75,degree=2,family='gaussian',control=loess.control(surface='interpolate',statistics = 'none',trace.hat = 'approximate'))
  #ls_fit<-loess(pred~obs,span=1,degree=1,family='symmetric',control=loess.control(surface='interpolate',statistics = 'none',trace.hat = 'approximate'))
  srt_obs<-sort(obs,index.return=T)
  ls_pred<-predict(ls_fit,srt_obs$x)
  diff_ls_pred<-diff(ls_pred)
  infl_idx<-length(diff_ls_pred[diff_ls_pred>=0])
  slope<-diff_ls_pred[infl_idx]/(srt_obs$x[infl_idx+1]-srt_obs$x[infl_idx])
  y<-ls_pred[infl_idx+1]
  intcpt<-y-slope*srt_obs$x[infl_idx+1]
  cutoff<-srt_obs$x[infl_idx+1]
  return(list(ls_fit,cutoff,slope,intcpt))
}

#this function uses the fitted parameters from hyb_loess_fit to predict conditional mean for new data
hyb_loess_out<-function(loess_mod,cutoff,slope,intcpt,x){
  xout<-sort(x,index.return=T)
  split_idx<-length(which(xout$x<=cutoff))
  if(split_idx==length(x)){
    yout<-predict(loess_mod,x)}
  if(split_idx<length(x)){
    y_loess<-predict(loess_mod,xout$x[1:split_idx])
    y_lm<-lin_mod(intcpt,slope,xout$x[(split_idx+1):length(x)])
    yout<-rep(0,length(x))
    yout[xout$ix[1:split_idx]]<-y_loess
    yout[xout$ix[(split_idx+1):length(x)]]<-y_lm}
  return(yout)
}

bcf_fun = function(lambda){
  bcf = exp(mean(lambda)-((sd(lambda)^2)/2))
}

#linear model helper functions
lin_mod<-function(intcpt,slope,x){
  y = intcpt + slope * x
  return(y)
}

lin_fit<-function(pars,x,y){
  pred = pars[1] + pars[2] * x
  err = sum((y - pred)^2)
  return(err)
}

log10_minor_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(minor_breaks))
  }
}

log10_major_break = function (...){
  function(x) {
    minx         = floor(min(log10(x), na.rm=T))-1;
    maxx         = ceiling(max(log10(x), na.rm=T))+1;
    n_major      = maxx-minx+1;
    major_breaks = seq(minx, maxx, by=1)
    minor_breaks = 
      rep(log10(seq(1, 9, by=1)), times = n_major)+
      rep(major_breaks, each = 9)
    return(10^(major_breaks))
  }
}