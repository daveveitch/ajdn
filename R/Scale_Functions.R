#' Scale rule of thumbs
#'
#' @description Calculates rule of thumb values for \eqn{(\underline{s},\bar{s})}
#'
#' @param n length of time series
#' @param p dimension of time series
#' @return vector of length 2 \eqn{(\underline{s},\bar{s})}
#' @export
#' @examples
#' ajdn::scales_rot(1000,10)
scales_rot<-function(n,p){
  s_rot=rep(0,2)
  s_rot[1]=round((1/120)*n**(-1/3)*(log(p*n))**2*n)/n
  s_rot[2]=round((1/27)*n**(-1/6)*log(p*n)*n)/n

  if(s_rot[2]<s_rot[1]){s_rot[1]=s_rot[2]}

  return(s_rot)
}

#' Sparsify scales
#'
#' @description Calculates a sparse sequence of scales \eqn{s_{r,1},\dots,s_{r,\delta_n}}
#'
#' @param s_low lower scale \eqn{\underline{s}_r}
#' @param s_high upper scale \eqn{\bar{s}_r}
#' @param n length fo time series
#' @param p dimension of time series
#' @param max_scales OPTIONAL argument, maximum number of scales to return
#' @param epsilon default value = 0.51, used in calculation of \eqn{\delta_n}
#' @return vector of length \eqn{\delta_n} \eqn{s_{r,1},\dots,s_{r,\delta_n}}
#' @export
#' @examples
#' ajdn::sparsify_scales(0.05,0.08,1000,10)
sparsify_scales<-function(s_low,s_high,n,p,max_scales=NULL,epsilon=.51){

  if(!is.null(max_scales)){
    delta_n=max_scales
  }else{
    delta_n=ceiling((1/2000)*log(n)^(1+epsilon)*log(p*n)^(5/2))
  }

  sparse_scales=rep(0,delta_n)

  for(i in seq(1,delta_n)){
    # Calculate raw scale
    sparse_scales[i] = 2^( log(s_low,base=2) + (i-1)*(log(s_high,base=2)-log(s_low,base=2))/(delta_n-1))

    # Round the scale
    sparse_scales[i] = round(sparse_scales[i]*n)/n
  }

  # Get rid of duplicates
  sparse_scales=unique(sparse_scales)

  # In case where number of scales requested exactly equal to
  # number of possible ones
  if(!is.null(max_scales)){
    if((n*(s_high-s_low)+1)==max_scales){
      sparse_scales=seq(s_low*n,s_high*n)/n
    }
  }

  # case where no scales requested
  if(is.nan(sparse_scales[1])){
    sparse_scales=c(s_low,s_high)
  }

  return(sparse_scales)
}

#' Calculate LRV ratio
#'
#' @description This funtion calculates the desired LRV ratio based on n (see Appendix B.1 on selection of hyperparameters for further description).
#'
#' @param n length of time series
#' @return LRV ratio (a scalar)
#' @export
#' @examples
#' calc_desired_lrv_ratio(2000)
#
calc_desired_lrv_ratio<-function(n){
  bias=1+.85/log(n)
  return(bias)
}


#' Multiplier bootstrap long-run variance
#'
#' @description This function estimates the long run variance of the multiplier bootstrap
#'
#' @param s_prime_list vector of \eqn{s'} (block bootstrap sizes) that we are interested in
#' @param autocovariances vector of autocovariances, first entry is variance, will need a maximum autocovariance of ns'-1
#' @return a vector of long run variances of multiplier bootstrap corresponding to the different \eqn{s'} we tested
#' @export
#' @examples
#' print('enter example here')
bootstrap_lrv<-function(s_prime_list,autocovariances){
  lrv_list=rep(0,length(s_prime_list))

  # Fill in autocovariances that do not exist with 0
  if(length(s_prime_list)>length(autocovariances)){
    autocovariances=c(autocovariances,rep(0,length(s_prime_list)-length(autocovariances)))
  }

  count=1
  for(s_prime in s_prime_list){
    if(s_prime==1){
      lrv_list[count]=autocovariances[1]
    }else{
      summation=s_prime*(autocovariances[1])

      for(k in 1:(s_prime-1)){
        summation=summation+(s_prime-k)*(2*autocovariances[k+1])
      }
      lrv_list[count]=summation/s_prime
    }

    count=count+1
  }

  return(lrv_list)

}

#' Preprocesses scales
#'
#' @description Preprocesses scales into counts and \eqn{T_r} for each dimension
#'
#' @param s vector of scales (e.g. c(0.03,0.05,0.07)), or list of length p with vector of scales for each entry
#' @param p dimension of time series
#' @param n length of time series
#' @return a list of lists, and one vector,
#'                   1st entry - scale counts (scales in terms of number of observation,
#'                                     e.g. n=100 s=0.05, scale count=5),
#'                   2nd entry - t indexes (indexes where we test for jumps, \eqn{T_r} in paper)
#'                   NOTE - 1st and 2nd entry each contain p elements, corresponding to each dimension.
#'                   3rd entry - vector of unique scales across all dimensions. The reason we create this is to
#' @export
#' @examples
#' a = ajdn::preprocess_s(c(.02,.04),2,1000)
#' print(a)
preprocess_s<-function(s,p,n){
  # If only one vector of scales is given we turn this into a list where each
  # dimension has the same set of scales
  if(typeof(s)=="double"){
    scale_list=list()
    for(i in seq(1,p)){
      scale_list[[i]]=s
    }
  }else{
    if(length(scale_list!=p)){
        stop('Length of scale list must equal p')
    }
    scale_list=s
  }

  # Here we convert the scale_list into scale_count_list which is the number of
  # observations associated with the scale (this is easier to work with for computations).
  # Also we compile a vector of unique scales and their counts.
  scale_count_list=list()
  unique_scales=c()
  for(i in seq(1,p)){
    scale_count_list[[i]]=round(scale_list[[i]]*n)
    unique_scales=c(unique_scales,scale_list[[i]])
  }
  unique_scales=sort(unique(unique_scales))
  unique_scale_counts = unique(round(unique_scales*n))

  # Here we calculate the t_index_list which is a list of time indexes that
  # we test for each dimension this is based on the maximum scale ($T_r$ in paper for r=1,...,p)
  t_index_list=list()
  for(i in seq(1,p)){
    t_index_list[[i]]=seq(max(scale_count_list[[i]])+1,n-max(scale_count_list[[i]])-1,by=1)
  }

  return(list(scale_count_list,t_index_list,unique_scale_counts))
}
