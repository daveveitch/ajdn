scales_rot<-function(n,p){
  # This function takes in a n,p and returns the high and low scales associated
  # with the rule of thumb
  # INPUT
  # n           - sample size
  # p           - dimension
  # OUTPUT
  # s_rot       - vector with two entries (s_low,s_high)

  s_rot=rep(0,2)
  s_rot[1]=round((1/120)*n**(-1/3)*(log(p*n))**2*n)/n
  s_rot[2]=round((1/27)*n**(-1/6)*log(p*n)*n)/n

  if(s_rot[2]<s_rot[1]){s_rot[1]=s_rot[2]}

  return(s_rot)
}

sparsify_scales<-function(s_low,s_high,n,p,max_scales=NULL){
  # This function sparsifies the number of scales for a predeterimned
  # number of scales
  # INPUT
  # s_low         - low scale
  # s_high        - high scale
  # n             - sample size
  # p             - dimension
  # max_scales    - optional argument - maximum number of scales to return
  # OUTPUT
  # sparse_scales - a vector of sparsified scales

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

calc_desired_lrv_ratio<-function(n){
  # This funtion calculates the desired LRV ratio based on n (see Appendix on selection
  # of hyperparameters for further description).
  # INPUT
  # n          -  sample size
  # OUTPUT
  # bias       - amount of bias that is ideal
  bias=1+.85/log(n)
  return(bias)
}

s_prime_select<-function(dgp,n,p,max_s_prime=NULL){
  # This function computes the level of s_prime to use based on the
  # desired LRV ratio, for various DGP. It does this by estimating the autocovariance structure
  # of each of the simulated error processes.
  # INPUT
  # dgp                  - name of dgp
  # n                    - sample size
  # p                    - dimension
  # OUTPUT
  # s_prime              - value of s_prime

  #if(is.null(max_s_prime)){max_s_prime=20}
  if(is.null(max_s_prime)){max_s_prime=round(n**(1/3))}

  desired_bias=calc_desired_lrv_ratio(n)

  # In the below if statements we estimate the autocovariance structure
  # of the simulated data, and then approximate its long run variance and compare
  # it to an approximation of the long run variance of the bootstrap.

  if(dgp%in%c('GS','GST')){
    s_prime_list=seq(1,max_s_prime)
    phi=.25
    th_div_bs_lrv_vec=rep(0,length(s_prime_list))

    for(i in 1:length(s_prime_list)){
      s_prime=s_prime_list[i]
      autocov <- ts.extend::ARMA.autocov(n = 2*(max_s_prime), ar = c(phi))
      LRV=autocov[1]+2*sum(autocov[2:(max_s_prime)])
      bootstrap_LRV=bootstrap_lrv(s_prime,autocov)
      th_div_bs_lrv_vec[i]=LRV/bootstrap_LRV
    }

    s_prime_unrounded=s_prime_list[which.min(abs(th_div_bs_lrv_vec-desired_bias))]
    s_prime=s_prime_unrounded/n
  }else if(dgp%in%c('PS','PST')){
    s_prime_list=seq(1,max_s_prime)

    ma_coef=c(0.5,0,0.5)
    th_div_bs_lrv_vec=rep(0,length(s_prime_list))

    for(i in 1:length(s_prime_list)){
      s_prime=s_prime_list[i]
      autocov <- ts.extend::ARMA.autocov(n = 2*(max_s_prime), ma = ma_coef)
      LRV=autocov[1]+2*sum(autocov[2:(max_s_prime)])
      bootstrap_LRV=bootstrap_lrv(s_prime,autocov)

      th_div_bs_lrv_vec[i]=LRV/bootstrap_LRV
    }

    # Select s_prime level
    s_prime_unrounded=s_prime_list[which.min(abs(th_div_bs_lrv_vec-desired_bias))]
    s_prime=s_prime_unrounded/n
  }else if(dgp%in%c('PLS','PLST')){
    coef_structure='x'
    ar_coefs=.25
    rho=.5
    sin_cycles=1

    # create coef_matrix
    coef_matrix=matrix(data=0,nrow=p,ncol=p)
    diag(coef_matrix)=rep(ar_coefs,p)
    if(coef_structure=='x'){
      for(i in 1:p){
        coef_matrix[p-i+1,i]=ar_coefs
      }
    }else if(coef_structure=='line'){
      for(i in 1:(p-1)){
        coef_matrix[i,i+1]=ar_coefs
      }
    }

    # Create the covariance matrix
    cov_matrix=matrix(data=0,nrow=p,ncol=p)
    for(ii in 1:p){
      for(jj in 1:p){
        cov_matrix[ii,jj]=rho**(abs(ii-jj))
      }
    }
    eigen_decomp=eigen(cov_matrix)
    cov_matrix_sqrt=eigen_decomp$vectors%*%diag(sqrt(eigen_decomp$values))

    # Create the sin coefficients that will be used
    sin_coefs=(sin(((1:n)/n)*2*pi*sin_cycles)+1)/2

    # Calculate autocovariances
    autocovariances=var_1_LRV(coef_matrix,cov_matrix,max_s_prime)
    true_lrv_vector=rep(0,p)
    for(jj in 1:p){
      true_lrv_vector[jj]=autocovariances[jj,jj,1]+2*sum(autocovariances[jj,jj,2:(dim(autocovariances)[3]-1)])
    }

    s_prime_list=seq(1,max_s_prime)
    s_prime_matrix=matrix(data=0,nrow=p,ncol=length(s_prime_list))

    bs_lrv_vec=rep(0,length(s_prime_list))

    for(ii in 1:p){
      for(jj in 1:length(s_prime_list)){
        s_prime=s_prime_list[jj]
        s_prime_matrix[ii,jj]=bootstrap_lrv(s_prime,autocovariances[ii,ii,])
      }
    }

    th_div_bs_lrv=true_lrv_vector/s_prime_matrix

    # Choose level of s' such that all dimensions have bias <= desired bias
    s_prime_unrounded=max(s_prime_list)
    for(ii in 1:length(s_prime_list)){
      if(sum(th_div_bs_lrv[,ii]<=desired_bias)==p){
        s_prime_unrounded=s_prime_list[ii]
        break
      }
    }

    s_prime=s_prime_unrounded/n
  }else if(dgp%in%c('IID','IIDT')){
    s_prime_unrounded=1
    s_prime=s_prime_unrounded/n
  }else if(dgp%in%c('LS','LST')){
    rho=.5
    ar_coefs=.25
    coef_structure='x'
    sin_cycles=1

    # create coef_matrix
    coef_matrix=matrix(data=0,nrow=p,ncol=p)
    diag(coef_matrix)=rep(ar_coefs,p)
    if(coef_structure=='x'){
      for(i in 1:p){
        coef_matrix[p-i+1,i]=ar_coefs
      }
    }else if(coef_structure=='line'){
      for(i in 1:(p-1)){
        coef_matrix[i,i+1]=ar_coefs
      }
    }

    # Create the covariance matrix
    cov_matrix=matrix(data=0,nrow=p,ncol=p)
    for(ii in 1:p){
      for(jj in 1:p){
        cov_matrix[ii,jj]=rho**(abs(ii-jj))
      }
    }
    eigen_decomp=eigen(cov_matrix)
    cov_matrix_sqrt=eigen_decomp$vectors%*%diag(sqrt(eigen_decomp$values))

    # Create the sin coefficients that will be used
    sin_coefs=(sin(((1:n)/n)*2*pi*sin_cycles)+1)/2

    # Calculate autocovariances
    autocovariances=var_1_LRV(coef_matrix,cov_matrix,max_s_prime)
    true_lrv_vector=rep(0,p)
    for(jj in 1:p){
      true_lrv_vector[jj]=autocovariances[jj,jj,1]+2*sum(autocovariances[jj,jj,2:(dim(autocovariances)[3]-1)])
    }

    s_prime_list=seq(1,max_s_prime)
    s_prime_matrix=matrix(data=0,nrow=p,ncol=length(s_prime_list))

    bs_lrv_vec=rep(0,length(s_prime_list))

    for(ii in 1:p){
      for(jj in 1:length(s_prime_list)){
        s_prime=s_prime_list[jj]
        s_prime_matrix[ii,jj]=bootstrap_lrv(s_prime,autocovariances[ii,ii,])
      }
    }

    th_div_bs_lrv=true_lrv_vector/s_prime_matrix

    # Choose level of s' such that all dimensions have bias <= desired bias
    s_prime_unrounded=max(s_prime_list)
    for(ii in 1:length(s_prime_list)){
      if(sum(th_div_bs_lrv[,ii]<=desired_bias)==p){
        s_prime_unrounded=s_prime_list[ii]
        break
      }
    }

    # Select s_prime level
    s_prime=s_prime_unrounded/n

  }

  return(s_prime)

}

bootstrap_lrv<-function(s_prime_list,autocovariances){
  # This function estimates the long run variance of the multiplier bootstrap
  # INPUT
  # s_prime_list    - list of s_primes (block bootstrap sizes) we are interested in
  # autocovariances - list of autocovariances, first one is variance, will need
  #                   a maximum autocovariance of ns'-1
  # OUTPU
  # lrv_list        - a list of long run variances

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

var_1_LRV<-function(coef_matrix,noise_cov,max_lag){
  # This function calculates the long run variances of a VAR_1 process using
  # a truncated sum
  # INPUTS
  # coef_matrix        - a pxp matrix of VAR coefficients
  # noise_cov          - covariance matrix of the noise
  # OUTPUTS
  # lrv_matrix         - a px21 matrix where each entry is equal to the diagonal
  #                      of the autocovariance matrix at lag i, i=1 is lag 0, i=21
  #                      is lag 20

  p=nrow(coef_matrix)

  # Psi matrices in causal representation X_t=\sum_j=0^\infty \Psi^j Z_{t-j}
  # which are equal to powers of coefficient matrices
  psi_array=array(data=0,dim=c(p,p,2*max_lag+1))
  psi_array[,,1]=diag(rep(1,p))
  psi_array[,,2]=coef_matrix
  for(j in 3:(2*max_lag)){
    psi_array[,,j]=psi_array[,,j-1]%*%coef_matrix
  }

  autocov_array=array(data=0,dim=c(p,p,max_lag+1))

  for(lag in 1:(max_lag+1)){
    autocov_array[,,lag]=psi_array[,,1+lag-1]%*%noise_cov%*%t(psi_array[,,1])

    for(j in 1:max_lag){
      autocov_array[,,lag]=autocov_array[,,lag]+psi_array[,,1+lag-1+j]%*%noise_cov%*%t(psi_array[,,1+j])
    }

  }

  return(autocov_array)

}

preprocess_s<-function(s,p,n){
  # This function takes the values of s provided to the AJDN function and preprocesses them
  # into counts, and also calcualtes $T_r$ for each dimension (time indexes to test for jumps)
  # INPUTS
  # s               - vector of scales (in 1/n form), or list of length p with vector of scales for each entry
  # p               - dimension
  # n               - time series length
  # OUTPUTS
  # preprocess_list - a list of lists, and one vector,
  #                   [[1]] scale counts (scales in terms of number of observation,
  #                                      e.g. n=100 s=0.05, scale count=5),
  #                   [[2]] t indexes (indexes where we test for jumps, $T_r$ in paper)
  #                   NOTE - preprocess_list[[1]] and [[2]] each contain p elements, corresponding to each dimension.
  #                   [[3]] vector of unique scales across all dimensions. The reason we create this is to
  #

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
