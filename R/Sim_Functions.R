#' Autocov of VAR(1) model
#'
#' @description This function calculates the autocovariance array of a VAR_1 process using a truncated sum.
#'
#' @param coef_matrix \eqn{p \times p} matrix of VAR coefficients
#' @param noise_cov covariance matrix of the noise
#' @param max_lag maximum lag to return of the autocovariance
#' @return array of autocovariance matrices
#' @export
#' @examples
#' autocov_array = ajdn::var_1_LRV(matrix(c(0.25,0.1,0.25,0.1),nrow=2),diag(c(1,1)),20)
#' print(autocov_array[,,1])
var_1_LRV<-function(coef_matrix,noise_cov,max_lag){

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


#' Choose \eqn{s'}
#'
#' @description This function computes \eqn{s'} based on the desired LRV ratio from Appendix B.1.2, for data generating processes used in the simulations.
#'
#' @param dgp data generating process, one of: IID, IIDT, GS, GST, PS, PST, LS, LST, PLS, PLST
#' @param n length of time series
#' @param p dimension of time series
#' @param max_s_prime OPTIONAL argument, the maximium \eqn{ns'} to consider (default is \eqn{n^{1/3}})
#' @return value of \eqn{s'}
#' @export
#' @examples
#' ajdn::s_prime_select('GS',1000,100)
 s_prime_select<-function(dgp,n,p,max_s_prime=NULL){

  if(is.null(max_s_prime)){max_s_prime=round(n**(1/3))}

  desired_bias=calc_desired_lrv_ratio(n)

  # In the below if statements we calculate the autocovariance structure
  # of the simulated data, and then approximate its long run variance and compare
  # it to an approximation of the long run variance of the bootstrap.

  if(dgp%in%c('GS','GST')){
    s_prime_list=seq(1,max_s_prime)
    phi=.25
    th_div_bs_lrv_vec=rep(0,length(s_prime_list))

    # Autocovariance for AR(1) process
    autocov = rep(0,2*max_s_prime)
    autocov[1] = 1/(1-phi**2) # variance ar(1)
    for(l in 2:length(autocov)){
      autocov[l] = autocov[1]*phi**(l-1)
    }

    for(i in 1:length(s_prime_list)){
      s_prime=s_prime_list[i]
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

    # Autocovariance MA(3) process with coefficients c(0.5,0,0.5)
    autocov = rep(0,2*max_s_prime)
    autocov[1] = 1.5
    true_autocov = c(1.5,0.5,0.25,0.5)

    for(l in 2:length(autocov)){
      if(l <= length(true_autocov)){
        autocov[l] = true_autocov[l]
      }
    }

    for(i in 1:length(s_prime_list)){
      s_prime=s_prime_list[i]
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


#' Generate simulated data
#'
#' @description This function generates simulated time series as used in the paper.
#'
#' @param dgp data generating process, one of: IID, IIDT, GS, GST, PS, PST, LS, LST, PLS, PLST
#' @param n length of time series
#' @param p dimension of time series
#' @param s either vector of scales (if same for all dimension), or list of scales if different across dimensions used for testing
#' @param x_matrix_split OPTIONAL argument, if true return list with first element the nosie matrix, and second
#' argument the potentially time varying mean matrix
#' @return either a \eqn{p \times n} matrix of the time series, or a list of the noise matrix and mean matrix
#' @export
#' @examples
#' ex_ts = ajdn::generate_x_matrix(1000,100,'GST',c(.05,.06))
generate_x_matrix<-function(n,p,dgp,s,x_matrix_split=NULL){

   mean_matrix=matrix(data=rep(0,n*p),nrow=p,ncol=n)

   # Generate X matrix according to data generating processes
   if(dgp=='IID'){
     # IID normal, independent dimensions
     noise_matrix = matrix(data=stats::rnorm(n*p),ncol=n)
   }else if(dgp=='IIDT'){
     noise_matrix = generate_x_matrix(n,p,'IID',s)
     mean_matrix =  generate_means(noise_matrix,s)
   }else if(dgp=='GS'){
     # AR(1) processes, independent dimensions
     noise_matrix = matrix(data=rep(0,n*p),nrow=p,ncol=n)
     for(i in 1:p){
       innovations=stats::rnorm(n+1)
       noise_matrix[i,1]=innovations[2]+0.25*innovations[1]
       for(j in 2:n){
         noise_matrix[i,j]=innovations[j+1]+0.25*noise_matrix[i,j-1]
       }
     }
   }else if(dgp=='GST'){
     noise_matrix = generate_x_matrix(n,p,'GS',s)
     mean_matrix =  generate_means(noise_matrix,s)
   }else if(dgp=='PS'){
     # VMA(3) process with jump in variance
     sd_jump=2
     rho=.5
     ma_coef=c(.5,0,.5)
     q=length(ma_coef)

     noise_matrix=matrix(data=0,nrow=p,ncol=n)

     # Here create an equicorrelation matrix for first half, and then one
     # with twice the variance in second half
     cov_matrix_1=matrix(data=rho,nrow=p,ncol=p)
     diag(cov_matrix_1)=1
     cov_matrix_2=diag(rep(sd_jump,p))%*%cov_matrix_1%*%diag(rep(sd_jump,p))

     eigen_1=eigen(cov_matrix_1)
     eigen_2=eigen(cov_matrix_2)

     cov_matrix_1_sqrt=eigen_1$vectors%*%diag(sqrt(eigen_1$values))
     cov_matrix_2_sqrt=eigen_2$vectors%*%diag(sqrt(eigen_2$values))

     innovations=matrix(data=stats::runif((n+q)*p,min=-1*sqrt(12)/2,max=sqrt(12)/2),nrow=p,ncol=(n+q))
     innovations[,1:(round(n/2)+q)]=cov_matrix_1_sqrt%*%innovations[,1:(round(n/2)+q)]
     innovations[,(round(n/2)+q+1):(n+q)]=cov_matrix_2_sqrt%*%innovations[,(round(n/2)+q+1):(n+q)]

     for(l in 1:n){
       noise_matrix[,l]=innovations[,l+q]
       for(m in 1:q){
         noise_matrix[,l]=noise_matrix[,l]+innovations[,l+q-m]*ma_coef[m]
       }
     }

   }else if(dgp=='PST'){
     noise_matrix = generate_x_matrix(n,p,'PS',s)
     mean_matrix =  generate_means(noise_matrix,s)
   }else if(dgp=='LS'){
     # VAR(1) model with time varying coefficients
     rho=.5
     ar_coefs=.25
     coef_structure='x'
     sin_cycles=1

     # create coefficient matrix
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

     # Create the covariance matrix - Kac-Murdock-Szego
     cov_matrix=matrix(data=0,nrow=p,ncol=p)
     for(ii in 1:p){
       for(jj in 1:p){
         cov_matrix[ii,jj]=rho**(abs(ii-jj))
       }
     }
     eigen_decomp=eigen(cov_matrix)
     cov_matrix_sqrt=eigen_decomp$vectors%*%diag(sqrt(eigen_decomp$values))

     # Create the sin coefficients that will be used to make the coefficient matrix time varying
     sin_coefs=(sin(((1:n)/n)*2*pi*sin_cycles)+1)/2

     noise_matrix=matrix(data=0,nrow=p,ncol=n)

     # Create a matrix of innovations
     innovations_matrix=matrix(data=stats::rbinom((n+1)*p,10,.3)-3,nrow=p,ncol=n+1)
     x_0=cov_matrix_sqrt%*%innovations_matrix[,1]
     noise_matrix[,1]=(coef_matrix*sin_coefs[1])%*%x_0+cov_matrix_sqrt%*%innovations_matrix[,2]

     for(ii in 2:n){
       noise_matrix[,ii]=(coef_matrix*sin_coefs[ii])%*%noise_matrix[,ii-1]+cov_matrix_sqrt%*%innovations_matrix[,ii+1]
     }
   }else if(dgp=='LST'){
     noise_matrix = generate_x_matrix(n,p,'LS',s)
     mean_matrix =  generate_means(noise_matrix,s)
   }else if(dgp=='PLS'){
     # Same as dgp=='LS' with exception there is a jump in variance of innovations in middle
     # of time series
     sd_jump_number=2
     sd_jump=2
     rho=.5
     ar_coefs=.25
     coef_structure='x'
     sin_cycles=1

     # create a sd_vector with each entry equal to the sd of innovations at a given time
     sd_vector=rep(1,n)
     sd_segment_length=round(n/(sd_jump_number+1))+1

     count_n=1
     count_sd=1

     while(count_n<=n){
       if(count_sd%%2==0){
         sd_vector[count_n:(count_n+sd_segment_length-1)]=sd_jump
         count_sd=1
       }else{
         sd_vector[count_n:(count_n+sd_segment_length-1)]=1
         count_sd=2
       }
       count_n=count_n+sd_segment_length
     }

     # Add an extra 1 because the innovations vector has one extra
     # innovation for burn in
     sd_vector=c(1,sd_vector[1:n])

     # create coefficient
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

     noise_matrix=matrix(data=0,nrow=p,ncol=n)

     # Create an innovations matrix
     innovations_matrix=matrix(data=stats::rbinom((n+1)*p,10,.3)-3,nrow=p,ncol=n+1)

     for(ii in 1:p){
       innovations_matrix[ii,]=innovations_matrix[ii,]*sd_vector
     }

     x_0=cov_matrix_sqrt%*%innovations_matrix[,1]
     noise_matrix[,1]=(coef_matrix*sin_coefs[1])%*%x_0+cov_matrix_sqrt%*%innovations_matrix[,2]

     for(ii in 2:n){
       noise_matrix[,ii]=(coef_matrix*sin_coefs[ii])%*%noise_matrix[,ii-1]+cov_matrix_sqrt%*%innovations_matrix[,ii+1]
     }

   }else if(dgp=='PLST'){
     noise_matrix = generate_x_matrix(n,p,'PLS',s)
     mean_matrix =  generate_means(noise_matrix,s)
   }

   if(is.null(x_matrix_split)){
     return(noise_matrix+mean_matrix)
   }else{
     if(x_matrix_split){
       return(list(noise_matrix,mean_matrix))
     }else{
       return(noise_matrix+mean_matrix)
     }
   }

}

#' Generate time varying means
#'
#' @description This function generates sin curve means with slightly different phases in each dimension
#'
#' @param noise_matrix matrix of random noise being used
#' @param s either vector of scales (if same for all dimension), or list of scales if different across dimensions used for testing
#' @return either a \eqn{p \times n} matrix of the time varying mean
#' @export
#' @examples
#' noise_matrix = ajdn::generate_x_matrix(1000,100,'IID',c(.05,.06))
#' mean_matrix = ajdn::generate_means(noise_matrix,c(.05,.06))
generate_means<-function(noise_matrix,s){
   p=dim(noise_matrix)[1]
   n=dim(noise_matrix)[2]

   if(typeof(s)=="double"){
     scale_list=list()
     for(i in seq(1,p)){
       scale_list[[i]]=s
     }
   }else{
     scale_list=s
   }

   # Matrix with estimates of local standard deviation
   local_sd_matrix=calc_local_sd(noise_matrix,s)

   # Create matrix of the generated time varying mean, and then the time varying mean
   # adjusted for the local standard deviation at each time point
   mu_matrix=matrix(data=0,nrow=p,ncol=n)
   mu_matrix_scaled=matrix(data=0,nrow=p,ncol=n)

   for(i in seq(1,p)){
     mu_matrix[i,]=sin(((1:n)/n)*2*pi+(i/p)*2*pi)

     mu_matrix_scaled[i,1]=mu_matrix[i,1]

     for(j in seq(2,n)){
       mu_matrix_scaled[i,j]=mu_matrix_scaled[i,j-1]+(mu_matrix[i,j]-mu_matrix[i,j-1])*local_sd_matrix[i,j]
     }

     mu_matrix_scaled[i,] = mu_matrix_scaled[i,]

   }

   return(mu_matrix_scaled)
 }

#' Calculate local standard deviations
#'
#' @description This function calculates the local standard deviation matrix for a time series. This
# is used to normalize the means that are added to the data
#'
#' @param noise_matrix matrix of random noise being used
#' @param s either vector of scales (if same for all dimension), or list of scales if different across dimensions used for testing
#' @return \eqn{p \times n} matrix (same size as noise matrix) where each entry estimate of local standard deviation at that time
#' @export
#' @examples
#' noise_matrix = ajdn::generate_x_matrix(1000,100,'IID',c(.05,.06))
#' local_sd = ajdn::calc_local_sd(noise_matrix,c(.05,.06))
 calc_local_sd<-function(noise_matrix,s){
   # This function calculates the local standard deviation matrix for a time series. This
   # is used to normalize the test statistics the means that are added to the data
   # INPUT
   # noise_matrix               - matrix of time series
   # s                      - vector of scales

   # Output
   # local_sd_matrix        - matrix same size as noise_matrix where each entry esimate
   #                          of local standard deviation at each time
   p=dim(noise_matrix)[1]
   n=dim(noise_matrix)[2]

   if(typeof(s)=="double"){
     scale_list=list()
     for(i in seq(1,p)){
       scale_list[[i]]=s
     }
   }else{
     scale_list=s
   }

   scale_count_list=list()
   for(i in seq(1,p)){
     scale_count_list[[i]]=round(scale_list[[i]]*n)
   }

   # CALCULATE OPTIMAL VALUE OF S*
   # Define what s* is (here we set it by default to be the minimum scale of each dim)
   s_star_count_list=list()
   for(i in seq(1,p)){
     s_star_count_list[[i]]=min(scale_count_list[[i]])
   }

   # Here we calculate the t_index_list which is a list of time indexes that
   # we test for each dimension this is based on the maximum scale
   t_index_list=list()
   for(i in seq(1,p)){
     t_index_list[[i]]=seq(max(scale_count_list[[i]])+1,n-max(scale_count_list[[i]])-1,by=1)
   }

   # COMPUTE THE NORMALIZATION MATRIX
   # this is needed to normalize test/bootstrap statistics to account for
   # the time series' different varainces
   if(length(unique(scale_list))==1){
     normalization_df=time_series_normalization_fast(noise_matrix,t_index_list,scale_count_list,loud=NULL)
   }else{
     normalization_df = time_series_normalization(noise_matrix,t_index_list,scale_count_list)
   }

   full_normalization_matrix = matrix(data=NA,nrow=p,ncol=n)

   for(i in seq(1,p)){
     min_t_index=min(t_index_list[[i]])
     max_t_index=max(t_index_list[[i]])
     full_normalization_matrix[i,min_t_index:max_t_index]=normalization_df[((normalization_df$dimension==i)&
                                                                              (normalization_df$t_index%in%t_index_list[[i]])),'sd']

     # fill in first few, and last few observations
     full_normalization_matrix[i,1:min_t_index]=full_normalization_matrix[i,min_t_index]
     full_normalization_matrix[i,max_t_index:n]=full_normalization_matrix[i,max_t_index]
   }

   return(full_normalization_matrix)

 }

#' Calculates local standard deviation at poitn in time
#'
#' @description This function takes a matrix of time series and then calculates for a specific dimension and specific time
#' what the local standard deviation is the reason behind having this function is a quick way to calculate local
#' standard deviation for the purposes of adding changepoints of the correct scale to simulated data.
#'
#' @param x_matrix \eqn{p \times n} matrix of time series
#' @param s_low low scale for dimension adding jump
#' @param s_high high scale for dimension jumps is being added to
#' @param interested_dim which dimension we want local_sd for
#' @param interested_t which time (in (0,1)) we want local standard deviation for
#' @return local standard deviation
#' @export
#' @examples
#' x_matrix = ajdn::generate_x_matrix(1000,100,'IID',c(.05,.06))
#' sd = ajdn::single_local_sd(x_matrix,.05,.07,1,.5)
single_local_sd<-function(x_matrix,s_low,s_high,interested_dim,interested_t){
   n=dim(x_matrix)[2]

   t_index=round(interested_t*n)
   min_scale_count=round(s_low*n)
   max_scale_count=round(s_high*n)

   # Here if doing a single scale case, just use 1/2 the scale count as the window
   if(min_scale_count==max_scale_count){min_scale_count=round(max_scale_count/2)}

   window_width = max_scale_count-min_scale_count + 1 # +1 given we include the endpoints of the interval

   left_mean = mean(x_matrix[interested_dim,(t_index-max_scale_count):(t_index-min_scale_count)])
   right_mean = mean(x_matrix[interested_dim,(t_index+min_scale_count):(t_index+max_scale_count)])

   sample_sd = sqrt((sum((x_matrix[interested_dim,(t_index-max_scale_count):(t_index-min_scale_count)]-left_mean)**2)+
                       sum((x_matrix[interested_dim,(t_index+min_scale_count):(t_index+max_scale_count)]-right_mean)**2))/(2*window_width))

   return(sample_sd)
 }

#' Add scaled jump
#'
#' @description For the simulation scenarios used, add the scaled jumps to the time series
#'
#' @param scenario_number 1 : \eqn{\gamma \times p/2} dimensions undergo jump at t=.25, separate  \eqn{\gamma \times p/2} jump at t=0.75
#' 2: Evenly spaced jumps between t=0.2 t=-.8 for \eqn{\gamma} dimensions undergoing jumps
#' 10: first dimension undergoes a single jump at t=0.5
#' @param gamma proportion of dimensions undergoing a jump
#' @param signal_to_noise size of jump relative to local sd
#' @param x_matrix \eqn{p\times n} matrix of timme series
#' @param s vector of scales
#' @return list, first entry matrix which can be added to x_matrix to add jumps, second matrix is matrix with two columns 'dim', 't'
#' tracking where jumps occur.
#' @export
#' @examples
#' x_matrix = ajdn::generate_x_matrix(1000,100,'IID',c(.05,.06))
#' add_cp = ajdn::add_changepoints_scaled(1,0.5,1,x_matrix,c(.05,.06))
 add_changepoints_scaled<-function(scenario_number,gamma,signal_to_noise,x_matrix,s){
   # INPUTS
   # scenario_number        - the changepoint scenario we are adding CPs for
   # scenario_parameters    - proportion of dimensions undergoing a CP
   # signal_to_noise        - size of changepoint in relation to local standard deviation
   # x_matrix               - matrix of time series which we are applying CPs to
   # s                      - scales being considered
   # num_s_star_to_consider - how many s_star to use in calculation of local standard deviation matrix

   # OUTPUTS
   # jump_matrix            - a matrix which can be added to the x_matrix for changepoints
   # changepoint_locations  - a list of changepoint locations, note that a cp occuring at t=.5
   #                        - means that the mean increases at t=.5+1/n

   # CONSTANTS
   # the probability a jump is positive or negative is 0.5

   n=dim(x_matrix)[2]
   p=dim(x_matrix)[1]

   # This is a matrix containing the jumps which we will add to the time series
   jump_matrix = matrix(data=rep(0,n*p),nrow=p,ncol=n)

   # This is a matrix which will keep track of the locations of the changepoints will
   changepoint_locations=matrix(ncol=2,nrow=0)
   colnames(changepoint_locations)=c('dim','t')

   if(scenario_number==1){
     # Here gamma*p/2 dimensions undergo changepoints t=.25
     # and another gamma*p/2 dimensions undergo changepoints t=.75
     first_cp_location=round(n*.25)+1
     second_cp_location=round(n*.75)+1
     dims_undergoing_1_cp=ceiling(p*gamma/2-.00000001)

     # First CP
     for(i in seq(1,dims_undergoing_1_cp)){
       jump_matrix[i,first_cp_location:n]=(single_local_sd(x_matrix,min(s),max(s),i,.25)*
                                             signal_to_noise*
                                             ((stats::runif(1)>.5)*2-1))

       changepoint_locations=rbind(changepoint_locations,c(i,.25))
     }

     # Second CP
     for(i in seq(dims_undergoing_1_cp+1,2*dims_undergoing_1_cp)){
       jump_matrix[i,second_cp_location:n]=(single_local_sd(x_matrix,min(s),max(s),i,.75)*
                                              signal_to_noise*
                                              ((stats::runif(1)>.5)*2-1))

       changepoint_locations=rbind(changepoint_locations,c(i,.75))
     }

   }else if(scenario_number==2){

     # Here the first changepoints occurs at t=.2, the last at t=.8
     # and all others are evenly spaced in between this
     first_cp_location=round(n*.2)+1
     last_cp_location=round(n*.8)+1
     other_dims_undergoing_cp=ceiling(p*gamma-.000001)-2

     # Locations of other changepoints
     other_cp_locations=round(seq(first_cp_location,
                                  last_cp_location,
                                  length.out=other_dims_undergoing_cp+2)[2:(2+other_dims_undergoing_cp-1)])

     # Add first changepoint
     dim_counter=1

     jump_matrix[dim_counter,first_cp_location:n]=(single_local_sd(x_matrix,min(s),max(s),dim_counter,.2)*
                                                     signal_to_noise*
                                                     ((stats::runif(1)>.5)*2-1))
     changepoint_locations=rbind(changepoint_locations,c(dim_counter,.2))

     # Add middle changepoints
     for(i in other_cp_locations){
       dim_counter=dim_counter+1
       t_of_cp=(i-1)/n

       jump_matrix[dim_counter,i:n]=(single_local_sd(x_matrix,min(s),max(s),dim_counter,t_of_cp)*
                                       signal_to_noise*
                                       ((stats::runif(1)>.5)*2-1))
       changepoint_locations=rbind(changepoint_locations,c(dim_counter,t_of_cp))
     }

     # Add last changepoint
     dim_counter=dim_counter+1
     jump_matrix[dim_counter,last_cp_location:n]=(single_local_sd(x_matrix,min(s),max(s),dim_counter,.8)*
                                                    signal_to_noise*
                                                    ((stats::runif(1)>.5)*2-1))
     changepoint_locations=rbind(changepoint_locations,c(dim_counter,.8))

   }else if(scenario_number==10){
     # Here the first dimension undergoes a changepoint at t=.5
     # this is used in the simulated rejection probabilities section of the
     # paper
     first_cp_location=round(n*.5)+1

     dims_undergoing_1_cp=1

     jump_matrix[1,first_cp_location:n]=(single_local_sd(x_matrix,min(s),max(s),1,.5)*
                                           signal_to_noise*
                                           ((stats::runif(1)>.5)*2-1))

     changepoint_locations=rbind(changepoint_locations,c(1,.5))


   }

   return(list(jump_matrix,changepoint_locations))
 }

#' Evaluate detected jumps
#'
#' @description This functions takes a set of estimated jump locations, and a set of
#' true jump locations, and returns a number of performance metrics
#' @param estimated_changepoint_locations a dataframe with jump locations, this dataframe must have columns 'dim' and 't'
#' representing the dimension the jump occured in and the time, in  \[0,1\] that it occurred. Note it can contain other columns
#' but this will not affect the evaluation
#' @param true_changepoint_locations the locations of the jump (as they were added when the data was generated), contains column 'dim' and 't'
#' @param margin_of_error tolerance for how far away from true value, an estimated value can be for it to be considered a successful
#' idenification, units in time in \[0,1\] (e.g. could be .0001)
#' @param p dimension of times series
#' @param n length of time series
#' @return a vector with a number of performance metrics
#' @export
evaluate_changepoints<-function(estimated_changepoint_locations,
                                 true_changepoint_locations,
                                 margin_of_error,
                                 p,n){

   # Setup Experiment Results DF
   results_vector <- data.frame(m_n=numeric(),      # true number of CP
                                m_hat_n=numeric(),      # number of identified CP
                                correct_cp=numeric(),   # number of correctly identified CP
                                dim_cp=numeric(),       # number of dimensions with CP
                                dim_hat_cp = numeric(), # number of dimensions where we correctly ID'd number of CP
                                MAD = numeric(),        # mean abolute deviation of estimated vs  actual value of CPs
                                stringsAsFactors=FALSE)

   # Experiment 'correct' values
   if(length(true_changepoint_locations)==2){m_n=1
   }else{
     m_n = length(true_changepoint_locations[,1])
   }

   results_vector[1,'m_n']=m_n # actual number jumps
   results_vector[1,'m_hat_n']=dim(estimated_changepoint_locations)[1] # number of jumps estimated

   # Vector with the number of true jumps per dimension
   true_cp_per_dim = rep(0,p)

   # Case where more than one jump
   if(length(true_changepoint_locations)==2){
     true_cp_per_dim[true_changepoint_locations[,'dim'][1]]=1
   }else{
     for(j in c(true_changepoint_locations[,1])){
       true_cp_per_dim[j]=true_cp_per_dim[j]+1
     }
   }

   est_cp_per_dim=rep(0,p)
   cor_cp_per_dim=rep(0,p)
   MAD_total=0

   # Go through each dimension and determine how many jumps were correctly estimated

   # Case where only one jump
   if(length(true_changepoint_locations)==1){
     correct_cp_count = 0

     cp_dim=true_changepoint_locations[,'dim']
     cp_t=true_changepoint_locations[,'t']
     est_cp_per_dim[cp_dim]=1

     list_of_cp_to_check=estimated_changepoint_locations[estimated_changepoint_locations[,'dim']%in%cp_dim,'t_index']/n

     for(h in list_of_cp_to_check){
       if(n*abs(h-cp_t)<=(round(n*margin_of_error))){
         correct_cp_count=correct_cp_count+1
         MAD_total=MAD_total+abs(h-cp_t)
       }
     }

     cor_cp_per_dim[cp_dim]=correct_cp_count
   }else{

     for(j in seq(1,p)){
       correct_cp_count = 0 # correct number of jumps in a dimension
       list_of_cp_to_check=estimated_changepoint_locations[estimated_changepoint_locations[,'dim']==j,'t_index']
       est_cp_per_dim[j]=length(list_of_cp_to_check)

       # if in a given dimension there are actually jump then check if we correctly estimated them
       if(length(true_changepoint_locations[true_changepoint_locations[,'dim']==j,'t'])>0){
         for(h in list_of_cp_to_check){
           # for an estimated jump, find the closest true jump to it
           closest_true_cp = which.min(abs(h/n-true_changepoint_locations[true_changepoint_locations[,'dim']==j,'t']))
           closest_true_cp = true_changepoint_locations[true_changepoint_locations[,'dim']==j,'t'][closest_true_cp]

           # if for a given jump it is within the margin of error of a true jump
           # we count it as a correct jump
           # note I use round here because if I did not use round it was giving me strange errors,
           # I suspect related to how margin_of_error is stored in memory
           if(round(n*abs(h/n-closest_true_cp))<=round(n*margin_of_error)){
             correct_cp_count=correct_cp_count+1
             # add estimation error to MAD total
             MAD_total=MAD_total+abs(h/n-closest_true_cp)
           }

         }
         cor_cp_per_dim[j]=correct_cp_count
       }
     }
   }

   results_vector[1,'correct_cp']=sum(cor_cp_per_dim)
   results_vector[1,'dim_cp']=sum(true_cp_per_dim>0)
   results_vector[1,'dim_hat_cp']=sum((cor_cp_per_dim==true_cp_per_dim)&(true_cp_per_dim>0))
   if(sum(cor_cp_per_dim)==0){
     results_vector[1,'MAD']=0
   }else{
     results_vector[1,'MAD']=MAD_total/sum(cor_cp_per_dim)
   }

   return(results_vector)

 }

#' Power summary statistics
#'
#' @description This function takes the result from K power experiments, and outputs summary statistics to a CSV file
#' @param results_df dataframe of resutls of power simulations
#' @param K number of experiments
#' @param scenario scenario from simulation
#' @param test_type type of test (power, type 1)
#' @param gamma gamma used
#' @param s scales used
#' @param p dimension of time series
#' @param jump_size size of jump for time series
#' @param sim_result_dir where to output results
#' @export
 power_summary_statistics<-function(results_df,
                                    K,
                                    scenario,
                                    test_type,
                                    gamma,
                                    s,
                                    p,
                                    jump_size,
                                    sim_result_dir){

   # This function takes the result from K power experiments, and outputs summary
   # statistics to a CSV file

   # INPUTS
   # results_df - a dataframe of the results of the power simulations (with inputs
   #              from function evaluate_changepoints)
   # K - num experiments
   # scenario, test_type, gamma,s,p,jump_size - characteristics of an experiment
   # sim_result_dir - where the CSV with the summary statistics is written to

   # OUTPUT
   # CSV file written to sim_result_dir
   # After all experiments for a given set of parameters has been run
   pct_perfect = sum((results_df[,'m_n']==results_df[,'m_hat_n'])&(results_df[,'m_n']==results_df[,'correct_cp']))/K
   pct_perfect_dim = sum(results_df[,'dim_cp']==results_df[,'dim_hat_cp'])/K
   avg_correct_cp = mean(results_df[results_df[,'m_n']>0,'correct_cp']/results_df[results_df[,'m_n']>0,'m_n'])
   avg_correct_dim = mean(results_df[results_df[,'m_n']>0,'dim_hat_cp']/results_df[results_df[,'m_n']>0,'dim_cp'])
   avg_MAD = mean(results_df[,'MAD'])
   avg_m_hat_n = mean(results_df[,'m_hat_n'])

   # placeholder values of -1 for avg correct cp/avg correct dim (if never a CP)
   if(sum(results_df[,'m_n']>0)==0){
     avg_correct_cp=-1
     avg_correct_dim=-1
   }

   results_file_name = paste('P',
                             'sce',scenario,
                             'param',gamma,
                             test_type,
                             'p',sprintf("%04d", p),
                             min(s),max(s),
                             'jmp',jump_size,
                             '%per',round(pct_perfect,2),
                             '%perdim',round(pct_perfect_dim,2),
                             'avgcorcp',round(avg_correct_cp,2),
                             'avgcordim',round(avg_correct_dim,2),
                             'avgMAD1k',round(avg_MAD*1000,2),
                             'm_hat_n',round(avg_m_hat_n,0))

   setwd(sim_result_dir)
   utils::write.csv(results_df,file=paste(results_file_name,'.csv'),row.names=FALSE)
   print(results_file_name)

 }

#' Evaluate detected jumps
#'
#' @description This function takes a list of jumps identified by the four methods tested and compares them to the true jump
#' locations, and then updates results_df appropriately based on the scenario note signal_to_noise an input to deal with case where jump_size=0
#' @param true_changepoint_locations a df with the dim and time of the true jumps
#' @param margin_of_error a value within \[0,1\] about how far away a detected jump can be to be labelled as correct
#' @param results_df a dataframe with parameters from the simulation, and what results we are recording
#' @param changepoints_by_method_list a list which has the jumps detected by different methods the entries are as follows
#' 1st entry INSPECT , 2nd method DBLCUSUM, 3rd entry LOCLIN, 4th entry AJDN. First 3 entries are vectors, 4th entry a dataframe
#' with dimension and location of jumps.
#' @param scenario power scenario
#' @param n length of time series
#' @param K number of experiments
#' @param signal_to_noise SNR of jump
#' @param gamma proportion of dimensions undergoing jump
#' @return dataframe of results
#' @export
 evaluate_changepoints_compare <- function(true_changepoint_locations,
                                           margin_of_error,
                                           results_df,
                                           changepoints_by_method_list,
                                           scenario,n,K,
                                           signal_to_noise,gamma){
   # This function takes a list of jumps identified by the four methods tested and
   # compares them to the true jump locations, and then updates results_df
   # appropriately based on the scenario
   # note signal_to_noise an input to deal with case where jump_size=0

   # INPUT
   # true_changepoint_locations  - a df with the dim and time of the true jumps
   # margin_of_error             - a value within [0,1] about how far away a detected jump
   #                               can be to be labelled as correct
   # results_df                  - a dataframe with parameters from the simulation, and what results
   #                               we are recording
   # changepoints_by_method_list - a list which has the jumps detected by different methods
   #                               the entries are as follows [[1]] sparse, [[2]] dblcusum
   #                                                          [[3]] chenwangwu, [[4]] mjpdhd
   #                               [[1]],[[2]],[[3]] are vectors, whereas [[4]] is a dataframe with
   #                               both dimension and time where the cp was located
   # scenario, n, k, snr, gamma  - parameters from the simulation being run

   # OUTPUT
   # results_df
   #
   # NOTE
   # here for the calculation of MJPD MAD we are assuming that in a dimension jumps
   # are at least 1 margin of error apart (this is fine for our simulations)

   for(i in seq(1,4)){

     cp_located=changepoints_by_method_list[[i]]

     # for AJDN jumps need to take only the times they were located at
     if(i==4){
       cp_located=unique(cp_located[,'t_index'])
     }

     # Set MAD to 999, if no jumps correctly identified
     results_df[i,'MAD']=999

     if(scenario==1){
       # Scenario 1 #
       # CP in Gamma/2 dimensions at t=.25, t=.75

       true_cp_indexes=unique(true_changepoint_locations[,'t'])*n
       distances=rep(0,length(true_cp_indexes))
       final_cp_list=c()

       if(length(cp_located)>0){
         for(j in seq(1,length(true_cp_indexes))){
           distances[j]=min(abs(cp_located-true_cp_indexes[j]))
         }

         # Finds CPs correctly detected, and then uses this to calc MAD
         if(min(distances)<=round(margin_of_error*n)){
           # only count distances that are within margin of error
           MAD_distances=distances[distances<=round(margin_of_error*n)]
           results_df[i,'MAD']=mean(MAD_distances)
         }

         # If no CP, this makes every CP a large distance
         if(signal_to_noise==0){
           distances=distances+10000
         }

         # Determine if all jumps were identified correctly
         if(max(distances)<=round(margin_of_error*n)){
           results_df[i,'correct']=1

           # Determine if all CP identified within margin of error
           # of some true jump
           min_distances=rep(10000,length(cp_located))
           for(j in 1:length(min_distances)){
             min_distances[j]=min(abs(cp_located[j]-true_cp_indexes))
           }

           if(max(min_distances)<=round(margin_of_error*n)){
             results_df[i,'perfect']=1
           }
         }

         # Check if we have identified each of the jumps
         for(j in seq(1,length(true_cp_indexes))){
           if(distances[j]<=round(margin_of_error*n)){
             final_cp_list=c(final_cp_list,true_cp_indexes[j])
           }
         }

         # jumps that were not close to true CP
         incorrect_cp = c()

         for(j in seq(1,length(cp_located))){
           if(min(abs(cp_located[j]-true_cp_indexes))>round(margin_of_error*n)){
             incorrect_cp=c(incorrect_cp,cp_located[j])
           }
         }

         final_cp_list=c(final_cp_list,incorrect_cp)

         results_df[i,'cp_found']=length(final_cp_list)

       }

     }else if(scenario==2){
       # Scenario 2
       # First cp at t=.2, last at t=.8
       # Others spaced in between

       true_cp_indexes=unique(true_changepoint_locations[,'t'])*n
       distances=rep(0,length(true_cp_indexes))
       final_cp_list=c()

       if(length(cp_located)>0){
         # Here because the jumps are potentially close together
         # we need to take care to not double count jumps being
         # identified
         #
         # To do this, we assume that the jumps found closest
         # to the true jumps

         for(j in seq(1,length(true_cp_indexes))){
           distances[j]=min(abs(cp_located-true_cp_indexes[j]))
         }

         # If no CP, this makes every CP a large distance
         if(signal_to_noise==0){
           distances=distances+10000
         }

         # Finds CPs correctly detected, and then uses this to calc MAD
         if(min(distances)<=round(margin_of_error*n)){
           # only count distances that are within margin of error
           MAD_distances=distances[distances<=round(margin_of_error*n)]
           results_df[i,'MAD']=mean(MAD_distances)
         }

         # Determine if all jumps were identified correctly
         if(max(distances)<=round(margin_of_error*n)){
           results_df[i,'correct']=1

           # Determine if all CP identified within margin of error
           # of some true jump
           min_distances=rep(10000,length(cp_located))
           for(j in 1:length(min_distances)){
             min_distances[j]=min(abs(cp_located[j]-true_cp_indexes))
           }

           if(max(min_distances)<=round(margin_of_error*n)){
             results_df[i,'perfect']=1
           }
         }

         # Check if we have identified each of the jumps
         for(j in seq(1,length(true_cp_indexes))){
           if(distances[j]<=round(margin_of_error*n)){
             final_cp_list=c(final_cp_list,true_cp_indexes[j])
           }
         }

         # jumps that were not close to true CP
         incorrect_cp = c()

         for(j in seq(1,length(cp_located))){
           if(min(abs(cp_located[j]-true_cp_indexes))>round(margin_of_error*n)){
             incorrect_cp=c(incorrect_cp,cp_located[j])
           }
         }

         final_cp_list=c(final_cp_list,incorrect_cp)

         results_df[i,'cp_found']=length(final_cp_list)

       }

     }

     # MAD Adjustment for MJPD
     # here we are calculating MAD for MJPD in a different way than for the other methods
     # mostly what we are doing is averaging the absolute deviation of correctly identified
     # jumps across all dimensions (so if there are two CP identified within the
     # margin of error than we average the absolute deviation of both)
     if(i==4){
       absolute_deviation_sum=0
       number_of_cps=0
       mjpd_cp_df=changepoints_by_method_list[[i]]
       mjpd_cp_df$dim=as.numeric(mjpd_cp_df$dim)
       mjpd_cp_df$t_index=as.numeric(mjpd_cp_df$t_index)

       # a dataframe which has each true CP and finds if any CP detected it
       true_cp_compare_detected=data.frame(true_changepoint_locations)
       true_cp_compare_detected[,'detected_cp']=-1

       for(l in seq(1,dim(true_cp_compare_detected)[1])){
         true_cp_index=round(true_cp_compare_detected$t[l]*n)

         # which dimension this specific CP was found in
         cp_detected_in_dim=mjpd_cp_df[mjpd_cp_df$dim==true_cp_compare_detected$dim[l],]

         # if at least one CP detected in the dimension we care about
         if(nrow(cp_detected_in_dim)>0){
           # decide if a detected jump falls within the margin of error
           # of this jump
           min_distance=min(abs(cp_detected_in_dim$t_index-true_cp_index))

           if(min_distance<=round(margin_of_error*n)){
             number_of_cps=number_of_cps+1
             absolute_deviation_sum=absolute_deviation_sum+min_distance
           }
         }

       }

       if(number_of_cps>0){
         results_df[i,'MAD']=absolute_deviation_sum/number_of_cps
       }else{
         results_df[i,'MAD']=999
       }
     }

   }
   return(results_df)

 }

#' Calculate BIC
#'
#' @description This function calculates BIC on an individual dimension of a time series based on the changepoints that were detected
#' @param time_series a vector which is a time series
#' @param cps_detected dataframe of detected jumps (with column t_index)
#' @param s_prime \eqn{s'}
#' @return numerical calcluation for BIC
#' @export
 calc_bic<-function(time_series,cps_detected,s_prime){

   n=length(time_series)
   cps_detected=sort(cps_detected$t_index)

   # case of no changepoints just calculate error variance
   if(length(cps_detected)==0){

     # Fit a nonparametric trend
     x=seq(1,length(time_series))
     y=as.numeric(time_series)
     auto_bandwidth=KernSmooth::dpill(x,y)
     np_mean_est=KernSmooth::locpoly(x,y,bandwidth = auto_bandwidth,gridsize=length(x))$y

     segment_residuals_raw=y-np_mean_est

     bic=n*log(sum(segment_residuals_raw**2)/n)+0

   }else{
     start_indexes=c(1,cps_detected+1)
     end_indexes=c(cps_detected,n)
     error_segments=rep(0,length(start_indexes))
     residuals=c(0)

     bandwidth_candidates=rep(0,length(start_indexes))
     for(k in 1:(length(start_indexes))){
       start_index=start_indexes[k]
       end_index=end_indexes[k]
       segment_length=end_index-start_index+1
       x=seq(1,segment_length)
       y=as.numeric(time_series[start_index:end_index])

       # Can throw errors if segment size very small, if this is the case
       # use a default bandwidth of 1
       bandwidth_candidates[k]=try(KernSmooth::dpill(x,y),silent=TRUE)

     }

     # Choose middle bandwidth (if odd number take slightly bigger one)
     bandwidth_selected=stats::median(sort(bandwidth_candidates))

     for(k in 1:(length(start_indexes))){
       start_index=start_indexes[k]
       end_index=end_indexes[k]
       segment_length=end_index-start_index+1

       x=seq(1,segment_length)
       y=as.numeric(time_series[start_index:end_index])

       np_mean_est=KernSmooth::locpoly(x,y,bandwidth = bandwidth_selected,gridsize=length(x))$y

       residuals=c(residuals,y-np_mean_est)

     }

     # Calculate BIC by Zhou Frequency Detection and CP of Complex Oscillation
     bic=n*log(sum(residuals**2)/n)+length(cps_detected)*log(n)

   }

   return(bic)

 }
