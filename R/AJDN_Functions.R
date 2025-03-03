#' Create W filter
#'
#' @description Creates a vector of the optimal jump pass filter for a \eqn{n} (length of time series), \eqn{s} (scale), and \eqn{t} (time)
#'
#' @param n length of full time series
#' @param scale_count rounded value of \eqn{ns} where \eqn{s} is the scale
#' @param t_index the index associated with the time at the center of the filter (\eqn{t_jn})
#' @return a vector of length \eqn{n} with weights that the filter takes
#' @export
#'
#' @examples
#' w = ajdn::create_W_filter(1000,100,300)
#' plot(w)
create_W_filter <- function(n,scale_count,t_index){
  # OPTIMAL WAVELET FILTER
  W_row = rep(0,n)

  # Values for filter on LHS of t_round
  LHS_indicies = seq(t_index-scale_count,t_index-1) # indicies we are using
  LHS_wave_inputs = (LHS_indicies/n-t_index/n)/(scale_count/n)# indicies converted to scale of [-1,0]

  # Filter equation from paper
  W_row[(t_index-scale_count):(t_index-1)] = ((-1294.2222)*(LHS_wave_inputs**6) +
                                                4246.6667*(abs(LHS_wave_inputs)**5) +
                                                (-5320)*(LHS_wave_inputs**4) +
                                                3188.8889*(abs(LHS_wave_inputs)**3) +
                                                (-933.3333)*(LHS_wave_inputs**2) +
                                                (112)*(abs(LHS_wave_inputs)))*sign(LHS_wave_inputs)

  # right at the boundary the filter should equal exactly 0
  W_row[(t_index-scale_count)] = 0
  W_row[(t_index+1):(t_index+scale_count)] = rev(W_row[(t_index-scale_count):(t_index-1)]*(-1))

  return(W_row)
}

#' Generate bootstrap matrix
#'
#' @description Generates the matrix used in the high dimensional block multiplier bootstrap
#'
#' @param x_matrix \eqn{p} x \eqn{p} (\eqn{p} = # dimensions, \eqn{n} = # obsrevations) matrix of time series
#' @param s_prime block size parameter
#' @return \eqn{p} x \eqn{n} matrix to be subsequently used in the block multiplier bootstrap
#' @export
#' @examples
#' x_matrix = matrix(stats::rnorm(1000),nrow=10)
#' upsilon_matrix = ajdn::generate_upsilon_matrix(x_matrix,1/100)
generate_upsilon_matrix <- function(x_matrix,s_prime){
  n = dim(x_matrix)[2]
  p = dim(x_matrix)[1]

  upsilon_matrix = matrix(0,nrow=p,ncol=n)

  # Since s_prime*n may be a fraction, and we need a whole
  # number to construct the bootstrap so we round s_prime * n to do so
  s_prime_round = round(s_prime*n)

  for(i in seq(1,n)){
    if(i<=(s_prime_round+1)){
      upsilon_matrix[,i] = (rowSums(x_matrix[,1:s_prime_round,drop=FALSE])-rowSums(x_matrix[,(s_prime_round+1):(2*s_prime_round),drop=FALSE]))/sqrt(2*s_prime_round)
    } else if((i>(s_prime_round+1))&(i<=(n-s_prime_round))){
      upsilon_matrix[,i] = (rowSums(x_matrix[,(i-s_prime_round):(i-1),drop=FALSE])-rowSums(x_matrix[,i:(i+s_prime_round-1),drop=FALSE]))/sqrt(2*s_prime_round)

    } else if(i>(n-s_prime_round)){
      upsilon_matrix[,i] = (rowSums(x_matrix[,(n-2*s_prime_round+1):(n-s_prime_round),drop=FALSE])-rowSums(x_matrix[,(n-s_prime_round+1):n,drop=FALSE]))/sqrt(2*s_prime_round)
    }
  }

  return(upsilon_matrix)
}

#' Second stage refinement of detected jumps
#'
#' @description Conducts the second stage procedure to refine estimation of jumps
#'
#' @param first_stage_changepoints a list where each entry of the list represents a jump with dim (dimension),
#' t_index (\eqn{t_j n}), scale_count (\eqn{ns} where s scale used to find jump)
#' e.g. fifth entry of list, attribute t_index is the time index the 5th jump identified was found at
#' @param x_matrix \eqn{p} x \eqn{n} (\eqn{p} = # dimensions, \eqn{n} = # obsrevations) matrix of time series
#' @param scale_count_list a list of length \eqn{p}, where the \eqn{p}th entry of the list is a vector of \eqn{ns_{r,1},...,ns_{r,\delta_n}} where
#' \eqn{\delta_n} is the number of scales being used in dimension \eqn{r}
#' @param ss_refine_alpha_tilde constant from paper
#' @return a list of same length as first_stage_changepoints with the refined estimate of the time index where the jump occurs
#' @export
second_stage_refine <-function(first_stage_changepoints,x_matrix,scale_count_list,ss_refine_alpha_tilde){

  n = dim(x_matrix)[2]
  second_stage_changepoints = first_stage_changepoints
  count=1

  # Iterate through each jump found in the first stage and conduct the second
  # stage procedure to refine the estimation of the jump
  for(cp in first_stage_changepoints){
    cp_t = cp['t_index']
    cp_dim = cp['dim']
    z_n = min(scale_count_list[[cp_dim]])/2

    # Calculate l_r,u_r,l_r_tilde,u_r_tilde for a given n
    reg_window = round((2+ss_refine_alpha_tilde)*z_n)
    tilde_window = z_n

    if(cp_t+reg_window>n){
      reg_window=n-cp_t
      z_n=reg_window
    }else if(cp_t-reg_window<1){
      reg_window=cp_t-1
      z_n=reg_window
    }

    l_r = (cp_t-reg_window)
    u_r = (cp_t+reg_window)

    l_r_tilde = (cp_t-tilde_window)
    u_r_tilde = (cp_t+tilde_window)

    cusum_stats = rep(0,n)

    S_lr_ur = sum(x_matrix[cp_dim,l_r:u_r])
    lambda_lr_ur = u_r-l_r+1

    s_lj_t=cumsum(as.numeric(x_matrix[cp_dim,l_r:u_r]))
    lambda_frac=(seq(l_r,u_r)-l_r+1)/lambda_lr_ur
    cusum_stats[seq(l_r,u_r)]=s_lj_t-lambda_frac*S_lr_ur

    max_evidence = max(abs(cusum_stats[l_r_tilde:u_r_tilde]))
    max_evidence_index = which(abs(cusum_stats)==max_evidence)[1]

    second_stage_changepoints[[count]]['t_index'] = max_evidence_index

    count=count+1
  }

  return(second_stage_changepoints)
}


#' Calculate estimate of local standard deviation of time series
#'
#' @description This function calcualtes an estimate of local standard deviation at each time point
#'
#' @param x_matrix p x n (p = # dimensions, n = # obsrevations) matrix of time series
#' @param t_index_list t_index_list a list of time indexes (\eqn{n t_j}) we need to calculate local standard deviations for for each dimension
#' @param scale_count_list a list of length \eqn{p}, where the \eqn{p}th entry of the list is a vector of \eqn{ns_{r,1},...,ns_{r,\delta_n}} where
#' \eqn{\delta_n} is the number of scales being used in dimension \eqn{r}
#' @param num_cores OPTIONAL argument, number of cores to use for parallelization
#' @param loud OPTIONAL, if TRUE print message when dimensions are normalized
#' @return a dataframe, with the local sd, time index, and dimension
#' @importFrom foreach %dopar%
#' @export
time_series_normalization <-function(x_matrix,t_index_list,scale_count_list,num_cores=NULL,loud=NULL){
  p = dim(x_matrix)[1]
  n = dim(x_matrix)[2]
  dimension = NULL # solves a r package undefined global variable problem

  if(is.null(num_cores)){
    normalization_df=data.frame(sd=numeric(),t_index=numeric(),dimension=numeric())

    for(i in seq(1,p)){
      min_scale_count=min(scale_count_list[[i]])
      max_scale_count=max(scale_count_list[[i]])
      window_width = max_scale_count-min_scale_count + 1 # +1 given we include the endpoints of the interval

      for(t_index in t_index_list[[i]]){

        left_mean = mean(x_matrix[i,(t_index-max_scale_count):(t_index-min_scale_count)])
        right_mean = mean(x_matrix[i,(t_index+min_scale_count):(t_index+max_scale_count)])
        sample_var = (sum((x_matrix[i,(t_index-max_scale_count):(t_index-min_scale_count)]-left_mean)**2)+
                        sum((x_matrix[i,(t_index+min_scale_count):(t_index+max_scale_count)]-right_mean)**2))/(2*window_width)

        normalization_df=rbind(normalization_df,
                               data.frame(sd=sqrt(sample_var),t_index=t_index,dimension=i))
      }
      if(!is.null(loud)){print(paste('Normalized dimension',i))}
    }
  }else{
    cl <- parallel::makeCluster(num_cores,outfile="")
    doParallel::registerDoParallel(cl)

    normalization_df=foreach::foreach(dimension=seq(1,p),.combine='rbind') %dopar% {

                               normalization_df_single_dim=data.frame(sd=numeric(),t_index=numeric(),dimension=numeric())
                               min_scale_count=min(scale_count_list[[dimension]])
                               max_scale_count=max(scale_count_list[[dimension]])
                               window_width = max_scale_count-min_scale_count + 1 # +1 given we include the endpoints of the interval

                               for(t_index in t_index_list[[dimension]]){

                                 left_mean = mean(x_matrix[dimension,(t_index-max_scale_count):(t_index-min_scale_count)])
                                 right_mean = mean(x_matrix[dimension,(t_index+min_scale_count):(t_index+max_scale_count)])
                                 sample_var = (sum((x_matrix[dimension,(t_index-max_scale_count):(t_index-min_scale_count)]-left_mean)**2)+
                                                 sum((x_matrix[dimension,(t_index+min_scale_count):(t_index+max_scale_count)]-right_mean)**2))/(2*window_width)

                                 normalization_df_single_dim=rbind(normalization_df_single_dim,
                                                                   data.frame(sd=sqrt(sample_var),t_index=t_index,dimension=dimension))
                               }
                               if(!is.null(loud)){print(paste('Normalized dimension',dimension))}
                               normalization_df_single_dim
                             }
    parallel::stopCluster(cl)
  }

  return(normalization_df)

}

#' Calculate estimates of local standard deviation of time series when scales same for all dimensions
#'
#' @description This function calculates an estimate of the standard deviation of a time series at each time point
#' when the scales \eqn{s} are the same for every dimension. This is generally significantly faster than time_series_normalization.
#'
#' @param x_matrix \eqn{p} x \eqn{n} (\eqn{p} = # dimensions, \eqn{n} = # obsrevations) matrix of time series
#' @param t_index_list a list of time indexes (\eqn{n t_j}) we need to calculate local standard deviations for for each dimension
#' @param scale_count_list a list of length \eqn{p}, where the \eqn{p}th entry of the list is a vector of \eqn{ns_{r,1},...,ns_{r,\delta_n}} where
#' \eqn{\delta_n} is the number of scales being used in dimension \eqn{r}
#' @param loud OPTIONAL, if TRUE print message when dimensions are normalized
#' @return a dataframe, with the local sd, time index, and dimension
#' @importFrom foreach %dopar%
#' @export
time_series_normalization_fast <-function(x_matrix,t_index_list,scale_count_list,loud=NULL){
  p = dim(x_matrix)[1]
  n = dim(x_matrix)[2]

  max_scale_count=max(scale_count_list[[1]])
  min_scale_count = min(scale_count_list[[1]])
  t_indexes=t_index_list[[1]] # since same s for all dimensions, these are same for all dimensions

  # If we only have one scale then min_scale_count==max_scale_count which means window size is 1 which is
  # not ideal, in this case we will min_scale_count to round(max_scale_count/2) to give us a reasonable window
  # for what the variances/means should be
  if(min_scale_count==max_scale_count){
    min_scale_count=round(max_scale_count/2)
  }

  window_size=max_scale_count-min_scale_count+1

  # Matrix we will input all of the variance entries into
  var_matrix=matrix(data=0,nrow=p,ncol=n)
  x_matrix_sq=x_matrix**2

  # In the below code we calculate all of the variances incrementally, since
  # the formula for variance at t, can be calculated by slightly adjusting the
  # formula for the variance at t-1 (i.e. dont have to sum over a bunch of values
  # in the x_matrix, only need to access a few elements)

  # Calculate means and sum of squares of time series, can then use these values to quickly
  # calculate variances
  t_index=t_indexes[1]
  left_mean=rowMeans(x_matrix[,(t_index-max_scale_count):(t_index-min_scale_count),drop=FALSE])
  right_mean=rowMeans(x_matrix[,(t_index+min_scale_count):(t_index+max_scale_count),drop=FALSE])
  left_sum_sq=rowSums(x_matrix_sq[,(t_index-max_scale_count):(t_index-min_scale_count),drop=FALSE])
  right_sum_sq=rowSums(x_matrix_sq[,(t_index+min_scale_count):(t_index+max_scale_count),drop=FALSE])

  var_matrix[,t_index]=((left_sum_sq-window_size*left_mean**2)+
                          (right_sum_sq-window_size*right_mean**2))/(2*window_size)

  for(t_index in t_indexes[2:length(t_indexes)]){
    # Update values of left_mean, right_mean, left_sum_sq, right_sum_sq
    # and then update var_matrix with them
    left_mean=left_mean+(x_matrix[,t_index-min_scale_count]-x_matrix[,t_index-1-max_scale_count])/window_size
    right_mean=right_mean+(x_matrix[,t_index+max_scale_count]-x_matrix[,t_index-1+min_scale_count])/window_size

    left_sum_sq=left_sum_sq+x_matrix_sq[,t_index-min_scale_count]-x_matrix_sq[,t_index-1-max_scale_count]
    right_sum_sq=right_sum_sq+x_matrix_sq[,t_index+max_scale_count]-x_matrix_sq[,t_index-1+min_scale_count]

    var_matrix[,t_index]=((left_sum_sq-window_size*left_mean**2)+
                            (right_sum_sq-window_size*right_mean**2))/(2*window_size)
  }

  # Create a dataframe with which to return the standard deviations that were found. Note we return this as a dataframe
  # rather than a matrix because it can then be merged onto the dataframe of raw test statistics
  normalization_df=data.frame(sd=numeric(),t_index=numeric(),dimension=numeric())

  min_t_index=min(t_indexes)
  max_t_index=max(t_indexes)

  # Note that the variance matrix may sometimes be slightly negative due to a numerical
  # precision error, therefore I will take its absolute value
  var_matrix=abs(var_matrix)

  for(i in 1:p){
    normalization_df=rbind(normalization_df,
                           data.frame(sd=sqrt(var_matrix[i,min_t_index:max_t_index]),
                                      t_index=t_indexes,dimension=i))
  }

  return(normalization_df)
}

#' Create filter matrix
#'
#' @description This function creates a filter matrix, as well as a filter matrix dictionary
#' which is then used to calculate both test statistics and bootstrap statistics for AJDN. The entires of this matrix
#' are equal to the \eqn{W((j/n-t)/s)} components of \eqn{H(t,s,r)} for various values of \eqn{j,t,s}.
#'
#' @param n \eqn{n} = # obsrevations of time series
#' @param unique_scale_counts vector of \eqn{s_jn} for all unique \eqn{s_j} that are being used across all dimensions
#' @return List of two items
#' 1st entry - sparse matrix that represents
#' 2nd entry - filter matrix dictionary which specifies what scale_count (i.e. \eqn{sn}) and \eqn{t_j n} each
#' row of the filter matrix corresponds to.
#' @export
#'
#' @examples
#' filter_matrix = ajdn::create_filter_matrix(1000,c(50,100))
#' plot(filter_matrix[[1]][100,])
#' print(filter_matrix[[2]][100,])
create_filter_matrix <- function(n,unique_scale_counts){
  # Create a filter matrix as a sparse matrix. It is necessary to use a sparsse matrix because the filter
  # matrix can be very large and contains lots of 0's therefore if it was not sparse it would take up too much memory.
  # Here we create a 0 row matrix and gradually append its rows.
  filter_matrix=matrix(ncol=n,nrow=0)
  filter_matrix = Matrix::Matrix(filter_matrix,sparse=TRUE)
  filter_matrix_dictionary=data.frame(scale_count=numeric(),t_index=numeric())

  # Loop through scales and create all necessary filters for that scale
  for(scale_count_i in unique_scale_counts){
    # Get a list of the (n x t_j) values we need a filter for for this scale
    t_indexes_to_test=seq(scale_count_i+1,n-scale_count_i,by=1)

    # Get the filter for the first (t_j x n)
    first_t_index_filter=create_W_filter(n,scale_count_i,t_indexes_to_test[1])/sqrt(2*scale_count_i)

    # get the entries of the filter corresponding to [t_j-s_i,t_j+s_i]
    filter_entries=first_t_index_filter[1:(t_indexes_to_test[1]+scale_count_i)]

    # Create a sparse filter matrix for just one scale (will eventually it to larger filter matrix eventually)
    single_scale_filter_matrix=Matrix::Matrix(data=0,
                                      nrow=length(t_indexes_to_test),
                                      ncol=n,
                                      sparse=TRUE)

    # This code fills single_scale_filter_matrix with the filter, shifted for each value of t_j. It does this by exactly
    # specifying the row and column indexes the non-zero values of the filter should go, and then filling these entries
    # of the sparse matrix with the corresponding filter value.
    row_indicies = rep(1:length(t_indexes_to_test),each=(2*scale_count_i+1))
    column_indices=c(mapply(seq,
                            from=(t_indexes_to_test[1]-scale_count_i):(t_indexes_to_test[length(t_indexes_to_test)]-scale_count_i),
                            length.out=round(2*scale_count_i+1)))
    single_scale_filter_matrix=Matrix::sparseMatrix(i=row_indicies,j=column_indices,x=rep(filter_entries,length(t_indexes_to_test)))

    # Create the dictionary for this filter matrix
    single_scale_filter_matrix_dictionary=data.frame(scale_count=rep(scale_count_i,length(t_indexes_to_test)),
                                                     t_index=t_indexes_to_test)

    # Bind this filter matrix for a single scale to the main filter matrix, and same with the dictionary
    filter_matrix=rbind(filter_matrix,single_scale_filter_matrix)
    filter_matrix_dictionary=rbind(filter_matrix_dictionary,single_scale_filter_matrix_dictionary)
  }

  return(list(filter_matrix,filter_matrix_dictionary))

}

#' Use AJDN to detect jumps
#'
#' @description Detects jumps using AJDN algorithm
#'
#' @param x_matrix \eqn{p} x \eqn{n} (\eqn{p} = # dimensions, \eqn{n} = # obsrevations) matrix of time series
#' @param s a vector of scales to test, or a list of scales where each item in
#' list corresponds to vector of scales for that dimension. Scales should be in (0,0.5)
#' @param s_prime block size parameter. \eqn{2ns'} is equal to number total number of observations
#' in each block of the block bootstrap
#' @param B number of bootstraps to use in calculation of critical value
#' @param alpha desired type I error
#' @param filter_matrix_list OPTIONAL argument, input should be the the list output from create_filter_matrix_function. This is used
#' to speed up instances where AJDN is being run multiple times with same \eqn{n,p,s} as creating the filter matrix can
#' be computationally intensive.
#' @param loud OPTIONAL argument, if loud=TRUE then print out progress of AJDN algorithm as it runs, and each jump detected
#' @param num_cores OPTIONAL argument, number of cores to parallelize computation across, if only using 1 core leave as NULL.
#' @param output_bootstraps_and_teststat OPTIONAL argument, if set to TRUE the function returns max test stat and all initial bootstraps
#' (used to examine behaviour of algorithm)
#' @param bootstrap_normal_seed OPTIONAL argument, a seed to set before generating the standard normal random variables for the bootstrap
#' @param cp_delete_c constant used when deleting a neighbourhood around an already detected jump
#' @param ss_refine_alpha_tilde constant related to the second stage refinement
#' @return dataframe of jumps detected
#' @export
#'
#' @examples
#' x_matrix = matrix(rnorm(2000),nrow=2)
#' x_matrix[1,501:1000] = x_matrix[1,501:1000] + 5
#' ajdn::ajdn_detect_jumps(x_matrix,c(.05),.001,1000,.05)
ajdn_detect_jumps<-function(x_matrix,
                            s,
                            s_prime,
                            B,
                            alpha,
                            filter_matrix_list=NULL,
                            loud=NULL,
                            num_cores=NULL,
                            output_bootstraps_and_teststat=NULL,
                            bootstrap_normal_seed=NULL,
                            cp_delete_c=.01,
                            ss_refine_alpha_tilde=1.5){

  # PRELIMINARIES
  p=dim(x_matrix)[1]
  n=dim(x_matrix)[2]
  i = NULL # solves a r package undefined global variable problem

  preprocess_list=preprocess_s(s,p,n)
  scale_count_list=preprocess_list[[1]]
  t_index_list=preprocess_list[[2]]
  unique_scale_counts=preprocess_list[[3]]

  if(!is.null(output_bootstraps_and_teststat)){bootstrap_teststat_list=list()}

  # Setup a filter matrix and a corresponding filter matrix dictionary. The filter matrix is a big
  # matrix which is then applied to the data to calculate the unnormalized test stats (H(t,s,r) in paper).
  # The filter matrix dictionary contains information on what time index and scale_count each row of the filter matrix
  # corresponds to. If filter matrix has already been created (this can save time if running algorithm repeadetly),
  # can ignore this step.
  if(!is.null(loud)){print('Creating Filter Matrix')}

  if(is.null(filter_matrix_list)){
    filter_matrix_list=create_filter_matrix(n,unique_scale_counts)
  }

  filter_matrix=filter_matrix_list[[1]]
  filter_matrix_dictionary=filter_matrix_list[[2]]

  # CALCULATE TEST STATISTICS
  if(!is.null(loud)){print('Calculating Test Statistics')}

  # There are 3 cases here we consider (in order of code written)
  # Case 1 - all dimensions have same scales, then calculating test statistics
  #          is a simple and efficient matrix multiplication
  # Case 2 - not all dimensions have same scales, then will need to loop through
  #          each dimension to calculate it, no paralleization
  # Case 3 - like case 2 but with parallelization

  if(length(unique(scale_count_list))==1){

    # Get only filter matrix rows we want to include (important since filter matrix could have rows corresponding
    # to t_j not in T_r...particularly for rows corresponding to smallest s)
    filter_rows_to_include=((filter_matrix_dictionary$scale_count%in%scale_count_list[[1]])&(
                             filter_matrix_dictionary$t_index%in%t_index_list[[1]]))
    all_dimension_filter_dictionary=filter_matrix_dictionary[filter_rows_to_include,]

    # a matrix of size f x p where f number of filters being used
    raw_test_statistics = abs(filter_matrix[filter_rows_to_include,]%*%t(x_matrix))

    # Create a test statistics dataframe, where each row has the test statistic and corresponding scale s, time index t,
    # and dimension r
    test_statistics = data.frame(test_stat=as.vector(raw_test_statistics),
                                 scale_count=rep(all_dimension_filter_dictionary$scale_count,p),
                                 t_index=rep(all_dimension_filter_dictionary$t_index,p),
                                 dimension=rep(1:p, each=dim(all_dimension_filter_dictionary)[1]))

  }else if(is.null(num_cores)){
    # Loop over each dimension and calculate the test statistics for each scale s, time t, and dimension r
    test_statistics = data.frame(test_stat=numeric(),
                                 scale_count=numeric(),
                                 t_index=numeric(),
                                 dimension=numeric())

    for(i in seq(1,p)){
      # Which rows of the filter matrix to include (list of true/false)
      filter_rows_to_include=((filter_matrix_dictionary$scale_count%in%scale_count_list[[i]])&(
                               filter_matrix_dictionary$t_index%in%t_index_list[[i]]))
      dimension_filter_dictionary=filter_matrix_dictionary[filter_rows_to_include,]

      # Calculate the unnormalized test statistics and create a df of them
      dimension_test_statistics=abs(filter_matrix[filter_rows_to_include,]%*%x_matrix[i,])[,1]
      dimension_test_df=cbind(data.frame(test_stat=dimension_test_statistics),
                              dimension_filter_dictionary,
                              data.frame(dimension=i))

      # Append df of test statistics for dimension to df of test statistics of all dimensions
      test_statistics=rbind(test_statistics,dimension_test_df)

      if(!is.null(loud)){print(paste('Calculating Test Statistics Dim',i))}
    }
  }else{
    # Setup cluster to parallelize across
    cl <- parallel::makeCluster(num_cores,outfile="")
    doParallel::registerDoParallel(cl)

    # For each dimension calculate a dataframe of test statistics, and bind them together into one dataframe
    test_statistics=foreach::foreach(i=seq(1,p),
                            .packages = c("Matrix"),
                            .combine='rbind') %dopar% {

                              # Which rows of the filter matrix to include (list of true/false)
                              filter_rows_to_include=((filter_matrix_dictionary$scale_count%in%scale_count_list[[i]])&(
                                                       filter_matrix_dictionary$t_index%in%t_index_list[[i]]))
                              dimension_filter_dictionary=filter_matrix_dictionary[filter_rows_to_include,]

                              # Calculate the unnormalized test statistics and create a df of them
                              dimension_test_statistics=abs(filter_matrix[filter_rows_to_include,]%*%x_matrix[i,])[,1]
                              dimension_test_df=cbind(data.frame(test_stat=dimension_test_statistics),
                                                      dimension_filter_dictionary,
                                                      data.frame(dimension=i))
                              dimension_test_df

                            }
    parallel::stopCluster(cl)
  }

  # NORMALIZE TEST STATISTICS BY LOCAL STANDARD DEVIATIONS
  if(!is.null(loud)){print('Normalizing Test Statistics')}

  # If all dimensions share the same scales we can use a fast method of normalizing the test statistics
  if((length(unique(scale_count_list))==1)){
    normalization_constants=time_series_normalization_fast(x_matrix,t_index_list,scale_count_list)
  }else{
    normalization_constants=time_series_normalization(x_matrix,t_index_list,scale_count_list,num_cores)
  }

  test_statistics = merge(test_statistics,normalization_constants,by=c('t_index','dimension'))
  test_statistics[,'test_stat_normalized']=test_statistics[,'test_stat']/test_statistics[,'sd']
  test_statistics=test_statistics[,c('dimension','t_index','scale_count','test_stat_normalized')]

  # PRELIMINARIES BEFORE DETECTING JUMPS
  # Record the maximum test statistic
  if(!is.null(output_bootstraps_and_teststat)){
    bootstrap_teststat_list[['test_stat']]=max(test_statistics$test_stat_normalized)
  }

  # Create a list that will keep track of jumps
  changepoints_identified = list()
  changepoint_number = 0

  # CALCULATE BOOTSTRAP STATISTICS
  if(!is.null(loud)){print('Conducting Bootstrap')}

  # Generate upsilon matrix from x's, this matrix is used in the block multiplier bootstrap
  upsilon_matrix = generate_upsilon_matrix(x_matrix,s_prime)

  # Create a matrix of dimension p x B where each entry is the maximum test statistic
  # of dimension i for bootstrap iteration j, we will fill this matrix and then take the
  # critical values of the bootstrap from it.
  bootstrap_evidence = matrix(data=NA,nrow=p,ncol=B)

  # Create random normals for use in the bootstrap, each row is for one of the K_0 bootstrap iteration, which has n random normals associated with it
  if(!is.null(bootstrap_normal_seed)){
    set.seed(bootstrap_normal_seed)
  }
  bootstrap_normals = matrix(data=stats::rnorm(n*B),nrow=B,ncol=n)

  # Here similarly to calculating test statistics we want to speed things up if there is only one set of scales.
  # If there is only one set of scales do not need to choose rows from filter matrix each time you look at a different dimension.
  if(length(unique(scale_count_list))==1){
    all_same_scales=1

    filter_matrix_specific_dim_indexes=((filter_matrix_dictionary$scale_count%in%scale_count_list[[1]])
                                        &(filter_matrix_dictionary$t_index%in%t_index_list[[1]]))

    filter_matrix_specific_dim=filter_matrix[filter_matrix_specific_dim_indexes,]
    filter_matrix_specific_dim_dictionary=filter_matrix_dictionary[filter_matrix_specific_dim_indexes,]
    filter_matrix_specific_dim_dictionary[,'row_order']=1:nrow(filter_matrix_specific_dim_dictionary)
  }else{
    all_same_scales=0
  }

  # Here we calculate the bootstrap statistics by either looping through each dimension or parallelizing across dimensions
  if(is.null(num_cores)){
    for(i in seq(1,p)){
      if(all_same_scales!=1){
        # get the filter matrix associated with this dimension
        filter_matrix_specific_dim_indexes=((filter_matrix_dictionary$scale_count%in%scale_count_list[[i]])
                                            &(filter_matrix_dictionary$t_index%in%t_index_list[[i]]))

        filter_matrix_specific_dim=filter_matrix[filter_matrix_specific_dim_indexes,]
        filter_matrix_specific_dim_dictionary=filter_matrix_dictionary[filter_matrix_specific_dim_indexes,]
        filter_matrix_specific_dim_dictionary[,'row_order']=1:nrow(filter_matrix_specific_dim_dictionary)
      }else{
        filter_matrix_specific_dim_dictionary = filter_matrix_dictionary
        filter_matrix_specific_dim_dictionary[,'row_order']=1:nrow(filter_matrix_specific_dim_dictionary)
      }

      # Calculate the upsilon_j * z_j for each time j and each bootstrap_iteration
      upsilon_times_stdnormal=sweep(bootstrap_normals,MARGIN=2,upsilon_matrix[i,],'*')

      # Apply the filter matrix to upsilon * z matrix for each time and scale
      bootstrap_stats_specific_dim = abs(filter_matrix_specific_dim%*%t(upsilon_times_stdnormal))

      # normalization constants for each filter
      normalization_constants_specific_dim=normalization_constants[normalization_constants$dimension==i,
                                                                   c('sd','t_index')]

      # match the normalization constant with the correct time in the same order of the filter matrix
      # for this specific dimension
      filter_matrix_specific_dim_dictionary = merge(filter_matrix_specific_dim_dictionary,normalization_constants_specific_dim,by='t_index')
      filter_matrix_specific_dim_dictionary = filter_matrix_specific_dim_dictionary[order(filter_matrix_specific_dim_dictionary$row_order),]

      # Apply normalization to all of the bootstrap statistics
      bootstrap_stats_specific_dim=sweep(bootstrap_stats_specific_dim,
                                         MARGIN=1,
                                         filter_matrix_specific_dim_dictionary$sd,
                                         '/')

      # Take maximum statistic for each bootstrap iteration and add it to the matrix
      # bootstrap_evidence
      bootstrap_evidence[i,]=apply(bootstrap_stats_specific_dim,2,max)

      if(!is.null(loud)){print(paste('Bootstrap Done on Dimension',i))}
    }

  }else{
    cl <- parallel::makeCluster(num_cores,outfile="")
    doParallel::registerDoParallel(cl)

    bootstrap_evidence=foreach::foreach(i=seq(1,p),
                               .packages = c("Matrix",'plyr'),
                               .combine='rbind') %dopar% {

                                 # get the filter matrix associated with this dimension
                                 filter_matrix_specific_dim_indexes=((filter_matrix_dictionary$scale_count%in%scale_count_list[[i]])
                                                                     &(filter_matrix_dictionary$t_index%in%t_index_list[[i]]))

                                 filter_matrix_specific_dim=filter_matrix[filter_matrix_specific_dim_indexes,]
                                 filter_matrix_specific_dim_dictionary=filter_matrix_dictionary[filter_matrix_specific_dim_indexes,]
                                 filter_matrix_specific_dim_dictionary[,'row_order']=1:nrow(filter_matrix_specific_dim_dictionary)

                                 # Calculate the upsilon_j * z_j for each time j and each bootstrap_iteration
                                 upsilon_times_stdnormal=sweep(bootstrap_normals,MARGIN=2,upsilon_matrix[i,],'*')

                                 # Apply the filter matrix to upsilon * z matrix for each time and scale
                                 bootstrap_stats_specific_dim = abs(filter_matrix_specific_dim%*%t(upsilon_times_stdnormal))

                                 # normalization constants for each filter
                                 normalization_constants_specific_dim=normalization_constants[normalization_constants$dimension==i,
                                                                                              c('sd','t_index')]

                                 # match the normalization constant with the correct time in the same order of the filter matrix
                                 # for this specific dimension
                                 filter_matrix_specific_dim_dictionary = merge(filter_matrix_specific_dim_dictionary,normalization_constants_specific_dim,by='t_index')
                                 filter_matrix_specific_dim_dictionary = filter_matrix_specific_dim_dictionary[order(filter_matrix_specific_dim_dictionary$row_order),]

                                 # Apply normalization to all of the bootstrap statistics
                                 bootstrap_stats_specific_dim=sweep(bootstrap_stats_specific_dim,
                                                                    MARGIN=1,
                                                                    filter_matrix_specific_dim_dictionary$sd,
                                                                    '/')

                                 # Take maximum statistic for each bootstrap iteration and add it to the matrix
                                 # bootstrap_evidence
                                 if(!is.null(loud)){print(paste('Bootstrap Done on Dimension',i))}
                                 apply(bootstrap_stats_specific_dim,2,max)
                               }

    parallel::stopCluster(cl)
  }

  # Record the maximums of each bootstrap iteration
  # Record the maximum test statistic
  if(!is.null(output_bootstraps_and_teststat)){
    bootstrap_teststat_list[['bootstraps']]=apply(bootstrap_evidence,2,max)
  }

  # Use a while loop, and in each iteration compare test statistic to critical value produced
  # by bootstrap. Any time a jump is detected, recalculate the bootstrap critical value (only
  # recalculate for dimension where jump detected to save computation), and maximum test statistic (
  # by ignoring values around the detected jump).
  exit_loop=0

  while(exit_loop==0){

    # The maximum evidence we have for each bootstrap iteration
    bootstrap_max_evidence_by_iteration = apply(bootstrap_evidence,2,max)

    # Decide whether to accept or reject null
    # Check to see if we reject the null hypothesis
    max_test_statistic=test_statistics[which.max(test_statistics$test_stat_normalized),]

    if(max_test_statistic$test_stat_normalized>stats::quantile(bootstrap_max_evidence_by_iteration,1-alpha)){

      changepoint_number = changepoint_number+1
      cp_dimension = max_test_statistic$dimension
      names(cp_dimension)='dim'
      cp_scale = max_test_statistic$scale_count
      cp_time = max_test_statistic$t_index

      if(!is.null(loud)){
        print(paste('CP',changepoint_number,'located in dimension',cp_dimension,'approx t index',cp_time,Sys.time()))
      }

      # Record jump
      changepoints_identified[[changepoint_number]]=c(cp_dimension, cp_scale,cp_time)
      names(changepoints_identified[[changepoint_number]])=c('dim','scale_count','t_index')

      # Set test statistics in neighbourhood of max(s) around this time index to
      # 0 so we do not detect another jump that is too close
      deletion_window_size=ceiling((1+cp_delete_c)*max(scale_count_list[[cp_dimension]]))
      t_index_to_delete_test_stat=seq(cp_time-deletion_window_size,
                                      cp_time+deletion_window_size)

      test_statistics[(test_statistics$dimension==cp_dimension)&
                        (test_statistics$t_index%in%t_index_to_delete_test_stat),
                      'test_stat_normalized']=0

      # Recalculate the bootstrap statistics for the dimension with the jump, what
      # this entails is getting the bootstrap statistics for this specific dimension but
      # excluding all of the statistics near jump that have already been identified
      # (note this code is nearly identical to the code that calculates the bootstrap
      # statistics initially)
      filter_matrix_specific_dim_indexes=((filter_matrix_dictionary$scale_count%in%scale_count_list[[cp_dimension]])
                                          &(filter_matrix_dictionary$t_index%in%t_index_list[[cp_dimension]]))

      filter_matrix_specific_dim=filter_matrix[filter_matrix_specific_dim_indexes,]
      filter_matrix_specific_dim_dictionary=filter_matrix_dictionary[filter_matrix_specific_dim_indexes,]
      filter_matrix_specific_dim_dictionary[,'row_order']=1:nrow(filter_matrix_specific_dim_dictionary)

      # Calculate the upsilon_j * z_j for each time j and each bootstrap_iteration
      upsilon_times_stdnormal=sweep(bootstrap_normals,MARGIN=2,upsilon_matrix[cp_dimension,],'*')

      # Apply the filter matrix to upsilon * z matrix for each time and scale
      bootstrap_stats_specific_dim = abs(filter_matrix_specific_dim%*%t(upsilon_times_stdnormal))

      # normalization constants for each filter
      normalization_constants_specific_dim=normalization_constants[normalization_constants$dimension==cp_dimension,
                                                                   c('sd','t_index')]

      # match the normalization constant with the correct time in the same order of the filter matrix
      # for this specific dimension
      filter_matrix_specific_dim_dictionary = merge(filter_matrix_specific_dim_dictionary,normalization_constants_specific_dim,by='t_index')
      filter_matrix_specific_dim_dictionary = filter_matrix_specific_dim_dictionary[order(filter_matrix_specific_dim_dictionary$row_order),]

      # Apply normalization to all of the bootstrap statistics
      bootstrap_stats_specific_dim=sweep(bootstrap_stats_specific_dim,
                                         MARGIN=1,
                                         filter_matrix_specific_dim_dictionary$sd,
                                         '/')

      # Determine what times to delete (since we already identified them as jumps
      # and don't want to consider small neighbourhood around them)
      cp_already_identified_in_dim=c()

      for(i in changepoints_identified){
        if(i['dim']==cp_dimension){
          cp_already_identified_in_dim=c(cp_already_identified_in_dim,i['t_index'])
        }
      }

      t_index_to_delete_all_cp=c()
      for(i in cp_already_identified_in_dim){
        t_index_to_delete_all_cp=c(t_index_to_delete_all_cp,seq(i-deletion_window_size,
                                                                i+deletion_window_size))
      }

      # Determine which rows of bootstrap_stats_specific_dim correspond to these times, and
      # set their bootstrap stats to 0
      t_index_to_delete_filter_matrix_rows=as.numeric((filter_matrix_specific_dim_dictionary$t_index%in%t_index_to_delete_all_cp)==FALSE)

      bootstrap_stats_specific_dim=sweep(bootstrap_stats_specific_dim,
                                         MARGIN=1,
                                         t_index_to_delete_filter_matrix_rows,
                                         '*')

      # Take maximum statistic for each bootstrap iteration and add it to the matrix
      # bootstrap_evidence
      bootstrap_evidence[cp_dimension,]=apply(bootstrap_stats_specific_dim,2,max)

    }else{
      exit_loop=1
    }

  }

  # SECOND STAGE REFINEMENT
  second_stage_changepoints = second_stage_refine(changepoints_identified,
                                                  x_matrix,scale_count_list,ss_refine_alpha_tilde)
  ss_cp_df = data.frame(dim=character(),
                        scale_count=character(),
                        t_index=character(),
                        stringsAsFactors=FALSE)

  if(length(second_stage_changepoints)>=1){
    for(i in seq(1,length(second_stage_changepoints))){ss_cp_df[i,]=second_stage_changepoints[[i]]}
  }

  ss_cp_df[,'t_index']=as.numeric(ss_cp_df[,'t_index'])

  # If we only want the bootstraps and test statistics return that, otherwise
  # return the changepoints detected via the second stage refinement method
  if(!is.null(output_bootstraps_and_teststat)){
    return(bootstrap_teststat_list)
  }else{
    return(ss_cp_df)
  }

}



