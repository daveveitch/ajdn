#' Create W filter
#'
#' @description Creates a vector of the optimal jump pass filter for a n (length of time series), s (scale), and t (time)
#'
#' @param n length of full time series
#' @param scale_count rounded value of n*s where s is the scale
#' @param t_index the index associated with the time at the center of the filter
#' @return a vector of length n with weights that the filter takes
#' @export
#'
#' @examples
#' w = create_W_filter(1000,100,300)
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
#' @param x_matrix p x n (p = # dimensions, n = # obsrevations) matrix of time series
#' @param s_prime block size parameter
#' @return p x n matrix to be subsequently used in the block multiplier bootstrap
#' @export
#'
#' @examples
#' x_matrix = matrix(rnorm(1000),nrow=10)
#' generate_upsilon_matrix(x_matrix,1/100)
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
#' t_index (time where jump found), scale_count (ns where s scale used to find jump)
#' e.g. first_stage_changepoints[[5]]['t_index'] is the time index the 5th jump identified was found at
#' @param x_matrix p x n (p = # dimensions, n = # obsrevations) matrix of time series
#' @param scale_count_list a list of length p, where the pth entry of the list is a vector of \eqn{ns_{r,1},...,ns_{r,\delta_n}} where
#' \eqn{\delta_n} is the number of scales being used in dimension r
#' @return a list of same length as first_stage_changepoints with the refined estimate of the time index where the jump occurs
#' @export
#'
second_stage_refine <-function(first_stage_changepoints,x_matrix,scale_count_list){
  # This function takes a list of jumps that were already identified and
  # applies a second stage procedure to enhance the estimation accuracy of the procedure
  # INPUT
  # first_stage_changepoints - the jumps identified during the first stage, a list with
  #                            dimension, scale, and time of each
  # x_matrix                 - matrix of x values
  # scale_count_list         - a list of scales (in count form) for each dimension

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


