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
#' w = create_W_filter(100,10,50)
#' plot(w)
create_W_filter <- function(n,scale_count,t_index){
  # This function creates a filter for a specific n,s,t
  # INPUT
  # scale_count - number of observations for this particular scale
  # t_index     - index associated with the time we are looking at
  # n           - number of total observations in time series
  # OUTPUT
  # W_row       - a vector of length n of the values the filter takes for each index i

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
#' @return p x n matrix to be subsequently used in the block multiplier
#' @export
#'
#' @examples
#' x_matrix = matrix(rnorm(1000),nrow=10)
#' generate_upsilon_matrix(x_matrix,1/100)
generate_upsilon_matrix <- function(x_matrix,s_prime){
  # This function calcuates the values used in the multiplier bootstrap, $\Upsilon_i$
  # INPUT
  # x_matrix       - p x n time series matrix
  # s_prime        - block size being used
  # OUTPUT
  # upsilon_matrix - p x n matrix of values used in the multiplier bootstrap

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

