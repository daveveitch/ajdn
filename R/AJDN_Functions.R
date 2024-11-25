#' Create W filter
#'
#' @description Creates a vector of the optimal jump pass filter for a specific scale and time.
#'
#' @param n length of full time series
#' @param scale_count rounded value of n*s where s is the scale
#' @param t_index the index associated with the time at the center of the filter
#' @return a vector
#' @export
#'
#' @examples
#' create_W_filter(100,10,50)
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
