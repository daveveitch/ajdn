#' 2023 Turkey-Syria Earthquake Data
#'
#' A subset of data from seismic sensors during the 2023 Turkey-Syria earthquake.
#' Data from the ORFEUS Data Center WebDC3 Web Interface http://orfeus-eu.org/webdc3/.
#' The data is centred and the log of the absolute value is taken.
#'
#' @format ## `earthquake_data`
#' A data frame with 28 rows and 5760 columns:
#' \describe{
#'   \item{row names}{The name of the sensor the observations are taken from}
#'   \item{columns}{Data sampled once every 10 seconds.}
#'   ...
#' }
#' @source <http://orfeus-eu.org/webdc3/>
"earthquake_data"

#' 2023 Turkey-Syria Earthquake Meta-Data
#'
#' Metadata on each of the sensors used to collect information during the 2023 Turkey-Syria earthquake.
#'
#' @format ## `earthquake_meta`
#' A data frame with 28 rows and 29 columns:
#' \describe{
#'   \item{sensor}{Short form of the sensor name}
#'   \item{maxgapseconds}{Maximum gap of downtime in seconds.}
#'   \item{numtraces}{number of traces in seismic data (1 trace = 0 downtime)}
#'   \item{starttime}{starttime of observations}
#'   \item{endtime}{endtime of observations}
#'   \item{filename}{name of file downloaded from ORFEUS Data Center WebDC3 Web Interfact}
#'   \item{latitude}{latitude of the sensor}
#'   \item{longitude}{longitude of the sensor}
#'   \item{location}{location name where the sensor is located}
#'   ...
#' }
#' @source <http://orfeus-eu.org/webdc3/>
"earthquake_meta"

#' Stocks Data
#'
#' Stocks data from 2021 to 2022 for META, AMZN, AAPL, NFLX, GOOG. The stock prices have been divided
#' by the value of the QQQ ETF to give a time series of the relative performance of the individual stocks
#' relative to the index.
#'
#' @format ## `stocks_data`
#' A data frame with 5 rows and 503 columns:
#' \describe{
#'   \item{row names}{ticker symbol of stock}
#'   \item{columns}{each column represents the daily price of a stock}
#'   ...
#' }
#' @source <https://www.quantmod.com/whatsnext/>
"stocks_data"
