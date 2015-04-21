#' Cleaned up version of the Mayo PBC data.
#'
#' A cleaned up version of the Mayo PBC data from \code{survival::pbc}. Only the
#' first 312 observations are used (the cases who participated in the
#' randomized trial). Only five of the covariates (listed below) are used.
#' Further, two of the covariates were log transformed.
#'
#' @format A data frame with 312 rows and 7 variables:
#' \describe{
#'   \item{time}{follow up time in days}
#'   \item{status}{1=death, 0=censored}
#'   \item{age}{age in years}
#'   \item{edema}{0=no edema, 0.5=untreated or successfully treated, 1=edema
#'     despite diuretic therapy}
#'   \item{bili}{log serum bilirubin level (original value from
#'     \code{survival::pbc} is unlogged)}
#'   \item{albumin}{serum albumin}
#'   \item{protime}{log standardized blood clotting time (original value from
#'     \code{survival::pbc} is unlogged)}
#' }
#' @source Cleaned up version of \code{survival::pbc}
"pbc5"
