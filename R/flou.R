flou <- function(x)
 {

   if (!is.vector(x)) stop("Only vectors handled currently as input.\n")

   log.x <- cbind(x, x)

   lo    <- quantile(x, .33)
   hi    <- quantile(x, .67)

   log.x[x >= hi, 1] <- 1   
   log.x[x >= hi, 2] <- 0   

   log.x[x <= lo, 1] <- 0
   log.x[x <= lo, 2] <- 1

   log.x[x > lo & x < hi, 1] <- (x[x > lo & x < hi] - lo) / (hi - lo)
   log.x[x > lo & x < hi, 2] <- (hi - x[x > lo & x < hi]) / (hi - lo)

   log.x

 }

