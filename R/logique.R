logique <- function(x)
 {

   if (!is.vector(x)) stop("Only vectors handled currently as input.\n")

   log.x <- cbind(x, x)
   log.x[x >= median(x), 1] <- 1   
   log.x[x <  median(x), 1] <- 0
   log.x[x >= median(x), 2] <- 0   
   log.x[x <  median(x), 2] <- 1

   log.x

 }

