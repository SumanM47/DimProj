#' @name prepDat
#' @title prepare the data for analysis
#'
#' @description
#' Prepares the data from different slices into the input format for the analysis function Slices_to_3D
#'
#' @usage prepDat(data,bma,gridsize)
#'
#' @param data list of ppp objects. Each entry of the list is a ppp object corresponding to a slice. Instead of a ppp object, it can be a named list with the following entries: x(numeric vector of x-coordinates of the locations of the observations), y(numeric vector of y-coordinates of the locations of the observations), and  marks(character vector of type or identifier of the observations that classifies them into the different possible species).
#' @param bma positive scalar. Thickness of each slice
#' @param gridsize numeric vector of size two with the entries being positive integers themselves. The desired dimensions of the computing grid.
#'
#' @import stats
#' @import grDevices
#'
#' @return list of data products needed for input to Slices_to_3D. Object of class X
#' @export

prepDat <- function(data,bma,
                    gridsize){
  I <- length(data)
  M <- gridsize[1]
  N <- gridsize[2]
  n <- M*N

  allmarks <- Sall <- NULL
  for(i in 1:I){
    allmarks <- c(allmarks,data[[i]]$marks)
    Sall <- rbind(Sall,cbind(data[[i]]$x,data[[i]]$y))
  }
  umarks <- unique(allmarks)
  J <- length(umarks)

  Y <- array(0,dim=c(I,J,n))
  S <- matrix(NA,n,2)
  xmin <- min(Sall[,1]); xmax <- max(Sall[,1])
  ymin <- min(Sall[,2]); ymax <- max(Sall[,2])

  dx <- (xmax - xmin)/(M-1); dy <- (ymax-ymin)/(N-1)
  xls <- xmin -0.5*dx + (1:M)*dx; xus <- xls + dx
  yls <- ymin -0.5*dx + (1:N)*dy; yus <- yls + dy

  for(i in 1:I){
    dat1 <- data[[i]]
    for(j in 1:J){
      eligind <- which(dat1$marks==umarks[j])
      if(length(eligind)>0){
      x <- dat1$x[eligind]
      y <- dat1$y[eligind]
      for(gi in 1:M){
        for(gj in 1:N){
          Y[i,j,(gj-1)*M+gi] <- length(intersect(which(x >= xls[gi] & x <= xus[gi]),which(y >= yls[gj] & y <= yus[gj])))
        }
      }
      }
    }
  }
    for(gi in 1:M){
      for(gj in 1:N){
        S[(gj-1)*M+gi,] <- c(xls[gi]+ 0.5*dx,yls[gj] + 0.5*dy)
      }
    }


  ll <- list("Y" = Y, "S" = S, "bma" = bma, "gridsize" = gridsize)
  return(ll)
}
