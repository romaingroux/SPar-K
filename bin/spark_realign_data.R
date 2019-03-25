#!/usr/bin/Rscript

library(optparse,     quietly=T)
library(RColorBrewer, quietly=T)

# functions
#' A function to realign data according to the results of the clustering
#' process.
#' @param data a matrix (m x n) of data to realign.
#' @param shifts.ref a vector (m) of int, corresponding to the shift.ref
#' field return by the clustering process.
#' @param shifts.dat a vector (m) of int, corresponding to the shift.dat
#' field return by the clustering process.
#' @param flips a vector (m) of int (1 or 0), corresponding to the flip to 
#' apply to each row of <DATA>. For instance, the flip values returned when 
#' running the clustering application with flip. By default flips=NULL which 
#' disables data flipping for the realignment.
#' @param n.shift the number of possible shift states, this 
#' is the value of the shift paramter used when running the clustering.
#' @return a matrix (m x n) of numerical corresponding the original data 
#' after realignment.
#' @author Romain Groux
realign.data = function(data, shifts.ref, shifts.dat, flips=NULL, n.shift)
{ 
  n.row        = nrow(data)
  n.col        = ncol(data)
  l.slice      = n.col - n.shift  + 1
  flip         = ifelse(is.null(flips), FALSE, TRUE)
  
  data.aligned = matrix(data=0, nrow=n.row, ncol=n.col)
  
  for(i in 1:n.row)
  { 
    from.ref = shifts.ref[i] + 1
    to.ref   = from.ref + l.slice - 1
    from.dat = shifts.dat[i] + 1
    to.dat   = from.dat + l.slice - 1
    
    if(flip && flips[i])
    { data.aligned[i,from.ref:to.ref] =  rev(data[i,from.dat:to.dat]) }
    else
    { data.aligned[i,from.ref:to.ref] =      data[i,from.dat:to.dat]  }
  }
  return(data.aligned)
}


# usage
usage = "Rscript realign_data.R [options]\n"
# description
description = "This program realigns an input matrix that has been clustered using SPar-K.
As SPar-K can shift and flip the data during the partitioning procedure, this program 
process the data in order to reproduce the alignment performed by SPar-K.\n"

# prologue
epilogue = "Written by Romain Groux, March 2019\n"
# options
option_list = list(
  make_option(c("--data"), action="store", default=NULL, type="character",
              help="The file containing the data which have been partitioned using SPar-K."),
  make_option(c("--partition"), action="store", default=NULL, type='character',
              help="The file containing the results of the SPar-K partitioning (as is)."),
  make_option(c("--shift"), action="store", default=1, type="numeric",
              help="The shifting freedom allowed when running SPar-K")
)

# parses options
opt = parse_args(OptionParser(usage=usage,
                              option_list=option_list,
                              description=description,
                              epilogue=epilogue))

# check options
if(is.null(opt$data))
{ stop("Error! no data given (--data)") }
if(is.null(opt$partition))
{ stop("Error! no partition given (--partition)") }
if(!file.exists(opt$data))
{ stop("Error! the given data file does not exist (--data)") }
if(!file.exists(opt$partition))
{ stop("Error! the given partition file does not exist (--partition)") }

# read data
data      = as.matrix(read.table(opt$data, header=F))
# read partition
partition = read.table(opt$partition, header=T, stringsAsFactors=F)
# if SPar-K was run without flipping, no flipping column, add one
if(is.null(partition$flip))
{ partition$flip = rep(0, nrow(partition)) }

# some extra checks
if(nrow(data) != nrow(partition))
{ stop(sprintf("Error! data and partition don't have the same number of rows : %d / %d",
               nrow(data), nrow(partition)))
}
if(opt$shift >= ncol(data))
{ stop(sprintf("Error! shift >= data column number : %d / %d",
               opt$shift, ncol(data)))
}

# realign the data
data = realign.data(data,
                    partition$shift_ref, 
                    partition$shift_dat,
                    partition$flip, opt$shift)
write.table(data,
            file="",
            quote=F,
            row.names=F,
            col.names=F,
            sep='\t')
