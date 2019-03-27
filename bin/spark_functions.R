#!/usr/bin/Rscript

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

#' Order the rows of a given matrix by similarity 
#' (correlation) to the aggregation (in descending 
#' order) and returns the order.
#' @param data the matrix of interest.
#' @return a vector of indices to reorder the 
#' original matrix.
#' @author Romain Groux
get.row.order = function(data)
{ if(is.vector(data))
  { return(c(1)) }
  else
  { ref    = colSums(data)
    scores = apply(data, 1, cor, ref)
    return(order(scores, decreasing=F))
  }
}

#' A function to condense a matrix. Every <rows> rows are averaged into 
#' one row.
#' @param MATRIX a matrix [m,n] of numerics of interest.
#' @param rows the number of rows which will be used for averaging
#' 20 by default.
#' @return a matrix [m',n] of numerics  where m' = ceiling(m/<rows>).
#' @author Romain Groux
condense.matrix = function(MATRIX, rows=20)
{ n.row = nrow(MATRIX)
  n.col = ncol(MATRIX)
  groups      = sort(rep(1:ceiling(n.row/rows), rows))[1:n.row]
  MATRIX.OUT  = matrix(ncol=n.col, nrow=max(groups))
  for(i in 1:n.col)
  { MATRIX.OUT[,i] = tapply(MATRIX[,i], groups, mean, na.rm=T) }
  return(MATRIX.OUT)
}

#' Create a label bar in a plot containing as many 
#' colored horizontal lines as they are elements in 
#' <labels>. Colors are choosen such that all the 
#' elements with a given value in <labels> have the
#' same color.
#' @param labels a vector containing labels to be 
#' displayed in a color format.
#' @param lwd the width of the lines.
#' @param colors a vector containing a color per 
#' <labels> element. The order of <colors> and 
#' <labels> is expected to be consistant.
#' @param each a number indicating that only one value every <each>
#' in <labels> should be plotted. Colors from <colors> will be 
#' used accordingly.
#' @return nothing.
#' @author Romain Groux
plot.label.bar = function(labels, lwd=1, colors=NULL, each=NULL)
{ # need to compress the labels vectors
  # only the 1st, 1+each th, 1+each+each th values will be plotted
  # if some values remains after the last step, they are not considered
  if(!is.null(each))
  { labels.tmp = vector(mode="numeric", length=ceiling(length(labels)/each))
    colors.tmp = labels.tmp
    i = j = 1
    while(j<length(labels))
    { labels.tmp[i] = labels[j]
      if(!is.null(colors))
      { colors.tmp[i] = colors[j] }
      i = i + 1
      j = j + each
    }
    labels = labels.tmp
    labels.tmp = NULL
    if(!is.null(colors))
    { colors = colors.tmp
      colors.tmp = NULL
    }
  }
  n = length(labels)
  plot(x=c(1:n), y=rep(1,n), xlim=c(0,1), ylim=c(1,n),
       # xaxs and yaxs expend the plotting area such that it stiks to the axis
       xaxt='n', yaxt='n', xaxs='i', yaxs='i',
       main='', ylab='', xlab='', xaxs='i', yaxs='i', col="white")
  if(is.null(colors))
  { colors = labels }
  segments(x0=rep(0,n), y0=c(1:n), x1=rep(1,n), y1=c(1:n), col=colors, lwd=lwd)
  
}

#' A function to reorder the rows of a matrix according to the overall
#' similarity (correlation) of each row to the column sums of the matrix. 
#' The top rows will have the highest similarity with the column sums and
#' the bottom rows the lowest.
#' @param MATRIX a matrix (m x n) of numerics of interest.
#' @return the reordered matrix (m x n) of numerics.
#' @author Romain Groux
order.rows = function(MATRIX)
{ # for some reason, MATRIX has only 1 row -> is a vector
  if(is.vector(MATRIX))
  { return(MATRIX) }
  else
  { REF      = colSums(MATRIX)
    SCORES   = apply(MATRIX, 1, cor, REF)
    ord      = order(SCORES, decreasing=F)
    return(MATRIX[ord,])
  }
}

#' @brief This function corrects the coordinates in a SGA file given some 
#' formated partitioning results, as returned by format.results().
#' The SGA coordinate correction works as follow : for each SGA coordinate, a 
#' correction offset is computed based on the corresponding shift and flip 
#' values.
#' The alignment procedure alignes a slice of each row of the data matrix 
#' against its reference. One of the bins of each row of the data matrix 
#' contains the position at which the original SGA coordinate was located.
#' Over the entire matrix, this is a given column. Once all the slices 
#' have been aligned to the references, these positions are shifted away 
#' from their original locations.
#' 
#'                       X                           orginal SGA coord in bin 6
#' | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |   row
#' | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |             shift 1 slice (shift_dat = 1)
#'     | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |        shift 2 slice (shift_dat = 2)
#'         | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |   shift 3 slice (shift_dat = 3)
#'
#' Let's align each possible slice WITHOUT FLIP at difference shift_ref along 
#' a reference.
#'                       X                           orginal SGA coord in bin 6
#' | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |   reference
#'                       X                           shift_dat    shift_ref
#' | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |                 1            1
#'                           X
#'         | 2 | 3 | 4 | 5 | 6 | 7 | 8 |  9 | 10 |       2            3
#'                   X
#'     | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11 |            3            2
#' 
#' The orignal SGA location in each slice corresponds to bin 6. We now want 
#' to align all bins corresponding to the 6th bin to the 6th bin in the reference. 
#' Thus, for the first row, no change is required, the 6th bin is aligned with 
#' the 6th bin of the reference. For the second row, the 6th bin is aligned with 
#' the 7th bin in the reference. So the SGA position should be pushed to the left 
#' by 1 which means a correction of -1 such that the corrected SGA position lies 
#' on the 5th bin. For the third row, this is the opposite. The 6th bin is 
#' aligned on the reference 5th bin which means that the SGA position should 
#' be moved to the right by one so corrected by +1.
#' 
#'                       X                           orginal SGA coord in bin 6
#' | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |   reference
#'                       X                           shift_dat    shift_ref
#' | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |                 1            1
#'                       X<--X
#'         | 2 | 3 | 4 | 5 | 6 | 7 | 8 |  9 | 10 |       2            3
#'                   X-->X
#'     | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11 |            3            2
#'     
#'     
#'To realign slices WITH FLIP, the principle is almost the same.
#'                       X                           orginal SGA coord in bin 6
#' | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |   reference
#'               X                                  shift_dat    shift_ref
#' | 9 | 8 | 7 | 6 | 5 | 4 | 3 | 2 | 1 |                 1 (rev)     1
#'                           X
#'         |10 | 9 | 8 | 7 | 6 | 5 | 4 | 3 |  2 |        2 (rev)     3
#'                           X
#'     |11 |10 | 9 | 8 | 7 | 6 | 5 | 4  |  3 |           3 (rev)     2
#'
#' Here each slice is aligned to the reference in the but in reverse.
#' The first slice 6th bin is aligned to the 4th reference bin. Thus 
#' it should be moved by 2 to the right (visually). Which means that 
#' the SGA position should be moved from the 6th bin to the 4th bin, 
#' corresponding to a correction of -2. For the slices 2 and 3, 
#' the SGA position should be moved by 1 to the left (visually) so 
#' the applied correction should be +1.
#' 
#'                       X                           orginal SGA coord in bin 6
#' | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |   reference
#'               X------>X                           shift_dat    shift_ref
#' | 9 | 8 | 7 | 6 | 5 | 4 | 3 | 2 | 1 |                 1 (rev)     1
#'                       X<--X
#'         |10 | 9 | 8 | 7 | 6 | 5 | 4 | 3 |  2 |        2 (rev)     3
#'                       X<--X
#'     |11 |10 | 9 | 8 | 7 | 6 | 5 | 4  |  3 |           3 (rev)     2
#' 
#' @param sga a dataframe containing the sga file to correct.
#' @param results the formatted partitioning results, as returned by format.results().
#' @param shift the shifting freedom allowed during the partitioning.
#' @param data.ncol the number of columns of the partitioned data matrix.
#' @param bin.size the number of base pairs which is covered by each bin in the 
#' data matrix which was partitioned.
#' @seealso format.results()
#' @autor Romain Groux
realign.sga = function(sga, results, shift, data.ncol, bin.size)
{
  n.row        = nrow(sga)
  n.col        = data.ncol
  # l.slice      = n.col - shift  + 1
  if(is.null(results$flips))
  { results$flips = rep(FALSE, n.row) }
  
  # 0-base to 1-based
  results$shift_dat = results$shift_dat + 1
  results$shift_ref = results$shift_ref + 1
  
  # update position according to shifts
  # not flipped
  index   = which(results$flip == FALSE)
  if(length(index))
  { sga[index,3] = sga[index,3] + (results$shift_dat[index] - results$shift_ref[index])*bin.size }
  # flipped
  index   = which(results$flip == TRUE)
  shift.central    = ceiling(shift / 2)                      # central shift state
  to.shift.central = shift.central - results$shift_ref      # diff to central shift for ref shift
  shift.ref.rev    = results$shift_ref + 2*to.shift.central # shift ref corresponding to the other side
  # of the ref
  if(length(index))
  { sga[index,3] = sga[index,3] + (results$shift_dat[index] - shift.ref.rev[index])*bin.size }
  
  # update strand according to flip
  sga[,4] = as.character(sga[,4])
  for(i in 1:n.row)
  { # sga is strandless
    if(sga[i,4] == "0")
    { if(results$flip[i])
      { sga[i,4] = "-" }
      else
      { sga[i,4] = "+" }
    }
    # sga has strand info
    else
    { if(results$flip[i] && sga[i,4] == "+")
      { sga[i,4] = "-" }
      else if(results$flip[i] && sga[i,4] == "-")
      { sga[i,4] = "+" }
    }
  }
  return(sga)
}

