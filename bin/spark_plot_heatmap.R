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



# usage
usage = "Rscript plot_heatmap.R [options]\n"
# description
description = "This program allows to create a nice heatmap of a given data matrix
according to a SPar-K partitioning.
Because SPar-K can shift and flip the data during the partitioning procedure, the 
data are realigned before being displayed. The the data are ordered by cluster 
indicated by a color ribbon on the left. The rows within each cluster are sorted 
by correlation to the cluster aggregation (the most ressembling at the top and the 
most dissimilar at the bottom). Finally the heatmap is compressed to make the 
image nicer and easier to handle by averaging the positions over the rows (for a  
given column, each 20 rows values are replaced by their mean)\n"
# prologue
epilogue = "Written by Romain Groux, November 2018\n"
# options
option_list = list(
  make_option(c("--data"), action="store", default=NULL, type="character",
              help="The file containing the data which have been partitioned using SPar-K."),
  make_option(c("--partition"), action="store", default=NULL, type='character',
              help="The file containing the results of the SPar-K partitioning (as is)."),
  make_option(c("--shift"), action="store", default=1, type="numeric",
              help="The shifting freedom allowed when running SPar-K"),
  make_option(c("--title"), action="store", default="", type="character",
              help="A title to display atop of the plot."),
  make_option(c("--from"), action="store", default=NULL, type="numeric",
              help="The minimum possible value on the x-axis."),
  make_option(c("--to"), action="store", default=NULL, type="numeric",
              help="The maximum possible value on the x-axis."),
  make_option(c("--dim"), action="store", default="8,10", type="character",
              help="The coma separated image dimensions (width,height), in inches, by default 8,10."),
  make_option(c("--res"), action="store", default=NULL, type="numeric",
              help="The image resolution in ppi."),
  make_option(c("--output"), action="store", default="myplot.png", type="character",
              help="The file were the image will be saved in png format.")
)

# parses options
opt = parse_args(OptionParser(usage=usage,
                              option_list=option_list,
                              description=description,
                              epilogue=epilogue))

# parse dimensions
# if not specified, image is 8x10 inches, has 72 ppi of resolution and a pointsize of 12
opt$dim   = as.numeric(unlist(strsplit(opt$dim, split=",")))
opt$res   = ifelse(is.null(opt$res), 72, opt$res)
img.param = list(width=opt$dim[1], height=opt$dim[2], res=opt$res, pointsize=12)

# setwd("/local/groux/Kmeans_chipseq")
# opt = list()
# opt$data = "data_8.txt"
# opt$partition = "results_8.txt"
# opt$shift = 71
# opt$from = -1000
# opt$to = 1000
# opt$title = "TSS with H3K4me3"
# opt$dim = c(8,10)
# opt$res = 72
# img.param = list(width=opt$dim[1], height=opt$dim[2], res=opt$res, pointsize=12)

# check options
if(is.null(opt$data))
{ stop("Error! no data given (--data)") }
if(is.null(opt$partition))
{ stop("Error! no partition given (--partition)") }
if(!file.exists(opt$data))
{ stop("Error! the given data file does not exist (--data)") }
if(!file.exists(opt$partition))
{ stop("Error! the given partition file does not exist (--partition)") }
if(is.null(opt$from))
{ stop("Error! invalid minimum x-axis value (--from)") }
if(is.null(opt$to))
{ stop("Error! invalid maximum x-axis value (--to)") }
if(is.na(img.param$width) || is.na(img.param$height) || 
  (img.param$width <= 0)  || (img.param$height <= 0))
{ stop("Error! invalid image dimensions (--dim)") }
if(img.param$res <= 0)
{ stop("Error! invalid image resolution (--res)") }

# heatmap colors
color = colorRampPalette(c("white",  "red"),  space = "rgb")(100)
# cluster colors
color.lab = brewer.pal(8, "Set1")

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

# number of clusters
n.cluster = length(unique(partition$cluster))


# realign the data and by cluster
d             = matrix(nrow=nrow(data), ncol=ncol(data))
p             = d
data.aligned  = d
labels        = vector(mode="character", length=nrow(data))
from = 1; to = from ;
for(j in 1:n.cluster)
{ index                   = which(partition$cluster == j)
  to                      = from + length(index) -1
  d[from:to,]             = order.rows(data[index,])
  data.aligned[from:to,]  = realign.data(data[index,, drop=F], 
                                         partition$shift_ref[index], 
                                         partition$shift_dat[index], 
                                         partition$flip[index], opt$shift)
  order                   = get.row.order(data.aligned[from:to,, drop=F])
  data.aligned[from:to,]  = data.aligned[from:to,, drop=F][order,]
  labels[from:to]  = color.lab[j]
  from = to + 1
}

# compute labels for x-axis
x.at = seq(from=0,         to=1,      length.out=5)
x.lab = seq(from=opt$from, to=opt$to, length.out=5)

# plot
png(filename=opt$output, width=img.param$width, height=img.param$height, 
    units="in", res=img.param$res, pointsize=img.param$pointsize)
  # layout construction
  lab = c(1, 2)
  lay = layout(matrix(data=lab, nrow=1, ncol=2, byrow=T), widths=c(1,10))
  layout.show(lay)
  # custer labels
  p = par(mar=c(5, 0, 4, 1) + 0.1, oma=c(0,0,0,2))
  plot.label.bar(partition$cluster, lwd=1, colors=labels)
  # heatmap
  p = par(mar=c(5, 0, 4, 1) + 0.1)
  image(t(condense.matrix(data.aligned)), col=color, xaxt='n', yaxt='n', 
        main=opt$title, ylab="", xlab="Approx. pos. (bp)",
        cex.main=2.5, cex.lab=2.5)
  axis(side=1, at=x.at, labels=x.lab, cex.axis=2.5, cex.lab=2.5)
dev.off()


