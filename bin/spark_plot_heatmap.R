#!/usr/bin/Rscript

library(optparse,     quietly=T)
library(RColorBrewer, quietly=T)

# get script dir
initial.options = commandArgs(trailingOnly = FALSE)
file.arg.name = "--file="
script.name = sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dir = dirname(script.name)

# functions
source(file.path(script.dir, "spark_functions.R"))

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


