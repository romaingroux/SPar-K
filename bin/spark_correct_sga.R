#!/usr/bin/Rscript

library(optparse, quietly=T)

# get script dir
initial.options = commandArgs(trailingOnly = FALSE)
file.arg.name = "--file="
script.name = sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dir = file.path(getwd(), dirname(script.name))

# functions
print(script.dir)
source(file.path(script.dir, "spark_functions.R"))

# usage
usage = "Rscript make_sga.R [options]\n"
# description
description = "This program allows to correct the coordinates contained in a SGA file 
according to the results of a SPar-K partition.
The SGA coordinate correction works as follow : for each SGA coordinate, a 
correction offset is computed based on the corresponding shift and flip 
values.
The alignment procedure alignes a slice of each row of the data matrix 
against its reference. One of the bins of each row of the data matrix 
contains the position at which the original SGA coordinate was located.
Over the entire matrix, this is a given column. Once all the slices 
have been aligned to the references, these positions are shifted away 
from their original locations.

                      X                           orginal SGA coord in bin 6
| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |   row
| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |             shift 1 slice (shift_dat = 1)
    | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |        shift 2 slice (shift_dat = 2)
        | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |   shift 3 slice (shift_dat = 3)

Let's align each possible slice WITHOUT FLIP at difference shift_ref along 
a reference.
                      X                           orginal SGA coord in bin 6
| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |   reference
                      X                           shift_dat    shift_ref
| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |                 1            1
                          X
        | 2 | 3 | 4 | 5 | 6 | 7 | 8 |  9 | 10 |       2            3
                  X
    | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11 |            3            2

The orignal SGA location in each slice corresponds to bin 6. We now want 
to align all bins corresponding to the 6th bin to the 6th bin in the reference. 
Thus, for the first row, no change is required, the 6th bin is aligned with 
the 6th bin of the reference. For the second row, the 6th bin is aligned with 
the 7th bin in the reference. So the SGA position should be pushed to the left 
by 1 which means a correction of -1 such that the corrected SGA position lies 
on the 5th bin. For the third row, this is the opposite. The 6th bin is 
aligned on the reference 5th bin which means that the SGA position should 
be moved to the right by one so corrected by +1.

                      X                           orginal SGA coord in bin 6
| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |   reference
                      X                           shift_dat    shift_ref
| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |                 1            1
                      X<--X
        | 2 | 3 | 4 | 5 | 6 | 7 | 8 |  9 | 10 |       2            3
                  X-->X
    | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10| 11 |            3            2
    
    
To realign slices WITH FLIP, the principle is almost the same.
                      X                           orginal SGA coord in bin 6
| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |   reference
              X                                  shift_dat    shift_ref
| 9 | 8 | 7 | 6 | 5 | 4 | 3 | 2 | 1 |                 1 (rev)     1
                          X
        |10 | 9 | 8 | 7 | 6 | 5 | 4 | 3 |  2 |        2 (rev)     3
                          X
    |11 |10 | 9 | 8 | 7 | 6 | 5 | 4  |  3 |           3 (rev)     2

Here each slice is aligned to the reference in the but in reverse.
The first slice 6th bin is aligned to the 4th reference bin. Thus 
it should be moved by 2 to the right (visually). Which means that 
the SGA position should be moved from the 6th bin to the 4th bin, 
corresponding to a correction of -2. For the slices 2 and 3, 
the SGA position should be moved by 1 to the left (visually) so 
the applied correction should be +1.

                      X                           orginal SGA coord in bin 6
| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |   reference
              X------>X                           shift_dat    shift_ref
| 9 | 8 | 7 | 6 | 5 | 4 | 3 | 2 | 1 |                 1 (rev)     1
                      X<--X
        |10 | 9 | 8 | 7 | 6 | 5 | 4 | 3 |  2 |        2 (rev)     3
                      X<--X
    |11 |10 | 9 | 8 | 7 | 6 | 5 | 4  |  3 |           3 (rev)     2\n"
# prologue
epilogue = "Written by Romain Groux, November 2018\n"
# options
option_list = list(
  make_option(c("--sga"), action="store", default=NULL, type='character',
              help="A SGA file to update according to the given partitioning. The 
              SGA coordinates and strand will be corrected according to the shifting 
              and flipping values assigned by SPar-K."),
  make_option(c("--partition"), action="store", default=NULL, type='character',
              help="A file containing the result of a partitioning done by SPar-K. 
              This partioning should be linked to the SGA file."),
  make_option(c("--cluster"), action="store", default=0, type='numeric',
              help="Allows to correct the regions assigned to a given cluster only. By 
              default (--cluster=0) all the regions are corrected."),
  make_option(c("--shift"), action="store", default=1, type='numeric',
              help="The shifting freedom allowed when partitioning the data with 
              SPar-K."),
  make_option(c("--ncol"), action="store", default=0, type='numeric',
              help="The number of columns in the original data matrix partitioned 
              by SPar-K."),
  make_option(c("--binSize"), action="store", default=0, type='numeric',
              help="The number of base pairs contained within each column of the 
              original data matrix partitioned by SPar-K (the bin size).")
)

# parses options
opt = parse_args(OptionParser(usage=usage,
                              option_list=option_list,
                              description=description,
                              epilogue=epilogue))

# check options
if(is.null(opt$partition))
{ stop("Error! no partition given (--partition)") }
if(is.null(opt$sga))
{ stop("Error! no SGA file given (--sga)") }
if(!file.exists(opt$partition))
{ stop("Error! the given partition file does not exist (--partition)") }
if(!file.exists(opt$partition))
{ stop("Error! the given SGA file does not exist (--sga)") }
if(opt$cluster < 0)
{ stop(sprintf("Error! Invalid cluster label given (--cluster) : %d", opt$cluster)) }
if(opt$shift <= 0)
{ stop(sprintf("Error! Invalid shift given (--shift) : %d", opt$shift)) }
if(opt$binSize <= 0)
{ stop(sprintf("Error! Invalid bin size given (--binSize) : %d", opt$binSize)) }

# read partition
partition = read.table(opt$partition, header=T, stringsAsFactors=F)
# read SGA
sga = read.table(opt$sga,             header=F, stringsAsFactors=F)

# check that partitioning results and SGA file have the same number of rows
if(nrow(sga) != nrow(partition))
{ stop(sprintf("Error! the partition and the SGA file don't have the same number of rows : %d / %d",
               nrow(partition), nrow(sga)))
}

# realign
if(opt$cluster != 0)
{ index = which(partition$cluster == opt$cluster)
} else { 
  index = 1:nrow(sga)
}
sga  = realign.sga(sga[index,],
                   partition[index,], 
                   opt$shift,
                   opt$data_ncol,
                   opt$binSize)

# print
write.table(x=sga,
            file="",
            quote=FALSE,
            sep='\t',
            eol='\n',
            row.names=F,
            col.names=F)
