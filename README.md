# SPar-K

SPar-K (Signal Partitioning using K-means) is a modified version of a standard K-means algorithm designed to cluster vectors containing a sequence of signal (that is, the order in which the elements appear in the vectors is meaningful). In order to detect a possible phase shift or orientation inversion between two vectors, this program allows computing distances between two vectors by shifting and flipping them (see below).


## SPar-K partitioning procedure
SPar-K implements a modified version of the K-means algorithm. In brief, it iteratively partitions a set of genomic regions based on their signal profiles by optimizing K gap-free alignments of the signal in the regions.

## Input

The data should be stored as a numerical matrix in a simple text file. Each row of the matrix represents a region and each element of a row represents a position along the region. A given element in the matrix represents the amount of signal present at a given position, in a given region. Each row should be stored within the file as a single line. Each row element should be separated from the others by a blank character (space or tab). Finally, each row should have the same length and the matrix should only contains numerical values. No row nor column name are allowed! Here is an example of a valid input (you can find another one in data.txt) :

```
0 0 0 1 1 2 3 2 1 1 0 0 0
0 1 1 2 3 2 1 1 0 0 0 0 0
0 0 4 4 3 2 2 1 1 0 0 0 0
0 0 0 1 1 2 2 3 3 4 4 0 0
```

This matrix can contain the results of a spatial correlation between two sets of genomic features, for instance ChIP-seq reads (the targets) around +/- 1kb of a set of 1000 TSSs (the references). In that case, the matrix is expected to have 1000 rows (one per TSS) and one column per possible position around these references (here 2001 : 1000 downstream of each TSS, 1 where the TSSs are, 1000 upstream of each TSS). Then, each value of the matrix represents the number of targets (ChIP-seq reads) at a given position (the column) relative to a given reference (TSS). It is also possible to use bins, that is, to summarize several positions within each column, for instance to count the target every 10bp instead of every bp. In this case, each column would represent a bin of 10bp.

### Partitioning procedure

First, an optional data pre-processing step, to smooth out outliers, row-wise, is available. It allows to minimize the impact of extreme values on the correlation computations and to limit their driving force on the data alignment process.
The partition is initialized using a random seeding strategy or using the K-means++ strategy. Each cluster is composed of an alignment of regions (rows) assigned to this cluster and is summarized by the aggregation of the data alignment. The aggregation is a vector having a length equal to the number of columns of the input matrix. It represents the average signal profile over the regions assigned in this cluster. The aggregations are obtained by taking the mean signal at each position (column) in the alignment.
Then, the partition, is iteratively optimized. At each iteration, each region is compared to each cluster aggregation, using a modified correlation distance allowing shifting and flipping. Both parameters are defined by the user. In brief, the aim is to detect a possible signal shift of inversion between the two vectors. With a shifting freedom S, each region and cluster aggregation, both of lengths L, are broken down into S sub-parts of length L − S + 1.  To compare a region to a cluster, each sub-part of a given region is compared to each sub- part of the given cluster aggregation. The comparison with the lowest correlation distance is stored as well as the offsets at which the region and the cluster aggregation sub-parts started. Flipping is handled by allowing an additional set of comparisons with the reversed (flipped) region sub-part. The region is then assigned to the least dissimilar cluster. Eventually, the K alignments have been updated and allow to recompute the cluster aggregations.
This procedure is repeated, optimizing the partition until convergence or until reaching the maximum number of iterations.


### Output
SPar-K returns a table through the stdout. It contains a header row and as many rows as the input had. Each row contains several parameters for the corresponding reference. It contains 1) the cluster assignment, 2) the shift and flip values describing how the row and the corresponding cluster reference were aligned - that is the coordinate of the 1st element of the cluster reference and of the matrix row used in the comparison leading to assigning this row to the given cluster, whether one of the slice was flipped or not - and 3) the distance computed between these two slices. If flipping is not allowed, then no flipping information is returned.


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. To run SPar-K, you have three options. You can either choose to download a release source code, to download a Docker image or a Singularity image. All procedures are detailed below.

### From the release source code
#### Prerequisites

To compile and run SPar-K, the following programs and libraries need to be installed on your computer for a proper compilation :

	1) Scons v2.3.0 or higher to compile all the program listed above(https://scons.org/pages/download.html) 
	2) boost v1.4.1 or higher (https://www.boost.org/)
	3) UnitTest++ v2 (https://github.com/unittest-cpp/unittest-cpp)

The Scons configuration files SConstruct and SConscript are configured such that they will look for :

	1) boost static libaries in /usr/local/lib/boost and boost header files in /usr/local/include/boost
	2) UnitTest++ static libaries in /usr/local/lib/UnitTest++ and UnitTest++ header files in /usr/local/include/UnitTest++

so it is highly recommanded to install these libraries here. Another solution is to modify the SConscript file (located in src/)
to adapt the library paths (modify the lib_unittest_path and lib_boost_path variable values).

The following softwares and libraries are required to run the auxiliary scipts spark_correct_sga.R and spark_plot_heatmap.R :

	1) R version 3.X and Rscript to run these scripts in batch mode
	2) the R libraries optparse and RColorBrewer


#### Compiling

Once all the libraries are installed, download the source, unzip the archive, cd at the root of the repository (where the SConstruct file is located) and compile using Scons:

```
unzip SPar-K-release.zip
cd SPar-K-release
scons
```

The SPar-K exectuable should be located in bin/. To get SPar-K help, run :

```
bin/spark --help
```
To run SPar-K, run :
```
bin/spark <options>
```

#### Running the tests

At compilation, a test suite is also compiled and placed in bin/. To run it and test the different components of the code, use :

```
bin/unittests
```


### From the Singularity image
The Singularity image is build using [the latest release](https://github.com/romaingroux/SPar-K/releases) source code.

#### Prerequisites
You need to have Singularity installed on your machine. Check [this link](https://singularity.lbl.gov/install-linux) for more informations.

#### Pulling the image
Once you have a working version of Singularity, you can pull the image from Singularity Hub using this command
```
singularity pull --name spar-k.simg shub://romaingroux/SPar-K:latest
```

#### Running SPar-K from the image
Using SPar-K from Singularity is just as the same as using the compiled executable, excepted that the commands require to contain a call to Singularity. For instance, to get SPar-K help, use :
```
singularity exec spar-k.simg spark --help
```
To run SPar-K, use :
```
singularity exec spar-k.simg spark <options>
```


### From the Docker image
The Docker image is build using [the latest release](https://github.com/romaingroux/SPar-K/releases) source code.

#### Prerequisites
You need to have Docker installed on your machine. Check [this link](https://www.docker.com/get-started) for more informations. Depending on your installation, you may need root privileges.

#### Pulling the image
Once you have a working version of Docker, you can pull the image from Docker Hub using this command
```
docker pull rgroux/spar-k:latest
```

#### Running SPar-K from the image
Using SPar-K from Docker only requires to deploy a container and call SPar-K. For simplicity, let's tag the image as 'spar-k' (this will be assumed in all the following Docker related documentation). For instance, to get SPar-K help, use :
```
docker tag  rgroux/spar-k:latest spar-k

docker run -i spar-k spark --help
```
You noticed that we called the image by its tag name (spar-k), inside which we ran SPar-K executable (spark).


## Programs
The following three programs are distributed :

### spark
This is the main program. spark is the partitioning software.

#### options
  | short | long&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | description |
  |:------|:--------------------|:---------------------------|
  | -h    | \-\-help            | Produces the help message |
  | -v    | \-\-version         | Prints the version number  |
  | -p    | \-\-parallel  arg   | The number of threads dedicated to the computations, by default 1. |
  | -d    | \-\-data  arg       | The data file address. |
  | -r    | \-\-references  arg | The cluster reference pattern file address. |
  | -i    | \-\-iter  arg       | The maximum number of iterations. |
  | -c    | \-\-cluster  arg    | The number of cluster to find. |
  | -s    | \-\-shift  arg      | Enables this number of column of shifting freedom. By default, shifting is disabled (equivalent to --shift 1). This option and --width are mutually exclusive |
  | -w    | \-\-width           | Enables shifting by searching signal profiles of the given width. Setting --width L' is equivalent to set --shift L-L'+1 where L is the length of each region (the number of columns in the input matrix). By default, the profile width is equal to region width (L). This option and --shift are mutually exclusive.
  |       | \-\-flip            | Enables flipping. |
  |       | \-\-nooutlier       | Pre-pcocess the data to smooth out outliers from the data in a row-wise manner. Each row is searched for outliers which are defined as any value bigger/smaller than the row mean +/- 3*row standard deviation. If a value is an outlier it is replaced by the mean of its left and right neighbours. If a has only a left/right neighbour, then the two left/right neighbours are averaged. If several outliers follow each other the above process is applied to the values in a left to right order and at the end the new averaged values may still be outliers. |
  |       | \-\-dist arg        | Specify which distance should be used during the clustering. It should be 'corr' (by default) or 'normcorr'. |
  |       | \-\-seeding arg     | Specify which method should be used to initialise the cluster references. It should be 'random' or 'kmean++'. 'random' will sample k datum as the initial references, with uniform probabilities (by default).'kmean++' selects k datum using the kmean++ algorithm. It select a first center at random and iteratively select a new center with a probability proportional to the distance of each point to their nearest already choosen center, until k centers have been selected. |
  |       | \-\-seed arg        | A value to seed the random number generator. |
  |       | \-\-debug           | Enables debuggin verbosity. |


### Running an example

The file data.txt contains the number reads, from a H3K4me3 ChIP-seq experiment performed in CD4+ cells, that are mapped +/- 1kb around 23360 human TSSs within bins of 100bp (to reproduce this matrix, run [ChIP-Extract](https://ccg.vital-it.ch/chipseq/chip_extract.php) example by clicking the "Example" button and "Submit". Then, remove the first line and column of the resulting matrix which are headers). There are 99 bins per row (49 bins of 100bp upstream the TSS + the central bins containing the TSS + 49 bins upstream of the TSS). Here are the 4 first lines :

```
4 2 1 2 6 2 3 6 3 0 0 0 1 1 1 0 1 3 1 4 0 0 0 1 0 1 2 1 1 0 0 1 0 0 1 1 0 1 0 0 0 2 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 1 0 1 1 1 1 1 0 1 1 0 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 0 2 0 0 1 1 0 1 0 0 0 0 1 0 0 1
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 3
0 0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 1 1 1 0 0 0 0 1 1 2 2 6 0 1 0 0 1 1 3 8 4 1 0 0 1 0 0 1 0 0 0 0 0 0 3 4 1 18 25 13 4 2 3 0 1 2 7 12 3 5 4 2 2 9 10 8 8 9 1 1 1 0 5 5 7 3 2 0 3 0 1 0 0 0 0 0 0 0 0 0 0 0 0
1 0 0 1 3 5 2 1 1 0 0 0 0 1 1 1 0 0 0 1 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1 2 2 0 0 0 0 2 0 1 1 2 0 12 5 5 0 1 1 0 1 0 1 0 1 1 0 2 0 3 2 1 1 1 1 1 0 0 2 0 0 1 1 0 0 0 1 0 1 0 0 6 1 0 0 0 0 0 0

```
The data were organized such that all the TSS are oriented in the same direction. The first bin downstream each TSS is located in column 51. To partition these data into 3 clusters, based on the ChIP-seq profile in these regions, set a reasonable shifting freedom (7 bins, meaning +/-3*20bp) but no flipping (the TSSs are already oriented in the same direction), run :

```
bin/spark --data data.txt --cluster 4 --shift 7 --iter 30 --seeding kmean++ --seed 1234
```

If you are using the Singularity image, use :
```
singularity exec spar-k.simg spark --data data.txt --cluster 4 --shift 7 --iter 30 --seeding kmean++ --seed 1234
```
If you are using the Docker image, use :
```
docker run -i spar-k spark --data data.txt --cluster 4 --shift 7 --iter 30 --seeding kmean++ --seed 1234
```
The data 'data.txt' are contained within the image so there is no need to create a mount point between the host file system and the container file system yet. However, for cases where the data are outside the container (on the host file system) or when the results have to be sent outside the container (to the host file system), this will be required. It can be done as follows :

```
docker run -i -v <current dir>:/mount spar-k spark --data /mount/data_from_host.txt --cluster 4 --shift 7 --iter 30 --seeding kmean++ --seed 1234
```
where \<current dir\> is the absolute path to the current directory. On linux plateforms, you can use '$(pwd)'. Examples with a mount points can be found below.

As SPar-K implementation is fully multi-threaded, you can speed up the partitioning processes by dispatching the computations on several CPU cores. To do so, you need to use the -p option. For instance, to use 4 concurrent threads :
```
bin/spark --data data.txt --cluster 4 --shift 7 --iter 30 --seeding kmean++ --seed 1234 -p 4 > results.txt
```

With the Singularity image : 
```
singularity exec spar-k.simg spark --data data.txt --cluster 4 --shift 7 --iter 30 --seeding kmean++ --seed 1234 -p 4 > results.txt
```

With the Docker image :
```
docker run -i spar-k spark --data data.txt --cluster 4 --shift 7 --iter 30 --seeding kmean++ --seed 1234 -p 4 > results.txt
```

### spark_plot_heatmap.R
This program is an R script. Once a dataset has been partitioned using SPar-K, this script produces a heatmap of the results. 

Let's follow again the previous example. Now that you have your partition, you would like to display a nice heatmap. You would like to have the regions grouped by cluster and realigned as SPar-K aligned them. You can produce a plot of the data, realigned and ordered by cluster using :

```
Rscript bin/spark_plot_heatmap.R --data data.txt --partition results.txt --shift 7 --from -1000 --to 1000 --title "TSS with H3K4me3" --output myplot.png
```

With the Singularity image :
```
singularity exec spar-k.simg spark_plot_heatmap.R --data data.txt --partition results.txt --shift 7 --from -1000 --to 1000 --title "TSS with H3K4me3" --output myplot.png
```

With the Docker image :
```
docker run -i -v <current dir>:/mount spar-k spark_plot_heatmap.R --data data.txt --partition /mount/results.txt --shift 7 --from -1000 --to 1000 --title "TSS with H3K4me3" --output /mount/myplot.png
```
You noticed here the use of a mount point to read data from the host file system and to send the results to the host file system.

To get the help, run :

```
Rscript bin/spark_plot_heatmap.R --help
```

With the Singularity image :
```
singularity exec spar-k.simg spark_plot_heatmap.R --help
```

With the Docker image :
```
docker run -i spar-k spark_plot_heatmap.R --help
```

### spark_correct_sga.R
This program is an R script. Once a dataset has been partitioned using SPar-K, this scripts allows to update the corresponding SGA file according to the shift and flip values reported by SPar-K. 

Let's use the previous partitioning example (again). You have partitioned a dataset containing 23360 rows of length 99 with a shifting freedom of 7 and without flipping, the results are stored in results.txt and the TSS positions in a SGA file named references.sga ([about the SGA file format](https://ccg.vital-it.ch/chipseq/sga_specs.php)). Here are the first 4 lines :
```
NC_000001.10	TSS	   861123	+	1	SAMD11_1	2
NC_000001.10	TSS	   874653	+	1	SAMD11_2	2
NC_000001.10	TSS	   894631	-	1	NOC2L_1	2
NC_000001.10	TSS	   895964	+	1	KLHL17_1	2
```

Then, you can update the positions according to what SPar-K found to be the optimal alignment using :

```
Rscript bin/spark_correct_sga.R --sga references.sga --partition results.txt --shift 7 --ncol 99 --binSize 20
```

With the Singularity image (note that references.sga is inside the image) :
```
singularity exec spar-k.simg spark_correct_sga.R --sga references.sga --partition results.txt --shift 7 --ncol 99 --binSize 20
```
With the Docker image (note that references.sga is inside the image) :
```
docker run -i -v <current dir>:/mount spar-k spark_correct_sga.R --sga references.sga --partition /mount/results.txt --shift 7 --ncol 99 --binSize 20
```

If you want to correct only the reference positions of regions which were assigned to a given cluster - let's say cluster 2 - then you can run :

```
Rscript bin/spark_correct_sga.R --sga references.sga --partition results.txt --shift 7 --ncol 99 --binSize 20 --cluster 2
```

With the Singularity image (note that references.sga is inside the image):
```
singularity exec spar-k.simg spark_correct_sga.R --sga references.sga --partition results.txt --shift 7 --ncol 99 --binSize 20 --cluster 2
```

With the Docker image (note that references.sga is inside the image) :
```
docker run -i -v <current dir>:/mount spar-k spark_correct_sga.R --sga references.sga --partition /mount/results.txt --shift 7 --ncol 99 --binSize 20 --cluster 2
```

For help, run :

```
Rscript bin/spark_correct_sga.R --help
```

With the Singularity image :
```
singularity exec spar-k.simg spark_correct_sga.R --help
```

With the Docker image :
```
docker run -i spar-k spark_correct_sga.R --help
```

## Authors

* **Romain Groux**


## License

This project is licensed under the GNU General Public License v3 - see the [LICENSE.md](LICENSE.md) file for details


## Acknowledgments

* Philipp Bucher
* René Dreos
* Giovanna Ambosini


