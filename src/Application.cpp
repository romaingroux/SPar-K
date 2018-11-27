#include "Application.hpp"

#include <iostream>
#include <vector>
#include <boost/program_options.hpp>
#include <Matrix/Matrix2D.hpp>
#include <stdexcept>                 // std::invalid_argument

#include <Clustering/ClusteringEngine.hpp>
#include <Clustering/KmeanEngine.hpp>
#include <Random/RandomNumberGenerator.hpp>
#include <Utility/Constants.hpp>      // Constants
#include <Utility/Matrix_utility.hpp> // smooth_outliers
#include <GUI/ConsoleProgressBar.hpp> // ProgressBar



// some global variables
std::string version("SPar-K v1.0") ;
// possible seeding mode options
static std::string seeding_random   = "random" ;
static std::string seeding_kmean_pp = "kmean++" ;
// possible distances
static std::string dist_cor         = "corr" ;
static std::string dist_normcor     = "normcorr" ;



namespace po = boost::program_options ;


Application::Application(int argn, char **argv)
    : exit_code(0)
{   // try to set the options from what is given on the command line
    this->setOptions(argn, argv) ;
}


Application::~Application()
{}


void Application::setOptions(int argn, char** argv)
{
    // initialize the options to some default values
    this->options.file_data = "" ;
    this->options.file_references = "" ;
    this->options.thread_n = 1 ;
    this->options.iteration_n = 1 ;
    this->options.cluster_n = 1 ;
    this->options.shift_n = 1 ;
    this->options.flip = false ;
    this->options.no_outlier = false ;
    this->options.debug = false ;
    this->options.seed = "" ;
    this->options.seeding = seeding_random.c_str() ;
    this->options.distance = dist_cor ;

    if(argv == nullptr)
    {   this->exit_code = -1 ;
        return ;
    }

    // initialize variables to parse options
    char desc_msg[4096] ;
    sprintf(desc_msg, "\n"
                          "SPar-K is a program to cluster a matrix of signal vectors (of any nature).\n"
                          "Signal will be classified into K clusters using a K-mean like method except that shifting\n"
                          "and flipping the vectors can be enabled to realign the data (see the corresponding\n"
                          "options below). More formally, SPar-K will partition the given dataset into K alignments\n"
                          "which are as many clusters. Each cluster is reprented by a signal density vector representing\n"
                          "the mean signal at each position along the corresponding alignment.\n\n"
                          "%s\nCompiled on %s\n\n", version.c_str(), __DATE__) ;

    boost::program_options::variables_map vm ;
    boost::program_options::options_description desc(desc_msg) ;

    std::string opt_help_msg       = "Produces this help message";
    std::string opt_version_msg    = "Prints the version number";
    std::string opt_parallel_msg   = "The number of threads dedicated to the "
                                     "computations, by default 1.";
    std::string opt_data_msg       = "The data file address.";
    std::string opt_references_msg = "The cluster reference pattern file address.";
    std::string opt_iter_msg       = "The maximum number of iterations." ;
    std::string opt_cluster_msg    = "The number of cluster to find." ;
    std::string opt_shift_msg      = "Enables this number of column of shifting "
                                     "freedom. By default, shifting is "
                                     "disabled (equivalent to --shift 1)." ;
    std::string opt_flip_msg       = "Enables flipping.";
    std::string opt_no_outlier_msg = "Pre-pcocess the data to smooth out outliers from the data in a "
                                     "row-wise manner. Each row is searched for outliers which are defined "
                                     "as any value bigger/smaller than the row mean +/- 3*row standard deviation. "
                                     "If a value is an outlier it is replaced by the mean of its left and right "
                                     "neighbours. If a has only a left/right neighbour, then the two left/right "
                                     "neighbours are averaged. If several outliers follow each other the above "
                                     "process is applied to the values in a left to right order and at the end "
                                     "the new averaged values may still be outliers." ;

    char seeding_msg[2048] ;
    sprintf(seeding_msg,
            "Specify which method should be used to initialise the "
            "cluster references. It should be '%s' or '%s'. "
            "'%s' will sample k datum as the initial references, with uniform probabilities (by default)."
            "'%s' selects k datum using the kmean++ algorithm. It select a first center at random and "
            "iteratively select a new center with a probability proportional to the distance of "
            "each point to their nearest already choosen center, until k centers have been selected.",
            seeding_random.c_str(),
            seeding_kmean_pp.c_str(),
            seeding_random.c_str(),
            seeding_kmean_pp.c_str()) ;
    std::string opt_seeding_msg = seeding_msg ;
    std::string opt_seed_msg       = "A value to seed the random number generator.";

    char dist_msg[2048] ;
    sprintf(dist_msg,
            "Specify which distance should be used during the clustering. It should be '%s' (by default) or '%s'. ",
            dist_cor.c_str(), dist_normcor.c_str()) ;
    std::string opt_dist_msg = dist_msg ;

    std::string opt_debug_msg      = "Enables debuggin verbosity.";

    desc.add_options()
            ("help,h",       opt_help_msg.c_str())
            ("version,v",    opt_version_msg.c_str())

            ("parallel,p",   po::value<std::size_t>(&(this->options.thread_n)),  opt_parallel_msg.c_str())
            ("data,d",       po::value<std::string>(&(this->options.file_data)), opt_data_msg.c_str())
            ("references,r", po::value<std::string>(&(this->options.file_references)), opt_references_msg.c_str())

            ("iter,i",       po::value<size_t>(&(this->options.iteration_n)),    opt_iter_msg.c_str())
            ("cluster,c",    po::value<size_t>(&(this->options.cluster_n)),      opt_cluster_msg.c_str())
            ("shift,s",      po::value<size_t>(&(this->options.shift_n)),        opt_shift_msg.c_str())
            ("flip",         opt_flip_msg.c_str())
            ("nooutlier",    opt_no_outlier_msg.c_str())
            ("dist",         po::value<std::string>(&(this->options.distance)),  opt_dist_msg.c_str())

            ("seeding",      po::value<std::string>(&(this->options.seeding)),   opt_seeding_msg.c_str())
            ("seed",         po::value<std::string>(&(this->options.seed)),      opt_seed_msg.c_str())

            ("debug",        opt_debug_msg.c_str()) ;

    // parse
    try
    {   po::store(po::parse_command_line(argn, argv, desc), vm) ;
        po::notify(vm) ;
    }
    catch(std::invalid_argument& e)
    {   std::cerr << "Error! Invalid option given!" << std::endl
                  << e.what() << std::endl ;
        this->exit_code = EXIT_FAILURE ;
    }
    catch(...)
    {	std::cerr << "An unknown error occured while parsing the options" << std::endl ;
        this->exit_code = EXIT_FAILURE ;
    }

    // checks unproper option settings
    if(this->options.file_data == "" and (not vm.count("help")) and (not vm.count("version")))
    {   std::cerr << "Error! No data file was given (--data)!" << std::endl ;
        this->exit_code = EXIT_FAILURE ;
    }
    else if((this->options.seeding != seeding_random) and
            (this->options.seeding != seeding_kmean_pp))
    {   std::cerr << "Error! Unrecognized seeding method (--seeding)!" << std::endl ;
        this->exit_code = EXIT_FAILURE ;
    }
    else if((this->options.distance != dist_cor) and
            (this->options.distance != dist_normcor))
    {   std::cerr << "Error! Unrecognized distance (--dist)!" << std::endl ;
        this->exit_code = EXIT_FAILURE ;
    }

    if(vm.count("help"))      { std::cout << desc << std::endl ; this->exit_code = EXIT_FAILURE ; }
    if(vm.count("version"))   { std::cout << version << std::endl ; this->exit_code = EXIT_FAILURE ; }
    if(vm.count("flip"))      { this->options.flip  = true ; }
    if(vm.count("nooutlier")) { this->options.no_outlier  = true ; }
    if(vm.count("debug"))     { this->options.debug = true ; }
}


bool Application::isDebugOn()
{   return this->options.debug ; }


int Application::run()
{   try
    {
        if(this->isDebugOn())
        {   std::cerr << "data file      " << this->options.file_data   << std::endl ;
            std::cerr << "reference file " << this->options.file_data   << std::endl ;
            std::cerr << "iteration      " << this->options.iteration_n << std::endl ;
            std::cerr << "clusters       " << this->options.cluster_n   << std::endl ;
            std::cerr << "shifts         " << this->options.shift_n     << std::endl ;
            std::cerr << "flip           " << this->options.flip        << std::endl ;
            std::cerr << "nooutler       " << this->options.no_outlier  << std::endl ;
            std::cerr << "distance       " << this->options.distance    << std::endl ;
            std::cerr << "seed           " << this->options.seed        << std::endl ;
            std::cerr << "seeding        " << this->options.seeding     << std::endl ;
            std::cerr << "threads        " << this->options.thread_n    << std::endl ;
        }
        // something occured, cannot continue
        if(this->exit_code != 0)
        {   return this->exit_code ; }

        // run clustering
        else
        {   // read data
            Matrix2D<double> data(this->options.file_data) ;
            if(this->options.no_outlier)
            {   smooth_outliers(data) ; }

            // instantiate the clustering object
            ClusteringEngine* engine = nullptr ;
            Constants::dist_type dist_type ;
            if(this->options.distance == dist_cor)
            {   dist_type = Constants::dist_type::COR_DIST ; }
            if(this->options.distance == dist_normcor)
            {   dist_type = Constants::dist_type::NORM_COR_DIST ; }
            // instantiate the clustering object
            engine = new KmeanEngine(data,
                                     this->options.cluster_n,
                                     this->options.shift_n,
                                     this->options.flip,
                                     this->options.seeding,
                                     this->options.seed,
                                     this->options.thread_n,
                                     dist_type) ;

            size_t i = 0 ;
            ConsoleProgressBar bar(std::cerr, this->options.iteration_n,  20, "clustering") ;
            Constants::clustering_codes code = Constants::clustering_codes::SUCCESS ;
            while((i < this->options.iteration_n) and
                  (code != Constants::clustering_codes::FAILURE) and
                  (code != Constants::clustering_codes::CONVERGENCE))
            {   code = engine->cluster() ;
                bar.update() ; bar.display();
                i++ ;
            }
            switch(code)
            {   case Constants::clustering_codes::CONVERGENCE:
                {   this->exit_code = EXIT_SUCCESS ;
                    for(size_t j=i; j<=this->options.iteration_n; j++)
                    {   bar.update() ; }
                    bar.display() ;
                    std::cerr << std::endl ;
                    std::cerr << "Clustering converged after "
                               << i << " iterations"
                               << std::endl ;
                    break ;
                }
                case Constants::clustering_codes::SUCCESS:
                {   std::cerr << std::endl ;
                    std::cerr << "Clustering normally terminated after "
                              << i << " iterations"
                              << std::endl ;
                    break ;
                }
                case Constants::clustering_codes::FAILURE:
                {   this->exit_code = EXIT_FAILURE ;
                    for(size_t j=i; j<=this->options.iteration_n; j++)
                    {   bar.update() ; }
                    bar.display() ;
                    std::cerr << std::endl ;
                    std::cerr << "Clustering failed after "
                               << i << " iterations"
                               << std::endl ;
                    break ;
                }
                default:
                {   std::cerr << std::endl ;
                    std::cerr << "If you see this message, this is a bug!" << std::endl
                              << "ClusteringEngine::cluster returned an unexpected value."
                              << std::endl ;
                    break ;
                }
            } ;
            // if the clustering fails, print the data as they were at the last iteration
            // anyway, may be informative
            engine->print_results(std::cout) ;
            delete engine ;
        }
    }
    catch(std::logic_error& e)
    {   std::cerr << "Something occurred : " << std::endl
                  << e.what() << std::endl ;
        this->exit_code = EXIT_FAILURE ;
    }
    return this->exit_code ;
}



using namespace std ;



int main(int argn, char** argv)
{
    Application application(argn, argv) ;
    int exit_code = application.run() ;
    return exit_code ;
}

