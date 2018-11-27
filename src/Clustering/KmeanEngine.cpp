#include "KmeanEngine.hpp"

#include <cmath>      // ceil()
#include <functional> // bind(), ref()
#include <utility>    // move()
#include <algorithm>  // min_element()
#include <iomanip>
#include <limits>     // numeric_limits<size_t>::max(), numeric_limits<double>::max()
#include <Clustering/CorDistanceComputer.hpp>
#include <Clustering/NormCorDistanceComputer.hpp>
#include <Statistics/Distances.hpp>
#include <Statistics/Statistics.hpp>
#include <Utility/Utility.hpp>
#include <Utility/Constants.hpp>     // Constants::corr_fail, Constants::dist_type
#include <Random/Random.hpp>
#include <Utility/Utility.hpp>       // isNan() and flip enum
#include <Utility/Sorting.hpp>       // order()


KmeanEngine::~KmeanEngine()
{   // join the threads
    this->joinThreads() ;
    // free memory
    for(auto& ptr : this->dist_computers)
    {   if(ptr != nullptr)
        {   delete ptr ;
            ptr = nullptr ;
        }
    }
}


KmeanEngine::KmeanEngine(const Matrix2D<double>& data,
                         size_t cluster_n, size_t shift_n, bool flip,
                         const std::string& seeding_method,
                         const std::string& seed,
                         size_t n_threads,
                         Constants::dist_type dist_type)
    throw (std::invalid_argument)
    : distances(), distances_previous(), cluster_sizes(), threads(n_threads),
      dist_computers(n_threads, nullptr)
{   // this can throw a std::invalid_argument
    this->constructorInit(data, cluster_n, shift_n, flip, seeding_method, seed, dist_type) ;
}

KmeanEngine::KmeanEngine(const std::string& path_data,
                         size_t cluster_n, size_t shift_n, bool flip,
                         const std::string& seeding_method,
                         const std::string& seed,
                         size_t n_threads,
                         Constants::dist_type dist_type)
    throw (std::invalid_argument)
    // this can throw a std::invalid_argument
    : KmeanEngine::KmeanEngine(Matrix2D<double>(path_data),
                               cluster_n, shift_n, flip,
                               seeding_method, seed,
                               n_threads,
                               dist_type)
{}


KmeanEngine::KmeanEngine(const Matrix2D<double>& data,
                         const Matrix2D<double>& references,
                         size_t shift_n, bool flip,
                         size_t n_threads,
                         Constants::dist_type dist_type)
    throw (std::invalid_argument)
    : distances(), distances_previous(), cluster_sizes(), threads(n_threads),
      dist_computers(n_threads, nullptr)
{   // this can throw a std::invalid_argument
    this->constructorInit(data, references, shift_n, flip, dist_type) ;
}


KmeanEngine::KmeanEngine(const std::string& path_data,
                         const std::string& path_references,
                         size_t shift_n, bool flip,
                         size_t n_threads,
                         Constants::dist_type dist_type)
    throw (std::invalid_argument)
    // this can throw a std::invalid_argument
    : KmeanEngine::KmeanEngine(Matrix2D<double>(path_data),
                               Matrix2D<double>(path_references),
                               shift_n, flip,
                               n_threads,
                               dist_type)
{}


void KmeanEngine::constructorInit(const Matrix2D<double>& data,
                                  size_t cluster_n, size_t shift_n, bool flip,
                                  const std::string& seeding_method,
                                  const std::string& seed,
                                  Constants::dist_type dist_type)
    throw (std::invalid_argument)
{
    size_t n_row = data.get_nrow() ;
    size_t n_col = data.get_ncol() ;

    // check that the number of cluster is smaller than the number of data points
    if(n_row <= cluster_n)
    {   this->joinThreads() ;
        char msg[2048] ;
        sprintf(msg, "Error! number of cluster is too high (%zd), there are " \
                     "only %zd data points!", cluster_n, n_row) ;
        throw std::invalid_argument(msg) ;
    }
    // check that the shifting freedom is smaller than the number of columns - 1
    // such that we always have slices of data of at least length 2
    if(n_col-1 <= shift_n)
    {   this->joinThreads() ;
        char msg[2048] ;
        sprintf(msg, "Error! shifting freedom is too high (%zd), there are " \
                     "only %zd columns in the data!", shift_n, n_col) ;
        throw std::invalid_argument(msg) ;
    }

    // common to all constructors
    this->constructorRoutine(data, cluster_n, shift_n, flip, dist_type) ;
    // specific to this constructor
    this->references         = Matrix2D<double>(cluster_n, n_col, 0.) ;
    this->distances_previous = std::vector<double>(n_row, 0.) ;
    this->clusters_previous  = std::vector<size_t>(n_row, 0) ;
    this->shifts_previous    = std::vector<std::vector<size_t>>(n_row, std::vector<size_t>(2,0)) ;
    this->flips_previous     = std::vector<int>(n_row, 0) ;

    // seeds the random number generator BEFORE USING IT
    if(seed != std::string(""))
    {   getRandomGenerator(seed) ; }
    else
    {   getRandomGenerator() ; }
    // seeds the references
    this->seeding(seeding_method) ;
}


void KmeanEngine::constructorInit(const Matrix2D<double>& data,
                                  const Matrix2D<double>& references,
                                  size_t shift_n, bool flip,
                                  Constants::dist_type dist_type)
throw (std::invalid_argument)
{
    size_t n_row     = data.get_nrow() ;
    size_t n_col     = data.get_ncol() ;
    size_t n_cluster = references.get_nrow() ;

    // check that the number of cluster is smaller than the number of data points
    if(n_row <= n_cluster)
    {   this->joinThreads() ;
        char msg[2048] ;
        sprintf(msg, "Error! number of cluster is too high (%zd), there are " \
                     "only %zd data points!", n_cluster, n_row) ;
        throw std::invalid_argument(msg) ;
    }
    // check that the shifting freedom is smaller than the number of columns - 1
    // such that we always have slices of data of at least length 2
    if(n_col-1 <= shift_n)
    {   this->joinThreads() ;
        char msg[2048] ;
        sprintf(msg, "Error! shifting freedom is too high (%zd), there are " \
                     "only %zd columns in the data!", shift_n, n_col) ;
        throw std::invalid_argument(msg) ;
    }
    // common to all constructors
    this->constructorRoutine(data, n_cluster, shift_n, flip, dist_type) ;
    // specific to this constructor
    this->references = Matrix2D<double>(references) ;
    // not needed but initialised anyway
    this->distances_previous = std::vector<double>() ;
    this->clusters_previous  = std::vector<size_t>(0) ;
    this->shifts_previous    = std::vector<std::vector<size_t>>(n_row, std::vector<size_t>(2,0)) ;
    this->flips_previous     = std::vector<int>(0) ;
}


std::vector<double> KmeanEngine::getDistances() const
{   return this->distances ; }


std::vector<size_t> KmeanEngine::getClusters() const
{   return this->clusters ; }


std::vector<std::vector<size_t>> KmeanEngine::getShifts() const
{   return this->shifts ; }


std::vector<int> KmeanEngine::getFlips() const
{   return this->flips ; }


Constants::clustering_codes KmeanEngine::cluster()
{
    // keep track of the previous results
    if(this->n_iter > 0)
    {   this->clusters_previous  = this->clusters ;
        this->shifts_previous    = this->shifts ;
        this->flips_previous     = this->flips ;
        this->distances_previous = std::vector<double>(this->distances) ;
    }

    // slice data for the threads, each one will run on [from,to)
    size_t n_thread = this->getThreads().getNThread() ;
    size_t nrow     = this->data.get_nrow() ;
    size_t by       = nrow / n_thread ;
    size_t from = 0, to = 0 ;
    std::vector<std::future<bool>>  v_ftr_done(n_thread) ;
    std::vector<std::promise<bool>> v_prms_done(n_thread) ;
    // assign the datum to the clusters
    // ----------------------------- threads starts -----------------------------
    // only threads access data, references, shifts, flips, clusters, distances
    for(size_t n=0; n<n_thread; n++)
    {   from = (to == 0) ? 0 : (to) ;
        to   = from + by ;
        (to + by > nrow) ? (to = nrow) : 0 ;
        v_ftr_done[n] = v_prms_done[n].get_future() ;
        // assign datum to clusters
        this->getThreads().addJob(std::move(std::bind(&KmeanEngine::assign_to_group, this, from, to, this->dist_computers[n], std::ref(v_prms_done[n])))) ;
    }
    // the functions in the threads should set promises to instruct the program they are done
    // working. If a thread is not done, the program waits here until all threads are done.
    for(auto& f : v_ftr_done)
    {   f.get() ; }
    // ----------------------------- threads stops -----------------------------

    // update the clusters
    this->updateClusterSizes() ;
    size_t n_cluster = this->cluster_sizes.size() ;
    std::vector<std::future<std::vector<double>>>  v_ftr_ref(n_cluster) ;
    std::vector<std::promise<std::vector<double>>> v_prms_ref(n_cluster) ;
    // ----------------------------- threads starts -----------------------------
    // only threads access data, references, shifts, flips, clusters, distances
    for(size_t i=0; i<n_cluster; i++)
    {   v_ftr_ref[i] = v_prms_ref[i].get_future() ;
        this->getThreads().addJob(std::move(std::bind(&KmeanEngine::computeClusterReference, this, i, std::ref(v_prms_ref[i])))) ;
    }
    // this function access the futures corresponding to the thread promises. If a thread is not done yet
    // the function will wait until it is done.
    this->updateReferences(v_ftr_ref) ;
    // ----------------------------- threads stops -----------------------------

    this->n_iter++ ;

    // checks whether a cluster is empty (this is a failure) and returns
    bool success     = not this->hasEmptyCluster() ;
    bool convergence = this->hasConverged() ;
    if(convergence)
    {   return Constants::clustering_codes::CONVERGENCE ; }
    else if(success)
    {   return Constants::clustering_codes::SUCCESS ; }
    else
    {   return Constants::clustering_codes::FAILURE ; }
}


void KmeanEngine::seeding(const std::string& method) throw(std::runtime_error)
{   if(method == "random")
    {   this->seeding_random() ; }
    else if(method == "kmean++")
    {   this->seeding_kmean_pp() ; }
    else
    {   char msg[1024] ;
        sprintf(msg, "Error! unrecognized seeding method given : %s", method.c_str()) ;
        throw std::runtime_error(msg) ;
    }
}


void KmeanEngine::print_results(std::ostream& o)
{
    Matrix2D<double> results ;
    size_t nrow = this->data.get_nrow() ;
    if(not this->flip)
    {   results    = Matrix2D<double>(nrow, 4, 0) ;
        for(size_t i=0; i<nrow; i++)
        {   results(i,0) = this->clusters[i] + 1 ;
            results(i,1) = this->shifts[i][0] ;
            results(i,2) = this->shifts[i][1] ;
            results(i,3) = this->distances[i] ;
        }
        o << "cluster   shift_ref   shift_dat   distance" << std::endl ;
    }
    else
    {   results    = Matrix2D<double>(nrow, 5, 0) ;
        for(size_t i=0; i<nrow; i++)
        {   results(i,0) = this->clusters[i] + 1;
            results(i,1) = this->shifts[i][0] ;
            results(i,2) = this->shifts[i][1] ;
            results(i,3) = this->flips[i] ;
            results(i,4) = this->distances[i] ;
        }
        o << "cluster   shift_ref   shift_dat   flip   distance" << std::endl ;
    }
    o << results << std::endl ;
}


void KmeanEngine::joinThreads()
{   this->getThreads().join() ; }


ThreadPool& KmeanEngine::getThreads()
{   return this->threads ; }


void KmeanEngine::seeding_random()
{   size_t n_row     = this->data.get_nrow() ;
    // size_t n_col     = this->data.get_ncol() ;
    size_t n_cluster = this->references.get_nrow() ;

    std::vector<bool> choosen(n_row, false) ;

    for(size_t i=0; i<n_cluster; )
    {
        size_t index =  rand_int_uniform(size_t(0), n_row-1) ;

        // already choosen as reference
        if(choosen[index])
        { ; }
        // not yet choosen as reference
        else
        {   this->references.set_row(i, this->data.get_row(index)) ;
            choosen[index] = true ;
            i++ ;
        }
    }
}


void KmeanEngine::seeding_kmean_pp()
{   size_t n_row     = this->data.get_nrow() ;
    // size_t n_col     = this->data.get_ncol() ;
    size_t n_cluster = this->references.get_nrow() ;

    // Matrix2D<double> dist_tmp(flip+1, this->n_shift, 0.) ;  // to store tmp distances
    std::vector<bool> choosen(n_row, false) ;         // is a datum a center
    std::vector<double> dist_to_closest(n_row, 3.) ;  // datum-closest center distance
    size_t index ;                                    // the index of the sampled center

    // sample first center, avoid sampling a first center with sd == 0
    // otherwise all points will have a distance of 1 to this (if the
    // corr cannot be computed, it is set to 0).
    do
    {    index = rand_int_uniform(static_cast<size_t>(0), n_row-1) ; }
    while(sd(this->data.get_row(index)) == 0.) ;
    this->references.set_row(0, this->data.get_row(index)) ;
    choosen[index] = true ;

    // sample other centers
    size_t shift_ref, shift_dat ;
    bool flip_dat ;
    for(size_t k=1; k<n_cluster; k++)
    {   // find the distance to the nearest center for all points
        for(size_t i=0; i<n_row; i++)
        {   std::vector<double> dist_to_centers(k, 3.0) ;
            // all centers
            for(size_t j=0; j<k; j++)
            {   double dist_tmp ;
                // shift only
                if(not this->flip)
                {   dist_tmp = this->dist_computers[0]->compute_distance(this->references.get_row(j),
                                                                         this->data.get_row(i),
                                                                         this->n_shift,
                                                                         shift_ref,
                                                                         shift_dat) ;
                }
                // shift and flip
                else
                {   dist_tmp = this->dist_computers[0]->compute_distance(this->references.get_row(j),
                                                                         this->data.get_row(i),
                                                                         this->n_shift,
                                                                         shift_ref,
                                                                         shift_dat,
                                                                         flip_dat) ;
                }
                // dist data_i - reference_j
                dist_to_centers[j] = dist_tmp ;
            }
            // dist data_i - closest center
            dist_to_closest[i] = *min_element(dist_to_centers.begin(), dist_to_centers.end()) ;
        }

        // turn distances into weights
        std::vector<size_t> weights(n_row) ;
        for(size_t i=0; i<weights.size(); i++)
        {   weights[i] = static_cast<size_t>(dist_to_closest[i]*1000000000.) ; }
        // sample a new center using the weights (resample until finding a free datum if necessary)
        std::discrete_distribution<size_t> d(weights.begin(), weights.end()) ;
        for(index = d(getRandomGenerator()) ;
            choosen[index] ;
            index = d(getRandomGenerator()))
        {}
        this->references.set_row(k, this->data.get_row(index)) ;
        choosen[index] = true ;
    }
}


void KmeanEngine::constructorRoutine(const Matrix2D<double>& data, size_t cluster_n, size_t shift_n, bool flip, Constants::dist_type dist_type)
{   size_t nrow      = data.get_nrow() ;
    this->data       = Matrix2D<double>(data) ;
    this->clusters   = std::vector<size_t>(nrow, 0) ;
    this->shifts     = std::vector<std::vector<size_t>>(nrow, std::vector<size_t>(2,0)) ;
    this->distances  = std::vector<double>(nrow, 0.) ;
    // this will only be modified if clustering with flipping is performed
    this->flip       = flip ;
    this->flips      = std::vector<int>(nrow, Constants::FORWARD) ;
    // set shift states
    this->n_shift = shift_n ;
    // to compute distances
    this->allocateDistComputers(dist_type) ;
    // others
    this->n_iter = 0 ;
    this->cluster_sizes  = std::vector<double>(cluster_n, 0.) ;
}

void KmeanEngine::allocateDistComputers(Constants::dist_type type) throw (std::invalid_argument)
{   // correlation distance
    if(type == Constants::dist_type::COR_DIST)
    {   for(auto& ptr : this->dist_computers)
        {   ptr = new CorDistanceComputer(this->n_shift) ; }
    }
    // normalized correlation distance
    else if(type == Constants::dist_type::NORM_COR_DIST)
    {   for(auto& ptr : this->dist_computers)
        {   ptr = new NormCorDistanceComputer(this->n_shift) ; }
    }
    // unrecognized
    else
    {   throw std::invalid_argument("unknown distance type!") ; }
}

void KmeanEngine::resetClusterSizes()
{   for(size_t i=0; i<this->cluster_sizes.size(); i++)
    {   this->cluster_sizes[i] = 0 ;   }
}


bool KmeanEngine::hasEmptyCluster() const
{   for(auto i : this->cluster_sizes)
    {   if(i == 0.)
        {   return true ; }
    }
    return false ;
}


void KmeanEngine::assign_to_group(size_t from, size_t to, DistanceComputer* cd, std::promise<bool>& prms_done)
{
    // some dimensions
    size_t n_cluster = this->cluster_sizes.size() ;

    for(size_t i=from; i<to; i++)
    {   double best_dist      = std::numeric_limits<double>::max() ;  // next distance can only be smaller
        size_t best_cluster   = std::numeric_limits<size_t>::max() ;  // visible if in results as such
        size_t best_shift_dat = std::numeric_limits<size_t>::max() ;  // visible if in results as such
        size_t best_shift_ref = std::numeric_limits<size_t>::max() ;  // visible if in results as such
        size_t best_flip_dat  = std::numeric_limits<size_t>::max() ;  // visible if in results as such
        for(size_t j=0; j<n_cluster; j++)
        {   // with flipping
            if(this->flip)
            {   size_t shift_ref, shift_dat ;
                bool flip_dat ;
                double dist_tmp ;
                dist_tmp = cd->compute_distance(this->references.get_row(j), this->data.get_row(i), this->n_shift,
                                                shift_ref, shift_dat, flip_dat) ;
                if(dist_tmp < best_dist)
                {   best_dist      = dist_tmp ;
                    best_shift_ref = shift_ref ;
                    best_shift_dat = shift_dat ;
                    best_flip_dat  = flip_dat ;
                    best_cluster   = j ;
                }
            }
            // without flipping
            else
            {   size_t shift_ref, shift_dat ;
                double dist_tmp ;
                dist_tmp = cd->compute_distance(this->references.get_row(j), this->data.get_row(i), this->n_shift,
                                                shift_ref, shift_dat) ;
                if(dist_tmp < best_dist)
                {   best_dist      = dist_tmp ;
                    best_shift_ref = shift_ref ;
                    best_shift_dat = shift_dat ;
                    best_cluster   = j ;
                }
            }
        }

        // update the parameters
        // cannot update the cluster size now, this function is run in different thread -> overwritting danger
        this->clusters[i]  = best_cluster ;
        this->shifts[i][0] = best_shift_ref ;
        this->shifts[i][1] = best_shift_dat ;
        this->distances[i] = best_dist ;
        if(this->flip)
        {   this->flips[i] = best_flip_dat ; }
    }
    prms_done.set_value(true) ;
}


void KmeanEngine::updateReferences(std::vector<std::future<std::vector<double>>>& v_ftr_ref)
{   size_t n_cluster = this->cluster_sizes.size() ;
    // if threads are not finished working, program waits here until
    // that all are done
    for(size_t i=0; i<n_cluster; i++)
    {   std::vector<double> v_ref = v_ftr_ref[i].get() ;
        for(size_t j=0; j<v_ref.size(); j++)
        {   this->references(i,j) = v_ref[j] ; }
    }
}


void KmeanEngine::computeClusterReference(size_t cluster_id, std::promise<std::vector<double>>& prms_ref)
throw (std::invalid_argument)
{
    // check that cluster label is valid
    if(not (cluster_id < this->cluster_sizes.size()))
    {   char msg[2048] ;
        sprintf(msg, "Error! Invalid cluster label : %zd", cluster_id) ;
        throw std::invalid_argument(msg) ;
    }

    // in this function, every size_t is turned to int because of iterator arithmetic operations
    // the results are not expected to be negative however, it leads to expression containing signed and
    // unsigned values -> turn everything to int to avoid it

    int n_row = this->data.get_nrow() ;
    int n_col    = this->data.get_ncol() ;
    int l_slice  = n_col - this->n_shift + 1 ;
    // will contain the reference for this cluster
    std::vector<double> reference(n_col, 0) ;

    // nothing to do and avoids division by 0
    if(this->cluster_sizes[cluster_id] == 0)
    {   // needs to fill the promise anyway for the main thread to see the end of this
        // method and to give it the reference for this cluster (a flat reference)
        prms_ref.set_value(reference) ;
        return ;
    }

    // the number of values which will be stacked at each position along the aggregated reference
    std::vector<int> n_aligned_here(n_col, 0) ;

    for(int i=0; i<n_row; i++)
    {
        if(this->clusters[i] != cluster_id)
        {   continue ; }
        else
        {   // forward orientation
            if(this->flips[i] == Constants::FORWARD)
            {
                int from_dat = this->shifts[i][1] ;
                int to_dat   = from_dat + l_slice ; // slice is [from_dat,to_dat)
                int from_ref = this->shifts[i][0] ;
                // int to_ref   = from_ref + l_slice ;
                for(int j_dat=from_dat, j_ref=from_ref; j_dat<to_dat; j_dat++, j_ref++)
                {   reference[j_ref] += this->data(i, j_dat) ;
                    n_aligned_here[j_ref] += 1 ;
                }
            }
            // reverse
            else
            {   int from_dat = this->shifts[i][1] ;
                int to_dat   = from_dat + l_slice ; // slice is [from_dat->to_dat) but it should be reversed
                int from_dat_rev = to_dat - 1 ;
                int to_dat_rev   = from_dat  - 1 ;   // slice is (to_dat_rev<-from_dat]
                int from_ref = this->shifts[i][0] ;
                for(int j_dat=from_dat_rev, j_ref=from_ref; j_dat>to_dat_rev; j_dat--, j_ref++)
                {   reference[j_ref] += this->data(i, j_dat) ;
                    n_aligned_here[j_ref] += 1 ;
                }
            }
        }
    }

    // take the mean
    for(size_t i=0; i<reference.size(); i++)
    {   // avoid div by 0
        if(n_aligned_here[i] != 0)
        {   reference[i] /= static_cast<double>(n_aligned_here[i]) ; }
    }

    prms_ref.set_value(reference) ;
}


void KmeanEngine::print_alignment(std::ostream& o) const
{
    // in this function, every size_t is turned to int because of iterator arithmetic operations
    // the results are not expected to be negative however, it leads to expression containing signed and
    // unsigned values -> turn everything to int to avoid it

    int n_row        = this->data.get_nrow() ;
    int n_col        = this->data.get_ncol() ;
    size_t n_cluster = this->references.get_nrow() ;

    // a constant for forward to reverse iterator conversion
    int c = ((n_col - this->n_shift + 1) % 2 == 0 ? -1 : 0 ) ;

    for(size_t k=0; k<n_cluster; k++)
    {
        // print the reference of the cluster
        o << "REF     " << std::setw(3) << k << " " ;
        for(auto e : this->references.get_row(k))
        {   o << std::setprecision(2) << std::setw(3) << e << ' ' ; }
        o << std::endl ;


        for(int i=0; i<n_row; i++)
        {   if(this->clusters[i] == k)
            {   o << "DAT " << (this->flips[i] == Constants::FORWARD ? "FW  " : "REV ") << std::setw(3) << i << " " ;
                int from_dat = this->shifts[i][1] ;
                int from_ref = this->shifts[i][0] ;
                int to_ref   = n_col - this->n_shift + from_ref ;

                for(int j_ref=0; j_ref<n_col; j_ref++)
                {   // remap reference iterator to data
                    int j_dat_fw = j_ref + (this->shifts[i][1] - this->shifts[i][0]) ;
                    // data was mapped in reverse
                    if(this->flips[i])
                    {   // if in reverse, the data iterator should be sent at the other extremity of the mapped region
                        int j_dat_rev = j_dat_fw + ((n_col-this->n_shift+1) -2*(j_dat_fw-from_dat)) + c ;
                        // we're out of the aligned region
                        if(j_ref < from_ref or j_ref > to_ref)
                        {   o << std::setprecision(2) << std::setw(3) << '*' << ' ' ; }
                        // inside the aligned region
                        else
                        {   o << std::setprecision(2) << std::setw(3) << (this->data(i,j_dat_rev) < 1e-10 ? 0 : this->data(i,j_dat_rev)) << ' ' ; }
                    }
                    // data was mapped in forward
                    else
                    {   // we're out of the aligned region
                        if(j_ref < from_ref or j_ref > to_ref)
                        {   o << std::setprecision(2) << std::setw(3) << '*' << ' ' ; }
                        // inside the aligned region
                        else
                        {   o << std::setprecision(2) << std::setw(3) << (this->data(i,j_dat_fw) < 1e-10 ? 0 : this->data(i,j_dat_fw)) << ' ' ; }
                    }
                }
                o << std::endl ;
            }
        }
    }
}


void KmeanEngine::updateClusterSizes()
{   this->resetClusterSizes() ;

    for(size_t i=0; i<this->clusters.size(); i++)
    {   size_t cluster = this->clusters[i] ;
        this->cluster_sizes[cluster] += 1  ;
    }
}

bool KmeanEngine::hasConverged() const
{   // obviously before the second iteration
    // convergence cannot be achieved (no previous results)
    if(this->n_iter < 1)
    {   return false ; }

    size_t n = this->data.get_nrow() ;

    bool converged = true ;
    for(size_t i=0; i<n; i++)
    {   if(this->clusters[i]  != this->clusters_previous[i] or
           this->shifts[i]    != this->shifts_previous[i] or
           this->distances[i] != this->distances_previous[i])
         {  converged = false ;
            break ;
         }
         if(this->flip and this->flips[i] != this->flips_previous[i])
         {  converged = false ;
            break ;
         }
    }
    return converged ;
}
