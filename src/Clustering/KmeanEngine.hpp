#ifndef KMEANENGINE_HPP
#define KMEANENGINE_HPP

#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <future>
#include <utility>
#include <unordered_map>
#include <Matrix/Matrix.hpp>
#include <Matrix/Matrix2D.hpp>
#include <Utility/Constants.hpp>
#include <Clustering/ClusteringEngine.hpp>
#include <Clustering/DistanceComputer.hpp>
#include <Parallel/ThreadPool.hpp> // ThreadPool

/*!
 * \brief The KmeanEngine class is the support to cluster the data using the K-mean algorithm.
 * Basically, the clustering will do two things at the same time :
 * i)  partition the data into the given number of clusters
 * ii) re-align the data inside each cluster in such a way that most ressembling part of the datum
 * overlap each other.
 * Thus, this class implements an algorithm which tries to partition the data into K alignments.
 * Each cluster/alignment is represented by a vector containing the mean signal at each position along
 * the alignment.
 * To cluster a dataset, a KmeanEngine should be instanciated and cluster() should called as
 * many time as needed/wanted. Then the results can be extracted using the relevant getters.
 */
class KmeanEngine : public ClusteringEngine
{
    public:
        // constructors and destructor
        KmeanEngine() = delete ;
        KmeanEngine(const KmeanEngine& other) = delete ;

        /*!
         * \brief Constructs a classifier to partition a dataset into a given number
         * of clusters. References are initialized using the given method and seed.
         * \param data the matrix of data to cluster. The datum (or points) should
         * be contained on the rows. The column number is the data dimensionality.
         * \param cluster_n the number of clusters to parition the data into.
         * \param shift_n the freedom of shifting in number of columns of the data
         * matrix when clustering. For instance, 5 means that each row of the data
         * matrix will be decomposed into 5 smaller sub-parts of equal lengths
         * which will be used for the comparisons with the other data rows during the
         * clustering.
         * \param flip enables flipping when clustering meaning that each sub-part
         * can be reverted during the comparisons.
         * \param seeding_method a method to initialise the references. It should be
         * "random" or "kmean++".
         * "random_" will sample k datum as the initial referenes, with uniform
         * probabilities (by default).
         * "kmean++" selects k datum using the kmean++ algorithm.
         * \param seed a sequence to seed the random number generator. If seed = ""
         * then the random number generator is not initialized using any seed.
         * \param n_threads the number of threads dedicated to the computations.
         * \param dist_type the distance to use.
         * \throw std::invalid argument if the shifting freedom is bigger than
         * the data number of column - 1 or if the number of clusters is not
         * smaller than the data number of rows.
         */
        KmeanEngine(const Matrix2D<double>& data,
                    size_t cluster_n, size_t shift_n, bool flip,
                    const std::string& seeding_method="random",
                    const std::string& seed="",
                    size_t n_threads=1,
                    Constants::dist_type dist_type=Constants::COR_DIST)
        throw (std::invalid_argument) ;

        /*!
         * \brief Constructs a classifier to partition a dataset into a given number
         * of clusters. References are initialized using the given method and seed.
         * \param path_data the address of the file containing the data matrix to cluster.
         * This file should be a text file containing a numerical matrix. The datum
         * (or points) should be contained on the rows. The column number is the data
         * dimensionality.
         * \param shift_n the freedom of shifting in number of columns of the data
         * matrix when clustering. For instance, 5 means that each row of the data
         * matrix will be decomposed into 5 smaller sub-parts of equal lengths
         * which will be used for the comparisons with the other data rows during the
         * clustering.
         * \param flip enables flipping when clustering meaning that each sub-part
         * can be reverted during the comparisons.
         * \param seeding_method a method to initialise the references. It should be
         * "random" or "kmean++".
         * "random_" will sample k datum as the initial referenes, with uniform
         * probabilities (by default).
         * "kmean++" selects k datum using the kmean++ algorithm.
         * \param seed a sequence to seed the random number generator. If seed = ""
         * then the random number generator is not initialized using any seed.
         * \param n_threads the number of threads dedicated to the computations.
         * \param dist_type the distance to use.
         * \throw std::invalid argument if the shifting freedom is bigger than
         * the data number of column - 1 or if the number of clusters is not
         * smaller than the data number of rows.
         */
        KmeanEngine(const std::string& path_data,
                    size_t cluster_n, size_t shift_n, bool flip,
                    const std::string& seeding_method="random",
                    const std::string& seed="",
                    size_t n_threads=1,
                    Constants::dist_type dist_type=Constants::COR_DIST)
        throw (std::invalid_argument) ;

        /*!
         * \brief Constructs a classifier to partition a dataset into a given number
         * of cluster. Contrarly to the previous constructors, this one allows to
         * specify the cluster references explicitly by providing a matrix containing
         * them (thus no seeding is required). The number of cluster of the partition
         * is defined by the number of references provided.
         * \param data the matrix of data to cluster. The datum (or points) should
         * be contained on the rows. The column number is the data dimensionality.
         * \param references the matrix containing the cluster references. The references
         * should be on the rows. The column number should fit the data dimensionality.
         * \param shift_n the freedom of shifting in number of columns of the data
         * matrix when clustering. For instance, 5 means that each row of the data
         * matrix will be decomposed into 5 smaller sub-parts of equal lengths
         * which will be used for the comparisons with the other data rows during the
         * clustering.
         * \param flip enables flipping when clustering meaning that each sub-part
         * can be reverted during the comparisons.
         * \param n_threads the number of threads dedicated to the computations.
         * \param dist_type the distance to use.
         * \throw std::invalid argument if the shifting freedom is bigger than
         * the data number of column - 1 or if the number of clusters is not
         * smaller than the data number of rows.
         */
        KmeanEngine(const Matrix2D<double>& data,
                    const Matrix2D<double>& references,
                    size_t shift_n, bool flip,
                    size_t n_threads,
                    Constants::dist_type dist_type=Constants::COR_DIST)
        throw (std::invalid_argument) ;


        /*!
         * \brief Constructs a classifier to partition a dataset into a given number
         * of cluster. Contrarly to the previous constructors, this one allows to
         * specify the cluster references explicitly by providing a matrix containing
         * them (thus no seeding is required). The number of cluster of the partition
         * is defined by the number of references provided.
         * \param path_data the address of the file containing the data matrix to cluster.
         * This file should be a text file containing a numerical matrix. The datum
         * (or points) should be contained on the rows. The column number is the data
         * dimensionality.
         * \param path_references the address of the file containing the cluster references.
         * This file should be a text file containing a numerical matrix. The references
         * should be contained on the rows. The column number should fit the data
         * dimensionality.
         * \param shift_n the freedom of shifting in number of columns of the data
         * matrix when clustering. For instance, 5 means that each row of the data
         * matrix will be decomposed into 5 smaller sub-parts of equal lengths
         * which will be used for the comparisons with the other data rows during the
         * clustering.
         * \param flip enables flipping when clustering meaning that each sub-part
         * can be reverted during the comparisons.
         * \param n_threads the number of threads dedicated to the computations.
         * \param dist_type the distance to use.
         * \throw std::invalid argument if the shifting freedom is bigger than
         * the data number of column - 1 or if the number of clusters is not
         * smaller than the data number of rows.
         */
        KmeanEngine(const std::string& path_data,
                    const std::string& path_references,
                    size_t shift_n, bool flip,
                    size_t n_threads,
                    Constants::dist_type dist_type=Constants::COR_DIST)
        throw (std::invalid_argument) ;

        // methods
        /*!
         * \brief Destructor. Take care of joining the threads.
         */
        virtual ~KmeanEngine() override ;

        // ------------------------- getters ---------------------------
        /*!
         * \brief Returns a vector containing the distance of each datum
         * to its cluster reference.
         * \return The current distances of each datum to its cluster.
         */
        std::vector<double> getDistances() const ;

        /*!
         * \brief Returns the cluster assignment for each datum.
         * \return The current cluster assignments.
         */
        std::vector<size_t> getClusters() const ;

        /*!
         * \brief Returns two vectors (a vector of two vectors)
         * containing the shift values. The first vector contains
         * the reference shift values. The second vector contains
         * the data shift values.
         * \return The current shift values.
         */
        std::vector<std::vector<size_t>> getShifts() const ;

        /*!
         * \brief Returns a vector containg the flip states for
         * each datum.
         * \return The current flip values.
         */
        std::vector<int> getFlips() const ;

        // defined and implemented at the ClusteringEngine level
        // Matrix2D<double> getReferences()

        // ------------------- clustering interface ---------------------
        /*! \brief Runs one clustering iteration. It assignes the datum to
         * the nearest cluster and updates the references, clusters, shifts,
         * flips, distances and cluster sizes.
         * \return A code indicating whether the clustering succeeded.
         * A failure is defined when there is a cluster without any datum a
         * assigned to it. A convergence is defined as the case when
         * the clusters, the shifts and the flips assigmnents remain stable
         * since the last iteration (thus this is impossible before the 2nd
         * iteration).
         */
        virtual Constants::clustering_codes cluster() override ;

        // ------------------ results display method ---------------------
        /*!
         * \brief sends a nicely formatted matrix containing the current
         * cluster assignment, shift and flip values to the given stream.
         * \param o the stream on which the representation should be sent.
         */
        virtual void print_results(std::ostream& o) override ;

    private:
        // methods
        // -------------------------- thread related --------------------------
        /*!
         * \brief Joins the threads in the pool.
         */
        void joinThreads() ;

        /*!
         * \brief Gets the thread pool dedicated to the computations.
         * \return a reference to the thread pool.
         */
         ThreadPool& getThreads() ;

        // --------------------- seeding methods ------------------------
        /*!
         * \brief An interface to the seeding methods. Initialises the references
         * using the method corresponding to the given argument.
         * \param method the method which should be used to initialise
         * the references : "random_simple", "random", "random_shift_flip",
         * "kmean++". To know more about these methods, look at the
         * corresponding private methods.
         * \throw std::runtime_error if the method is not recognized.
         */
       virtual void seeding(const std::string& method) throw(std::runtime_error) override ;

        // ---------------- constructor related method -------------------
        /*!
         * \brief This method is called by the constructor with the same signature. It
         * is a construction routine containing all the instructions to construct the
         * instance.
         * \param data the matrix of data to cluster. The datum (or points) should
         * be contained on the rows. The column number is the data dimensionality.
         * \param cluster_n the number of clusters to parition the data into.
         * \param shift_n the freedom of shifting in number of columns of the data
         * matrix when clustering. For instance, 5 means that each row of the data
         * matrix will be decomposed into 5 smaller sub-parts of equal lengths
         * which will be used for the comparisons with the other data rows during the
         * clustering.
         * \param flip enables flipping when clustering meaning that each sub-part
         * can be reverted during the comparisons.
         * \param seeding_method a method to initialise the references. It should be
         * "random" or "kmean++".
         * \param seed a sequence to seed to random number generator.
         * \throw std::invalid argument if the number of clusters is not smaller than
         * the data number of rows.
         */
        void constructorInit(const Matrix2D<double>& data,
                             size_t cluster_n, size_t shift_n, bool flip,
                             const std::string& seeding_method,
                             const std::string& seed,
                             Constants::dist_type dist_type)
        throw (std::invalid_argument) ;

        /*!
         * \brief This method is called by the constructor with the same signature. It
         * is a construction routine containing all the instructions to construct the
         * instance.
         * \param data the matrix of data to cluster. The datum (or points) should
         * be contained on the rows. The column number is the data dimensionality.
         * \param cluster_n the number of clusters to parition the data into.
         * \param references the matrix containing the cluster references. The references
         * should be on the rows. The column number should fit the data dimensionality.
         * \param shift_n the freedom of shifting in number of columns of the data
         * matrix when clustering. For instance, 5 means that each row of the data
         * matrix will be decomposed into 5 smaller sub-parts of equal lengths
         * which will be used for the comparisons with the other data rows during the
         * clustering.
         * \param flip enables flipping when clustering meaning that each sub-part
         * can be reverted during the comparisons.
         * \param dist_type the distance to use.
         * \throw std::invalid argument if the shifting freedom is bigger than
         * the data number of column - 1 or if the number of clusters is not
         * smaller than the data number of rows.
         */
        void constructorInit(const Matrix2D<double>& data,
                             const Matrix2D<double>& references,
                             size_t shift_n, bool flip,
                             Constants::dist_type dist_type)
        throw (std::invalid_argument) ;

        /*!
         * \brief A construction routine called by all constructors, containing their
         * common parts. It initializes the data, the clusters assignment values,
         * the cluster sizes, the shift values, the distance values, the flip values
         * and the number of time cluster() was called
         * \param data the matrix of data to cluster.
         * \param cluster_n the number of cluster.
         * \param shift_n the freedom of shifting in number of columns of the data
         * matrix when clustering. For instance, 5 means that each row of the data
         * matrix will be decomposed into 5 smaller sub-parts of equal lengths
         * which will be used for the comparisons with the other data rows during the
         * clustering.
         * \param flip enables flipping when clustering meaning that each sub-part
         * can be reverted during the comparisons.
         * \param di\param shift_n the freedom of shifting in number of columns of the data
         * matrix when clustering, for instance 5 will be interpretated as
         * -2,-1, 0,+1,+2 possible shift states.
         * \param flip whether flipping freedom is enabled.st_type the distance to use.
         */
        void constructorRoutine(const Matrix2D<double>& data,
                                size_t cluster_n, size_t shift_n, bool flip,
                                Constants::dist_type dist_type) ;

        /*!
         * \brief A constructor routine dynamically allocating the objects in charge of the
         * distance computations, according to the wanted distance.
         * \param dist_type the distance to use.
         * \throw std::invalid_argument if the distance type is not recognized.
         */
        void allocateDistComputers(Constants::dist_type dist_type) throw (std::invalid_argument);

        // ---------------- reference seeding methods -------------------
        /*!
         * \brief Selects k (known because given to the constructor) datum randomly
         * as references using a uniform probability distribution.
         */
        void seeding_random() ;

        /*!
         * \brief Selects k (known because given to the constructor) datum as
         * references using the K-means++ algorithm.
         */
        void seeding_kmean_pp() ;        

        // ---------------- clustering related methods ------------------

        /*!
         * \brief Checks whether at least one cluster is empty.
         * \return Whether at least one cluster is empty.
         */
        bool hasEmptyCluster() const ;

        /*!
         * \brief Resets the cluster sizes to zero values.
         */
        void resetClusterSizes() ;

        /*!
         * \brief Implements the clustering routine. It runs the data cluster
         * assignment for one iteration. More specifically it assigns datum (rows)
         * [from,to) to the least dissimilar cluster (by comparison with the
         * cluster references). Then it updates the clusters assignment, shifts,
         * flip (if on) and distances cluster size values.
         *
         *  This method is made to be run in threads by providing to several
         * concurrent threads non-overlapping [from,to) values (slices of the
         * data).
         *
         * \param from the index of first row (included) of data which should
         * be considered.
         * \param to the index of the row AFTER (excluded) the last row to be
         * considered.
         * \param dc a pointer to the object in charge of computing the distances.
         * \param prms_done a promise to be filled instead to instruct the calling
         * thread that the method returned (instead of directly returning, allowing
         * thread synchronization).
         */
         void assign_to_group(size_t from,
                              size_t to,
                              DistanceComputer* dc,
                              std::promise<bool>& prms_done) ;

        /*!
         * \brief Replaces the current reference patterns with the content of the given
         * futures. Each futures is expected to contain one entire reference pattern.
         * The first future is expected to contain the reference pattern of the first
         * cluster, the second future the reference pattern of the second cluster and
         * so on...
         * \param v_ftr_ref a vector of futures containing each cluster reference
         * pattern.
         */
        void updateReferences(std::vector<std::future<std::vector<double>>>& v_ftr_ref) ;

        /*!
         * \brief Computes the given cluster reference given the current cluster
         * assignment, shift, flip and cluster size values. This method realign the
         * data assigned to a given cluster and aggregate them to obtain the mean
         * signal for the cluster.
         *
         * This method is made to be run in threads by providing concurrent threads
         * different cluster indexes.
         *
         * \param cluster_id the label of the cluster of interest (it will be
         * used to match the cluster assignment values directly).
         * \param prms_ref a promise to be filled at the end of the execution
         * with the computed cluster reference.
         * \throw std::invalid_argument if the cluster_id is not recognized, that
         * is it is bigger than the number of cluster - 1 (cluster labels are
         * 0-based).
         */
        void computeClusterReference(size_t cluster_id,
                                     std::promise<std::vector<double>>& prms_ref)
        throw (std::invalid_argument) ;

        /*!
         * \brief Sends a representation of the data alignment (with data and references)
         * given the current cluster assignment, shift, flip values to the given stream.
         * This method is designed for debugging purposes mostly.
         */
        void print_alignment(std::ostream& o) const ;


        /*!
         * \brief Computes the current cluster sizes.
         */
        void updateClusterSizes() ;

       /*!
        * \brief Checks whether convergence is achieved. Convergence is achieved
        * when the cluster assignment, shift and flip values are the same as those
        * at of the previous iteration.
        * \return whether convergence is achieved.
        */
        bool hasConverged() const ;

        // fields
        /*!
         * \brief the distance of each datum to its cluster reference.
         */
        std::vector<double> distances ;

        /*!
         * \brief the distance of each datum to its cluster reference at the previous
         * iteration (the previous call of cluster(), obviously not before the 2nd
         * iteration).
         */
        std::vector<double> distances_previous ;

        /*!
         * \brief The current data cluster assignment values.
         */
        std::vector<size_t> clusters ;

        /*!
         * \brief The previous data cluster assignment values (at the previous call of
         * cluster(), obviously not before the 2nd iteration)
         */
        std::vector<size_t> clusters_previous ;

        /*!
         * \brief The current cluster sizes.
         */
        std::vector<double> cluster_sizes ;

        /*!
         * \brief The current data shifting values. The first vector contains the
         * reference shift values and the second the data shift values.
         */
        std::vector<std::vector<size_t>> shifts ;
        /*!
         * \brief The previous data shifting values (at the previous call of cluster(),
         * obviously not before the 2nd iteration). The first vector contains the
         * reference shift values and the second the data shift values.
         */
        std::vector<std::vector<size_t>> shifts_previous ;

        /*!
         * \brief The current data flip values.
         */
         std::vector<int> flips ;
        /*!
         * \brief The previous data flip values (at the previous call of cluster(),
         * obviously not before the 2nd iteration)
         */
        std::vector<int> flips_previous ;

        /*!
         * \brief The threads dedicated to the computations.
         */
        ThreadPool threads ;

        /*!
         * \brief a vector of pointers to objects in charge of computing distances
         * between any two vectors.
         */
        std::vector<DistanceComputer*> dist_computers ;

        // defined in ClusteringEngine
        // Matrix2D<double> data ;
        // Matrix2D<double> references ;
        // size_t  n_shift ;
        // bool flip ;
        // size_t n_iter ;
} ;

#endif // KMEANENGINE_HPP
