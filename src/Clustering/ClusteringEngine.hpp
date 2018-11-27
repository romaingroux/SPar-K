#ifndef CLUSTERINGENGINE_HPP
#define CLUSTERINGENGINE_HPP

#include <iostream>
#include <string>
#include <vector>
#include <thread>
#include <future>
#include <utility>
#include "Matrix/Matrix.hpp"
#include "Matrix/Matrix2D.hpp"
#include "Utility/Constants.hpp"

/*!
 * \brief The ClusteringEngine class is an abstract class providing an interface
 * to other classes implementing data clustering methods.
 */
class ClusteringEngine
{
    public:
        /*!
         * \brief Destructor.
         */
        virtual ~ClusteringEngine() ;

        // ------------------------- getters ---------------------------
        /*!
         * \brief Should return the cluster references (reprensentations).
         * The references should be contained in the rows. For a K-means, it
         * should return the cluster centroids.
         * \return The current cluster references.
         */
        virtual Matrix2D<double> getReferences() const ;

        // ------------------- clustering interface ---------------------
        /*!
         * \brief Should assign the data to the clusters. This is the main
         * functionality of this class. For a K-means, which is an iterative
         * partitioning method, this method should run one iteration at the time
         * (and thus this method should be called repetitively).
         * \return a code indicating whether this call resulted in a successful
         * data assignation, whether it failed for some reason or whether it
         * converged (for an iterative partitioning, this means that calling this
         * method an other time won't change anything, the results are stable).
         */
        virtual Constants::clustering_codes cluster() = 0 ;

        // ------------------ results display method ---------------------
        /*!
         * \brief This should send a representation of the results (whatever this means)
         * on the given stream.
         * \param o the stream of interest.
         */
        virtual void print_results(std::ostream& o) = 0 ;


    protected:

        // --------------------- seeding methods ------------------------
        /*!
         * \brief This should set the initial cluster references using the method
         * corresponding to the given argument.
         * \param method a string indicating one method among several to use. For
         * instance "random" which would lead to a random assignment (whatherver
         * this means) of the references.
         * \throw std::runtime_error if the method is not recognized.
         */
        virtual void seeding(const std::string& method) throw(std::runtime_error) = 0 ;

        // --------------------------- fields ---------------------------

        /*!
         * \brief the data to cluster
         */
        Matrix2D<double> data ;

        /*!
         * \brief the current cluster references (medoids, centroids, whatever).
         * There should be one reference per row. The number of column should
         * depend on the reference dimensionality.
         */
        Matrix2D<double> references ;

        /*!
         * \brief the current iteration number (the number of times cluster()
         * was called so far).
         */
        size_t n_iter ;

        /*!
         * \brief whether flipping is on.
         */
        bool flip ;

        /*!
         * \brief The number of shift states.
         */
        size_t n_shift ;


} ;

#endif // CLUSTERINGENGINE_HPP
