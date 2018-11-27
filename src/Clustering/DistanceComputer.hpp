#ifndef DISTANCECOMPUTER
#define DISTANCECOMPUTER

#include <iostream>
#include <vector>
#include <stdexcept>   // std::invalid_argument

/*!
 * \brief The DistanceComputer class is an abstract class providing an interface to
 * be implemented by any sub-class which needs to compute distances between vectors
 * allowing shifting and flipping.
 *
 * Shifting means that each vector will be decomposed into a set of smaller
 * sub-parts of equal length. For instance, a shifting of 5 means that each vector
 * will be decomposed into 5 smaller sub-parts of equal lengths. A shift of 1 means
 * that the entire vectors are considered. Then, all the possible pair-wise comparisons
 * between both sets will be performed and the comparison resulting in the smallest
 * distance will be considered as the optimal comparison and reported.
 * Flipping means that, in addition to all the pairwis comparisons described above,
 * extra comparisons with the y sub-parts being in reversed orientation will be
 * performed. In this case, the returned distance is also the minimal distance computed
 * from all the possible pair-wise comparisons.
 *
 * To make computations faster, any distance computation allowing shifting and flipping
 * be solved using a dynamic programming solution to avoid recomputing more than once
 * the same thing. However, this requires to use buffers.
 * A problem arise when the need for distance computation grows. The buffers need to be
 * reallocated each time. This operation can become quiet cumbersome. This problem can
 * be solved by the usage of an object computing the distances between any two vectors
 * through a dedicated method (defined in this interface) and internally taking care of
 * the buffer state and reallocating them on need only, not every time. Typically the
 * buffers needs to be reallocated when the length of the vectors involved in the
 * computations changes from one call to the next, but not otherwise.
 * This class offers an interface to other classes dedicated to this purposes.
 * The buffer sizes depends on the shifting value only and should only be reallocated
 * when the shifting freedom changes since the last computations.
 */
class DistanceComputer
{
    public:
        // constructors and destructor
        /*!
         * \brief Destructor.
         */
        virtual ~DistanceComputer() ;

        // methods
        /*!
         * \brief Should implement the computation of the distance between two
         * vectors x and y of equal lengths, allowing up to <shift> freedom of
         * shifting (in number of vector elements).
         * \param x the first vector.
         * \param y the second vector
         * \param shift the freedom of shifting (in number of vector elements).
         * \param shift_x reports the index of the element of x which is the first element
         * of the sub-part of x used in the optimal comparison.
         * \param shift_y reports the index of the element of y which is the first element
         * of the sub-part of y used in the optimal comparison.
         * \return the distance between x and y.
         * \throw std::invalid_argument if x and y are not of equal length or if the
         * shift is too high given the vector lengths.
         */
        virtual double compute_distance(const std::vector<double>& x,
                                        const std::vector<double>& y,
                                        size_t shift,
                                        size_t& shift_x,
                                        size_t& shift_y) throw (std::invalid_argument) = 0 ;

        /*!
         * \brief Should implement the computation of the distance between two
         *  vectors x and y of equal lengths, allowing up to <shift> freedom of
         * shifting (in number of vector elements) and flipping.
         * \param x the first vector.
         * \param y the second vector
         * \param shift the freedom of shifting (in number of vector elements).
         * \param shift_x reports the index of the element of x which is the first element
         * of the sub-part of x used in the optimal comparison.
         * \param shift_y reports the index of the element of y which is the first element
         * of the sub-part of y used in the optimal comparison.
         * \param flip_y reports whether the slice of <y> was reversed in the optimal
         * comparison. If true, it simply means that y[shift_y] was the last element of
         * the slice instead of the first.
         * \return the distance between x and y.
         * \throw std::invalid_argument if x and y are not of equal length or if the
         * shift is too high given the vector lengths.
         */
        virtual double compute_distance(const std::vector<double>& x,
                                        const std::vector<double>& y,
                                        size_t shift,
                                        size_t& shift_x,
                                        size_t& shift_y,
                                        bool& flip_y) throw (std::invalid_argument) = 0 ;


    protected:
        // methods
        /*!
         * \brief This method should reallocate all the buffers given the new shift value.
         * \param shift the new shift value.
         * \return whether a reallocation happened.
         * \throw std::invalid_argument if shift is 0.
         */
        virtual bool reallocate_buffers(size_t shift) throw (std::invalid_argument) = 0 ;
} ;


#endif // DISTANCECOMPUTER
