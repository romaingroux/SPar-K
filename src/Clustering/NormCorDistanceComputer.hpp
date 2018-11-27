#ifndef NORMCORDISTANCECOMPUTER_HPP
#define NORMCORDISTANCECOMPUTER_HPP

#include "CorDistanceComputer.hpp"

#include <vector>
#include <stdexcept>  // invalid_argument
#include <cmath>      // pow(), sqrt()
#include <limits>     // numeric_limits<double>::max(), numeric_limits<int>::max()
#include <Utility/Utility.hpp>       // isNaN()
#include <Utility/Constants.hpp>     // Constants

/*!
 * \brief The CorDistanceComputer class is a class allowing to rapidly compute
 * a normalised correlation distance with shift and flip between any two vectors
 * x and y.
 *
 * Shifting means that each vector will be decomposed into a set of smaller
 * slices of equal length. More formally, each vector will be decomposed into
 * <shift> slices of length (<l>-<shift>+1) where <l> is the vector length. For
 * instance, a shifting of 5 means that each vector will be decomposed into 5 smaller
 * slices of equal lengths. A shift of 1 means that the entire vectors are considered.
 * Then, all the possible pair-wise comparisons between both sets will be performed
 * and the comparison resulting in the smallest distance will be considered as the
 * optimal comparison and reported.
 * Flipping means that, in addition to all the pairwis comparisons described above,
 * extra comparisons with the y slices being in reversed orientation will be
 * performed. In this case, the returned distance is also the minimal distance computed
 * from all the possible pair-wise comparisons.
 *
 * Here the distance computed corresponds a normalized version of the correlation
 * distance : d(x,y) = min{1-cor(x',y') / sum(all other comparison distances but the one
 * of the current slices x' and y')}. If several comparisons results in the same optimal
 * distance, the most central comparison is preferred (each slice comparison is assigned
 * a score using : abs(l/2 - (l/2 + s1)) + abs(l/2 - (l/2 + s2) where <l> is the vector
 * size and <s1> and <s2> are the 1st positions of each slice) corresponding to the
 * absolute difference of position of the aligned slice to the middle of the complete
 * vectors.
 *
 * The computations are made fast by 1) implementing a dynamic programming solution and
 * 2) allocating buffers only when it is necessary, not every time a new distance is
 * computed. The buffers required by the dynamic programming approach are managed
 * internally by the instance and are not reallocated each time a distance between a
 * pair of vectors is computed. The buffers are reallocated only when the shift value
 * changes from one call to compute_distance() to the next but not otherwise.
 */
class NormCorDistanceComputer : public CorDistanceComputer
{
    public:
        // constructors and destructor
        /*!
         * \brief Defaut constructor.
         * Constructs a instance able to compute distances between vectors
         * and allocate the buffers such that they can handle distance
         * computations for the given shifting freedom value. It the shift
         * value at the next call of compute_distance() changes, then the
         * buffers will be reallocated.
         * \param vector_size the expected vector size.
         * \param shift the expected shifting freedom.
         * \throw std::invalid_argument if shift is 0.
         */
        NormCorDistanceComputer(size_t shift=1) throw (std::invalid_argument) ;

        /*!
         * \brief Destructor.
         */
        virtual ~NormCorDistanceComputer() override ;

    protected:
        // methods
        /*!
         * \brief Overrides the internal computation routine compute() of the
         * CorDistanceComputer class.
         * Computes the correlation distance between two vectors x and y
         * of equal lengths, allowing up to <shift> freedom of shifting (in number
         * of vector elements). The minimal distance computed (that is, in the
         * optimal comparison) from all the comparisons is returned and the 1st
         * position of each slice used in the optimal comparison are reported
         * using <shift_x> and <shift_y>.
         * In case where a correlation between both vectors cannot be computed
         * (for instance one vector has a variance of 0), the distance is set
         * to 1 (corresponding to a correlation of 0).
         * \param x the first vector.
         * \param y the second vector
         * \param shift the allowed shifting freedom (in number of vector elements).
         * \param shift_x reports the index of the element of x which is the first element
         * of the sub-part of x used in the optimal comparison.
         * \param shift_y reports the index of the element of y which is the first element
         * of the sub-part of y used in the optison.
         * \param flip_y serves a double purpose. First, it is used to instruct the
         * method whether flipping is enabled (true) or not (false). Second, it is used
         * to return  whether the slice of <y> was reversed in the optimal alignment. If
         * true (this only happens if the value is true when calling the function), it
         * simply means that y[shift_y] was the last element of the slice instead of the
         * first.
         * \return the distance between x and y.
         * \throw std::invalid_argument if x and y are not of equal length or if the
         * shift is too high given the vector lengths.
         */
        virtual double compute(const std::vector<double>& x,
                               const std::vector<double>& y,
                               size_t shift,
                               size_t& shift_x,
                               size_t& shift_y,
                               bool& flip_y) throw (std::invalid_argument) override ;

} ;

#endif // NORMCORDISTANCECOMPUTER_HPP
