
#ifndef CORDISTANCECOMPUTER_HPP
#define CORDISTANCECOMPUTER_HPP

#include "DistanceComputer.hpp"

#include <vector>
#include <iostream>

// a few typedef
using vector_double   = std::vector<double> ;
using matrix2d_double = std::vector<vector_double> ;
using matrix3d_double = std::vector<matrix2d_double> ;
using vector_int      = std::vector<int> ;
using matrix2d_int    = std::vector<vector_int> ;


/*!
 * \brief The CorDistanceComputer class is a class allowing to rapidly compute
 * the correlation distance with shift and flip between any two vectors x and y.
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
 * Here, the distance computed corresponds to : d(x,y) = min{1-cor(x',y')}, that is the
 * minimal distance obtained by comparing any two slices x' and y' of x and y. If several
 * comparisons results in the same optimal distance, the most central
 * comparison is preferred (each slice comparison is assigned a score using :
 * abs(l/2 - (l/2 + s1)) + abs(l/2 - (l/2 + s2) where <l> is the vector size and <s1> and
 * <s2> are the 1st positions of each slice) corresponding to the absolute difference of
 * position of the aligned slice to the middle of the complete vectors.
 *
 * The computations are made fast by 1) implementing a dynamic programming solution and
 * 2) allocating buffers only when it is necessary, not every time a new distance is
 * computed. The buffers required by the dynamic programming approach are managed
 * internally by the instance and are not reallocated each time a distance between a
 * pair of vectors is computed. The buffers are reallocated only when the shift value
 * changes from one call to compute_distance() to the next but not otherwise.
 */
class CorDistanceComputer : public DistanceComputer
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
        CorDistanceComputer(size_t shift=1) throw (std::invalid_argument) ;

        /*!
         * \brief Destructor.
         */
        virtual ~CorDistanceComputer() override ;

        // methods
        /*!
         * \brief Computes the correlation distance between two vectors x and y
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
                                        size_t& shift_y) throw (std::invalid_argument) override ;

        /*!
         * \brief Computes the correlation distance between two vectors x and y
         * of equal lengths, allowing up to <shift> freedom of shifting (in number
         * of vector elements) and flipping. The minimal distance computed (that is,
         * in the optimal comparison) from all the comparisons is returned and the
         * 1st position of each slice used in the optimal comparison are reported
         * using <shift_x> and <shift_y> and the flipping state of the y slice (that
         * is whether this slice was flipped in the comparison) using <flip_y>.
         * In case where a correlation between both vectors cannot be computed
         * (for instance one vector has a variance of 0), the distance is set
         * to 1 (corresponding to a correlation of 0).
         * \param y the second vector
         * \param shift the allowed shifting freedom (in number of vector elements)
         * \param shift_x the index of the starting element in x used in the optimal
         * the optimal comparison.
         * \param shift_y the index of the starting element in y used in the optimal
         * the optimal comparison.
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
                                        bool& flip_y) throw (std::invalid_argument) override ;

    protected:
        // methods
        /*!
         * \brief This method should reallocate all the buffers given the new shift
         * value.
         * \param shift the new shift value.
         * \return whether a reallocation happened.
         * \throw std::invalid_argument if shift is 0.
         */
        bool reallocate_buffers(size_t shift) throw (std::invalid_argument) ;

        /*!
         * \brief This is the internal routine computing the correlation distance
         * between two vectors x and y of equal lengths, allowing up to <shift>
         * number of shifting states and flippin. This function performs all
         * possible comparisons (all possible alignments) between each slice of
         * x and y of length (l-shift+1) where <l> is the vector length
         * and <shift> the value of the corresponding parameter. Each comparison is
         * scored with a distance of 1-cor(x',y') where x' and y' are two slices
         * of x and y.
         * The minimal distanshift the freedom of shifting (in number of vector elements).
         * \param shift_x reports the index of the element of x which is the first element
         * of the sub-part of x used in the optimal comparison.
         * \param shift_y reports the index of the element of y which is the first element
         * of the sub-part of y used in the optice computed (that is, in the optimal comparison) from all
         * the comparisons is returned and the parameters of the optimal alignment are
         * reported using <shift_x>, <shift_y> and <flip_y>.
         * If several alignments results in the same optimal distance, the most central
         * alignment is preferred (each alignment is assigned a score using :
         * abs(l/2 - (l/2 + s1)) + abs(l/2 - (l/2 + s2) where <l> is the vector size and
         * <s1> and <s2> are the 1st positions of each slice) corresponding to the
         * absolute difference of position of the aligned slice to the middle of the
         * complete vectors.
         * In case where a correlation between both vectors cannot be computed
         * (for instance one vector has a variance of 0), the distance is set
         * to 1 (corresponding to a correlation of 0).
         * The two public methods compute_distance (see above) are wrappers around
         * this method.
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
                               bool& flip_y) throw (std::invalid_argument) ;

        // static methods
        /*!
         * \brief This method corrects (modify) a correlation value in three cases :
         * case 1) rounding to 0, a value is 0 +/-2xepsilon machine because of rounding errors,
         * case 2) rounding to 1, a value is 1 +/-2xepsilon machine because of rounding errors,
         * case 3) turns to Constants::corr_fail because it is a NaN value.
         * \param x the correlation to check for correction.
         * \return whether the value has been corrected.
         */
        static bool correct_cor(double& x) ;

        /*!
         * \brief Finds the outliers in a given vector - any value which is
         * bigger/smaller than the vector mean plus/minus three standard
         * deviations - and replaces these values by the mean +/- three standard
         * deviations (if the value is an upper or lower outlier respectively).
         * \param x a vector of interest.
         * \return a copy of the orginal vector with corrected outliers.
         * \throw std::invalid_argument if the size of the vector is < 2.
         */
        static std::vector<double> replace_outliers(const std::vector<double>& x) throw (std::invalid_argument) ;


        // fields
        /*!
         * \brief The current shift value which this
         * instance can compute distances for (all the buffers are
         * allocated to fit this shifting freedom).
         */
        size_t _shift ;

        /*!
         * \brief stores the distance 1-cor(x',y') between any
         * two slices x' and y' of two vectors x and y, in both
         * forward and reverse (if flip is allowed).
         */
        matrix3d_double _distances ;

        /*!
         * \brief stores the partial sums of x*y for two
         * x and y vectors.
         */
        matrix3d_double _sum_xy ;

        /*!
         * \brief stores the partial sum of x.
         */
        vector_double _sum_x ;

        /*!
         * \brief stores the partial sum of squares of x.
         */
        vector_double _sum_x2 ;

        /*!
         * \brief stores the partial sum of y.
         */
        vector_double _sum_y ;

        /*!
         * \brief stores the partial sum of squares y.
         */
        vector_double _sum_y2 ;

        /*!
         * \brief stores scores related to how central the
         * slices used for each comparison were.
         */
        matrix2d_int _scores ;
} ;


#endif // CORDISTANCECOMPUTER_HPP
