#ifndef DISTANCES_HPP
#define DISTANCES_HPP

#include <iostream>
#include <vector>
#include <Utility/Constants.hpp>



/*!
 * \brief Compute the correlation distance between 2 vectors x
 * and y. The distance is defined as 1-cor(x,y). x and y
 * should have the same length.
 * If no correlation can be computed (because one vector has a
 * standard deviation equal to 0) then Constants::corr_fail
 * is returned (see Utility/Constants).
 * \param x the first vector.
 * \param y the second vector.
 * \return the distance
 */
double cor_distance(const std::vector<double>& x,
                    const std::vector<double>& y) ;


/*!
 * \brief Computes the correlation distance between two vectors
 * x and y of equal lengths, allowing up to <shift> number of
 * shifting states. This function performs all possible
 * comparisons (all possible alignments) between each slice of
 * x and y of length (l-shift+1) where <l> is the vector length
 * and <shift> the value of the corresponding parameter.
 * The minimal correlation distance is returned and the parameters
 * of the optimal alignment are reported using <shift1>, and
 * <shift2>.
 * If several alignments results in the same optimal
 * distance, the most central alignment is preferred (each alignment
 * is assigned a score using : abs(l/2 - (l/2 + s1)) +
 * abs(l/2 - (l/2 + s2) where <l> is the vector size and <s1> and <s2>
 * are the 1st positions of each slice) corresponding to the
 * absolute difference of position of the aligned slice to the middle
 * of the complete vectors.
 * In case where a correlation between both vectors cannot be computed
 * (for instance one vector has a variance of 0), the distance is set
 * to 1 (corresponding to a correlation of 0).
 * \param x the first vector.
 * \param y the second vector.
 * \param shift the number of allowed shift states (different alignments).
 * \param shift1 used to report the shift state used for <x>
 * in the optimal alignment (corresponds to the starting
 * position of the slice : x[shift1]).
 * \param shift2 used to report the shift state used for <y>
 * in the optimal alignment (corresponds to the starting
 * position of the slice : y[shift2]).
 * \return the distance between both vectors. If the correlation
 * between both vectors could not be computed (a vector has a
 * sd of 0), then the distance is set to 1 (the correlation
 * is considered to be 0).
 */
double cor_distance(const std::vector<double>& x,
                    const std::vector<double>& y,
                    size_t shift,
                    size_t& shift1,
                    size_t& shift2) ;


/*!
 * \brief Computes the correlation distance between two vectors
 * x and y of equal lengths, allowing up to <shift> number of
 * shifting states. This function performs all possible
 * comparisons (all possible alignments) between each slice of
 * x and y of length (l-shift+1) where <l> is the vector length
 * and <shift> the value of the corresponding parameter. Additionally,
 * one extra alignment with the reversed slice of y is performed.
 * If several alignments results in the same optimal
 * distance, the most central alignment is preferred (each alignment
 * is assigned a score using : abs(l/2 - (l/2 + s1)) +
 * abs(l/2 - (l/2 + s2) where <l> is the vector size and <s1> and <s2>
 * are the 1st positions of each slice) corresponding to the
 * absolute difference of position of the aligned slice to the middle
 * of the complete vectors.
 * In case where a correlation between both vectors cannot be computed
 * (for instance one vector has a variance of 0), the distance is set
 * to 1 (corresponding to a correlation of 0).
 * \param x the first vector.
 * \param y the second vector.
 * \param shift the number of allowed shift states (different alignments).
 * \param shift1 used to report the shift state used for <x>
 * in the optimal alignment (corresponds to the starting
 * position of the slice : x[shift1]).
 * \param shift2 used to report the shift state used for <y>
 * in the optimal alignment (corresponds to the starting
 * position of the slice : y[shift2]).
 * \param flip used to report whether the slice of <y> was reversed
 * in the optimal alignment. If true, it simply means that y[shift2] was
 * the last element of the slice.
 * \return the distance between both vectors. If the correlation
 * between both vectors could not be computed (a vector has a
 * sd of 0), then the distance is set to 1 (the correlation
 * is considered to be 0).
*/
double cor_distance(const std::vector<double>& x,
                    const std::vector<double>& y,
                    size_t shift,
                    size_t& shift1,
                    size_t& shift2,
                    bool& flip) ;


/*!
 * \brief Does exactly the same as cor_distance() excepted that this
 * implementation uses a dynamic programming approach to speed up the
 * computations.
 * \param x the first vector.
 * \param y the second vector.
 * \param shift the number of allowed shift states (different alignments).
 * \param shift1 used to report the shift state used for <x>
 * in the optimal alignment (corresponds to the starting
 * \param shift2 used to report the shift state used for <y>
 * in the optimal alignment (corresponds to the starting
 * position of the slice : y[shift2]).
 * \return the distance between both vectors. If the correlation
 * between both vectors could not be computed (a vector has a
 * sd of 0), then the distance is set to 1 (the correlation
 * is considered to be 0).
 */
double cor_distance_fast(const std::vector<double>& x,
                         const std::vector<double>& y,
                         size_t shift,
                         size_t& shift1,
                         size_t& shift2) ;


/*!
 * \brief Does exactly the same as cor_distance() excepted that this
 * implementation uses a dynamic programming approach to speed up the
 * computations.
 * \param x the first vector.
 * \param y the second vector.
 * \param shift the number of allowed shift states (different alignments).
 * \param shift1 used to report the shift state used for <x>
 * in the optimal alignment (corresponds to the starting
 * \param shift2 used to report the shift state used for <y>
 * in the optimal alignment (corresponds to the starting
 * position of the slice : y[shift2]).
 * \param flip used to report whether the slice of <y> was reversed
 * in the optimal alignment. If true, it simply means that y[shift2] was
 * the last element of the slice.
 * \return the distance between both vectors. If the correlation
 * between both vectors could not be computed (a vector has a
 * sd of 0), then the distance is set to 1 (the correlation
 * is considered to be 0).
 */
double cor_distance_fast(const std::vector<double>& x,
                         const std::vector<double>& y,
                         size_t shift,
                         size_t& shift1,
                         size_t& shift2,
                         bool& flip) ;

/*!
 * \brief This function corrects (modify) a correlation value in three cases :
 * case 1) rounding to 0, a value is 0 +/-2xepsilon machine because of rounding errors,
 * case 2) rounding to 1, a value is 1 +/-2xepsilon machine because of rounding errors,
 * case 3) turns to Constants::corr_fail because it is a NaN value.
 * \param x the correlation to check for correction.
 * \return whether the value has been modified.
 */
bool correctCor(double& x) ;

#endif //DISTANCES_HPP
