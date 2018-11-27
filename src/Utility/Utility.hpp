#ifndef UTILITY_HPP
#define UTILITY_HPP

#include <vector>
#include <ostream>
#include <assert.h>
#include <algorithm>

#include "Matrix/Matrix2D.hpp"
#include "Constants.hpp"
#include "Vector_utility.hpp"


/*!
 * \brief Creates a sequence of number starting with <from>,
 * stopping at <to> with increments of step <by>.
 * \param from the start value.
 * \param to the upper most limit.
 * \param by the increment.
 * \return a vector with the sequence.
 */
std::vector<int> seq(int from, int to, int by) ;


/*!
 * \brief Finds and return the coordinates of the maximum in the given
 * matrix.
 * \param m a matrix of interest.
 * \return a vector containing the x/y coordinates.
 */
std::vector<size_t> find_max(const Matrix2D<double>& m) ;


/*!
 * \brief Finds and return the coordinates of the minimum in the given
 * matrix.
 * \param m a matrix of interest.
 * \return a vector containing the x/y coordinates.
 */
std::vector<size_t> find_min(const Matrix2D<double>& m) ;


/*! Checks whether a value is NaN.
 * According to IEEE754 NaN != NaN.
 * This should be portable...
 * \arg x the value.
 * \return whether the value is a boolean.
 */
template<class T>
bool isNaN(T x)
{   return x != x ; }


/*!
 * \brief Checks if a double should be rounded to 0, that is if
 * it lies within [-2*espilon machine, +2*epsilon machine] range.
 * If so, x is rounded to 0.
 * +/-2 because some given values may have been obtained using
 * several other values subjected to rounding errors already.
 * \param x the number to round.
 * \return whether the given value was modified.
 */
bool roundToZero(double& x) ;


/*!
 * \brief Checks if a double should be rounded to 1, that is if
 * it lies within [1-2*espilon machine, 1+2*epsilon machine] range.
 * If so, x is rounded to 1.
 * +/-2 because some given values may have been obtained using
 * several other values subjected to rounding errors already.
 * \param x the number to round.
 * \return whether the given value was modified.
 */
bool roundToOne(double& x) ;

/*!
 * \brief Performs an equality check between two doubles using a
 * given precision parameter.
 * \param x the first double.
 * \param y the second double.
 * \param epsilon the precision parameters
 * \return returns true if |x-y| < epsilon.
 */
bool isEqual(double x, double y, double epsilon) ;


#endif // UTILITY_HPP
