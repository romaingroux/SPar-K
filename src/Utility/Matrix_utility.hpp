#ifndef MATRIX_UTILITY_HPP
#define MATRIX_UTILITY_HPP

#include <stdexcept>                 // std::invalid_argument
#include <cmath>
#include <Matrix/Matrix2D.hpp>

/*!
 * \brief Smooth the outliers (row-wise) of a given matrix - any value
 * which is bigger/smaller than the row mean plus/minus three row
 * standard deviations - and replace these values by the row mean +/-
 * three row standard deviations (if the value is a bottom or top outlier
 * respectively). If several outliers follow each other, then because the
 * procedure is done in a left to right manner, some values may still be
 * outliers after.
 * \param m a matrix of interest.
 * \throw std::invalid_argument if the size of the vector is < 2.
 */
template<class T>
void smooth_outliers(Matrix2D<T>& m) throw (std::invalid_argument)
{
    if(m.get_ncol() < 3)
    {   throw std::invalid_argument("Cannot replace outliers in a matrix with ncol < 3!") ; }

    // replace outliers by threshold
    for(size_t i=0; i<m.get_nrow(); i++)
    {   // compute mean and sd on this row
        T mean = 0. ;
        T n    = 0. ;
        for(size_t j=0; j<m.get_ncol(); j++)
        {   mean += m(i,j) ;
            n    += 1. ;
        }
        mean /= n ;
        T sum = 0. ;
        for(size_t j=0; j<m.get_ncol(); j++)
        {   sum += pow(m(i,j) - mean, 2) ; }
        T sd = sqrt(sum / (n-1.)) ;
        T upper = mean + 3*sd ; // upper threshold
        T lower = mean - 3*sd ; // lower threshold

        // replace outliers on this row
        for(size_t j=0; j<m.get_ncol(); j++)
        {
            if(m(i,j) > upper or
               m(i,j) < lower)
            {   // only neighbours on the right
                if(j == 0)
                {   m(i,j) = (m(i,j+1) + m(i,j+2)) / 2. ; }
                // only neighbours on the left
                else if(j == m.get_ncol()-1)
                {   m(i,j) = (m(i,j-1) + m(i,j-2)) / 2. ; }
                // neighbours on both sides
                else
                {   m(i,j) = (m(i,j-1) + m(i,j+1)) / 2. ; }
            }
        }
    }
}

#endif // MATRIX_UTILITY_HPP
