#ifndef STATISTCS_HPP
#define STATISTICS_HPP

#include <vector>
#include <numeric>
#include <cmath>     // pow(), sqrt()
#include <assert.h>
#include <iostream>
#include <Utility/Vector_utility.hpp>


/*!
 * \brief Computes the mean of a given vector. If from and to are
 * different from -1, computes the mean using the slice v[from,to).
 * \param v a vector of interest.
 * \param from the starting index, if different from -1.
 * \param to the ending index (not included), if different from -1.
 * \return the mean.
 */
template<class T>
inline double mean(const std::vector<T>& v, int from=-1, int to=-1) ;


/*!
 * \brief Computes the weighted mean of a vector of measures <x> given their
 * probability <p>. The sum of <p> is expected to be 1, if not it will
 * be normalized.
 * \param x a vector of measures.
 * \param p a vector of probability associated to each element of <x>.
 * \param from the starting index, if different from -1.
 * \param to the ending index (not included), if different from -1.
 * \return the mean.
 */
template<class T>
inline double mean(const std::vector<T>& x, const std::vector<double>& p) ;


/*! Computes the unbiased standard deviation of a given vector. If from and
 * to are different from -1, computes the standard deviation using the
 * slice v[from,to).
 * \param v a vector of interest.
 * \param from the starting index, if different from -1.
 * \param to the ending index (not included), if different from -1.
 * \return the unbiased standard deviation.
 */
template<class T>
inline double sd(const std::vector<T>& v, int from=-1, int to=-1) ;

/*!
 * \brief Computes the unbiased standard deviation of a vector of
 * measures <x> given their probability <p>. The sum of <p> is
 * expected to be 1, if not it will be normalized.
 * \param x a vector of measures.
 * \param p a vector of probability associated to each element of <x>.
 * \param from the starting index, if different from -1.
 * \param to the ending index (not included), if different from -1.
 * \return the mean.
 */
template<class T>
inline double sd(const std::vector<T>& x, const std::vector<double>& p) ;

/*!
 * \brief Computes the pearson correlation coefficient of two given vectors
 * using the slices v1[from1-->to1) and v2[from2-->to2). If from1, to1, from2,
 * to2 are equal to -1, they are not accounted and the first and last elements
 * considered are the beginning and the end of the vectors.
 * \param v1 the first vector.
 * \param v2 the second vector.
 * \param from1 the index of the first value to use in the first vector.
 * \param to1 the index after the last value to use in the fist vector.
 * \param from2 the index of the first value to use in the second vector.
 * \param to2 the index after the last value to use in the second vector.
 * \return the correlation coefficient
 */
template<class T>
inline double cor_pearson(const std::vector<T>& v1,
                          const std::vector<T>& v2,
                          int from1=-1, int to1=-1,
                          int from2=-1, int to2=-1) ;

/*!
 * \brief Computes the pearson correlation coefficient for two given vectors
 * using the slices v1[from1-->to1) and v2(to2-->from2] (the second vector is
 * read backward, from the end).
 * If from1, to1, from2, to2 are equal to -1, they are not accounted and the first
 * and last elements considered are the beginning and the end of the vectors (for
 * the last element considered in v2 to be the 1st of the vector (the 0th), to2
 * should be set to -1, which is ignored. However this enables the default behaviour
 * using the 0th element as the last one).
 * \param v1 the first vector.
 * \param v2 the second vector.
 * \param from1 the index of the first value to use in the first vector.
 * \param to1 the index after the last value to use in the fist vector.
 * \param from2 the index of the first value to use in the second vector.
 * \param to2 the index after the last value to use in the second vector.
 * \return the correlation coefficient
 */
template<class T>
inline double cor_pearson_rev(const std::vector<T>& v1,
                              const std::vector<T>& v2,
                              int from1=-1, int to1=-1,
                              int from2=-1, int to2=-1) ;

template<class T>
inline double mean(const std::vector<T>& v, int from=-1, int to=-1)
{
    typename std::vector<T>::const_iterator from_it = v.begin() ;
    typename std::vector<T>::const_iterator to_it   = v.end() ;

    int size = v.size() ;
    if(from >= 0 and from <= size)
    {   from_it += from ; }
    if(to >= 0 and to <= size)
    {   to_it = v.begin() + to ; }

    // 0 <= from < to <= v.size()
    // if from == to -> return nan
    // if from > to  -> segfault
    assert(from_it >= v.begin()) ;
    assert(from_it < to_it) ;
    assert(to_it <= v.end()) ;

    double   n = to_it - from_it ;
    double sum = std::accumulate(from_it, to_it, 0.) ;
    return sum / n ;
}


template<class T>
inline double mean(const std::vector<T>& x, const std::vector<double>& p)
{
    assert(x.size() == p.size()) ;

    double mean = 0. ;
    double total = 0. ;
    for(auto i : p)
    {   total += i ; }
    for(size_t i=0; i<x.size(); i++)
    {   mean += x[i] * (p[i]/total) ; }
    return mean ;
}


template<class T>
inline double sd(const std::vector<T>& v, int from=-1, int to=-1)
{
    int size = v.size() ;

    // can't compute the sd of a single value...
    if(size <= 1 or (to - from == 1))
    {   return std::nan("") ; }

    if(from == -1)
    {   from = 0 ; }
    if(to == -1)
    { to = size ; }

    // 0 <= from < to <= v.size()
    // if from == to -> no sense, no value : return nan
    // if from  > to -> no sense, backward : does a segfault
    assert(from >= 0) ;
    assert(from < to) ;
    assert(to <= size) ;

    double m   = mean(v,from, to) ;
    double sum = 0. ;
    double n   = to - from ;

    for(size_t i=from; i<to; i++)
    {   sum += pow((v[i] - m),2) ; }
    return sqrt(sum / (n-1)) ;
}


template<class T>
inline double sd(const std::vector<T>& x, const std::vector<double>& p)
{
    assert(x.size() == p.size()) ;

    double v1    = 0. ;
    double v2    = 0. ;
    double m     = mean(x,p) ;
    double sd    = 0. ;
    double total = 0. ;
    for(auto i : p)
    {   total += i ; }
    for(size_t i=0; i<x.size(); i++)
    {   double p_norm = p[i] / total ;
        sd += pow(x[i] - m, 2) * (p_norm) ;
        v1 += p_norm ;
        v2 += pow(p_norm,2) ;
    }
    return sqrt(sd / (v1 - (v2/v1))) ;
}


template<class T>
inline double cor_pearson(const std::vector<T>& v1,
                          const std::vector<T>& v2,
                          int from1=-1, int to1=-1,
                          int from2=-1, int to2=-1)
{
    if(from1 == -1)
    {   from1 = 0 ; }
    if(to1 == -1)
    {   to1 = v1.size() ; }
    if(from2 == -1)
    {   from2 = 0 ; }
    if(to2 == -1)
    {   to2 = v2.size() ; }

    // 0 <= from < to < v.size()
    // if from == to -> no sense, no value : return nan
    // if from  > to -> no sense, backward : does a segfault
    assert(from1 >= 0) ;
    assert(from1 < to1) ;
    assert(to1 <= v1.size()) ;
    assert(from2 >= 0) ;
    assert(from2 < to2) ;
    assert(to2 <= v2.size()) ;

    // the terms of the formula
    double sum_v1v2 = 0. ;
    double sum_v1   = 0. ;
    double sum_v2   = 0. ;
    double sum_v1_2 = 0. ;
    double sum_v2_2 = 0. ;
    // the effective number of values considered
    double n = to1 - from1 ;

    for(size_t i1=from1, i2=from2; i1<to1; i1++,i2++)
    {   sum_v1v2 += v1[i1]*v2[i2] ;
        sum_v1   += v1[i1] ;
        sum_v2   += v2[i2] ;
        sum_v1_2 += pow(v1[i1],2) ;
        sum_v2_2 += pow(v2[i2],2) ;

    }

    return (n*sum_v1v2 - (sum_v1*sum_v2)) /
            (sqrt(n*sum_v1_2 - pow(sum_v1, 2)) * sqrt(n*sum_v2_2 - pow(sum_v2, 2)) ) ;
}

template<class T>
inline double cor_pearson_rev(const std::vector<T>& v1,
                              const std::vector<T>& v2,
                              int from1=-1, int to1=-1,
                              int from2=-1, int to2=-1)
{
    if(from1 == -1)
    {   from1 = 0 ; }
    if(to1 == -1)
    {   to1 = v1.size() ; }
    if(from2 == -1)
    {   from2 = v2.size() - 1 ; }
    if(to2 == -1)
    {   to2 = 0 ; }

    // beware -> for v1 [from1, to1)
    //           for v2 [to2, from2]
    // use from1 and to1 for the loop conditions

    // 0 <= from < to    < v1.size()
    // 0 <= to2  < from2 < v2.size() because for v2 [to2, from2]
    // if from == to -> no sense, no value : return nan
    // if from  > to -> no sense, backward : does a segfault
    assert(from1 >= 0) ;
    assert(from1 < to1) ;
    assert(to1   <= v1.size()) ;
    assert(to2   >= 0) ;
    assert(to2   < from2) ;
    assert(from2 <  v2.size()) ;

    // the terms of the formula
    double sum_v1v2 = 0. ;
    double sum_v1   = 0. ;
    double sum_v2   = 0. ;
    double sum_v1_2 = 0. ;
    double sum_v2_2 = 0. ;
    // the effective number of values considered
    double n = to1 - from1 ;

    for(int i1=from1, i2=from2 ; i1<to1; i1++,i2--)
    {   sum_v1v2 += v1[i1]*v2[i2] ;
        sum_v1   += v1[i1] ;
        sum_v2   += v2[i2] ;
        sum_v1_2 += pow(v1[i1],2) ;
        sum_v2_2 += pow(v2[i2],2) ;
    }

    return (n*sum_v1v2 - (sum_v1*sum_v2)) /
            (sqrt(n*sum_v1_2 - pow(sum_v1, 2)) * sqrt(n*sum_v2_2 - pow(sum_v2, 2)) ) ;
}




#endif // STATISTICS_HPP
