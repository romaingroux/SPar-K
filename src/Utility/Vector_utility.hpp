#ifndef VECTOR_UTILITY_HPP
#define VECTOR_UTILITY_HPP

#include <vector>
#include <iostream>
#include <cassert>


/*!
 * \brief Overload the << operator for vectors.
 * \param out an output stream.
 * \param v a vector of interest.
 * \param sep a character to separate the values.
 * \return a reference to the output stream.
 */
template<class T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v)
{   for(auto i:v)
    {	out << i << " " ; }
    return out ;
}


/*!
 * \brief A method to print a vector. Unlike the << operator, this
 * function allows to specify the separator character.
 * \param out an output stream.
 * \param v a vector of interest.
 * \param sep a character to separate the values.
 */
template<class T>
void print_vector(std::ostream& out, const std::vector<T>& v, char sep=' ')
{   for(auto i:v)
    {	out << i << sep ; }
}


/*!
 * \brief Given a vector and an offset, this function pads (shifts)
 * the vector content by adding values (0 by default) on one side of
 * the vector but keeping the size of the vector constant.
 * \param v the vector.
 * \param offset the offset, negative means padding to the left, positive
 * means padding to the right.
 * \param value the value to use for the padding.
 */
template<class T>
void pad_vector(std::vector<T>& v, int offset, T value=0.)
{
    int from = 0 ;
    int to   = 0 ;
    size_t len  = v.size() ;

    assert((offset < v.size()) or (offset > v.size())) ;

    // padding to the left
    if(offset < 0)
    {   from = std::abs(offset) ;
        to   = static_cast<int>(len-1) ;
        for(int i=from; i<=to; i++)
        {   v[i+offset] = v[i] ; }
        // fill freed slots (extreme right) with values
        for(size_t i=to+offset+1; i<len; i++)
        {   v[i] = value ; }
    }
    // padding to the right
    else if(offset > 0)
    {   from = 0 ;
        to   = len - 1 - offset ;
        for(int i=to; i>=from; i--)
        {   v[i+offset] = v[i] ; }
        // fill freed slots (extreme left) with values
        for(size_t i=from; i<offset; i++)
        {   v[i] = value ; }
    }
    // no padding
    else
    {}
}


#endif // VECTOR_UTILITY_HPP
