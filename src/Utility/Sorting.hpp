#ifndef SORTING_HPP
#define SORTING_HPP

#include <vector>
#include <algorithm>

/*! \brief Given a vector of data, it returns a vector of index
 * corresponding to the index of each element after data sorting
 * in ascending (by default) order.
 * \param data a vector of interest.
 * \param ascending specifies that the index should be returned in
 * ascending order
 * \return the index after sorting.
 */
template<class T>
std::vector<std::size_t> order(const std::vector<T>& data, bool ascending=true)
{
    std::vector<std::size_t> index(data.size(), 0) ;

    for(std::size_t i=0; i<index.size(); i++)
    {   index[i] = i ; }

    sort(index.begin(), index.end(),
        [&](const int& a, const int& b)
        {   return (data[a] < data[b]) ; }
    ) ;

    if(not ascending)
    {   std::reverse(index.begin(), index.end()) ; }

    return index ;
}



#endif // SORTING_HPP
