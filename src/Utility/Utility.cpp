#include "Utility.hpp"

#include <limits>  // std::numeric_limits<double>::epsilon()
#include <cmath>   // std::abs()

std::vector<int> seq(int from, int to, int by)
{   std::vector<int> v ;
    for(int current=from; current<=to; current += by)
    {   v.push_back(current) ; }
    return v ;
}


std::vector<size_t> find_max(const Matrix2D<double>& m)
{   size_t x = 0, y = 0 ;
    // init to the lowest possible value
    double current_max = std::numeric_limits<double> ::min() ;
    for(size_t i=0; i<m.get_nrow(); i++)
    {   for(size_t j=0; j<m.get_ncol(); j++)
        {   if(m(i,j) > current_max)
            {   current_max = m(i,j) ;
                x = i ;
                y = j ;
            }
        }
    }
    return std::vector<size_t> {x, y} ;
}


std::vector<size_t> find_min(const Matrix2D<double>& m)
{   size_t x = 0, y = 0 ;
    // init to the highest possible value
    double current_min = std::numeric_limits<double> ::max() ;
    for(size_t i=0; i<m.get_nrow(); i++)
    {   for(size_t j=0; j<m.get_ncol(); j++)
        {   if(m(i,j) < current_min)
            {   current_min = m(i,j) ;
                x = i ;
                y = j ;
            }
        }
    }
    return std::vector<size_t> {x, y} ;
}


bool roundToZero(double& x)
{   if((x >= 2.*-std::numeric_limits<double>::epsilon()) and (x <= 2.*std::numeric_limits<double>::epsilon()))
    {   x = 0. ;
        return true ;
    }
    return false ;
}

bool roundToOne(double& x)
{   if((x >= 1-(2.*std::numeric_limits<double>::epsilon())) and
       (x <= 1+(2.*std::numeric_limits<double>::epsilon())))
    {   x = 1. ;
        return true ;
    }
    return false ;
}


bool isEqual(double x, double y, double epsilon)
{   return std::abs(x - y) < epsilon ; }
