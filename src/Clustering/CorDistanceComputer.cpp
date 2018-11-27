#include "CorDistanceComputer.hpp"

#include <vector>
#include <stdexcept>                 // std::invalid_argument
#include <cmath>
#include <limits>                    // numeric_limits<double>::epsilon(), numeric_limits<size_t>::max(), numeric_limits<double>::max()
#include <Utility/Utility.hpp>       // isNaN()
#include <Utility/Constants.hpp>     // Constants
#include <Statistics/Statistics.hpp> // mean(), sd()


CorDistanceComputer::CorDistanceComputer(size_t shift) throw (std::invalid_argument)
    : _shift(0)
{   this->reallocate_buffers(shift) ; }

CorDistanceComputer::~CorDistanceComputer()
{}

double CorDistanceComputer::compute_distance(const std::vector<double>& x,
                                             const std::vector<double>& y,
                                             size_t shift,
                                             size_t& shift_x,
                                             size_t& shift_y) throw (std::invalid_argument)
{   // computes without flipping
    bool flip_y = false ;
    try
    {   return this->compute(x,y,shift,shift_x,shift_y,flip_y) ; }
    catch(std::invalid_argument e)
    {   throw e ; }
}

double CorDistanceComputer::compute_distance(const std::vector<double>& x,
                                             const std::vector<double>& y,
                                             size_t shift,
                                             size_t& shift_x,
                                             size_t& shift_y,
                                             bool&   flip_y) throw (std::invalid_argument)
{   double d ;
    // computes with flipping
    bool flip = true ;
    try
    {   d = this->compute(x,y,shift,shift_x,shift_y,flip) ;
        flip_y = flip ;
        return d ;
    }
    catch(std::invalid_argument e)
    {   throw e ; }
}


double CorDistanceComputer::compute(const std::vector<double>& x,
                                    const std::vector<double>& y,
                                    size_t shift,
                                    size_t& shift_x,
                                    size_t& shift_y,
                                    bool& flip) throw (std::invalid_argument)
{
    // reallocate buffers if neeeded
    bool realloc = this->reallocate_buffers(shift) ;

    // check that both vectors have the same length
    if(x.size() != y.size())
    {   throw std::invalid_argument("vectors of different length!") ; }
    // the length of the vectors
    size_t l   = x.size() ;
    int l_half = l / 2 ;
    // the lenght of a vector slice taking the shift into account
    size_t n    = l - shift + 1 ;
    int n_half  = n / 2;
    if(not n >= 1)
    {   std::invalid_argument("shift too high/vectors too short!") ; }

    // some iterators on x and y -> always x[from_x,to_x) and y[from_y,to_y)
    size_t from_x, to_x ;
    size_t from_y, to_y, from_y_rev, to_y_rev ;

    // to track the minimum distance
    double min_distance  = 3. ;                      // the current minimum distance value
    std::vector<std::vector<size_t>> min_distances ; // the coordinates of all distance equal to the current min distance

    // initialize the sums of x and y elements for each shift
    if(not realloc)
    {   this->_sum_x[0] = this->_sum_x2[0] = this->_sum_y[0] = this->_sum_y2[0] = 0 + Constants::delta ;
        this->_sum_x[1] = this->_sum_x2[1] = this->_sum_y[1] = this->_sum_y2[1] = 0 + Constants::delta ;
    }
    for(size_t i=0; i<n; i++)
    {   // add Constants::delta to avoid zeroes and division by zero
        // for instance in cases where the x or the y vectors contain
        // only 0's (which case is often encountered)
        this->_sum_x[0]  += x[i]        + Constants::delta ;
        this->_sum_x2[0] += pow(x[i],2) + Constants::delta ;
        this->_sum_y[0]  += y[i]        + Constants::delta ;
        this->_sum_y2[0] += pow(y[i],2) + Constants::delta ;
    }
    for(size_t s=1, i=n; s<shift; s++, i++)
    {   this->_sum_x[s]  += this->_sum_x[s-1]  - x[i-n]        + x[i] ;
        this->_sum_x2[s] += this->_sum_x2[s-1] - pow(x[i-n],2) + pow(x[i],2) ;
        this->_sum_y[s]  += this->_sum_y[s-1]  - y[i-n]        + y[i] ;
        this->_sum_y2[s] += this->_sum_y2[s-1] - pow(y[i-n],2) + pow(y[i],2) ;
        if((not realloc) and (s+1 < shift))
        {   this->_sum_x[s+1] = this->_sum_x2[s+1] = this->_sum_y[s+1] = this->_sum_y2[s+1] = 0 ; }
    }

    // fill 1st row of dist -> compute the 1-cor with x having a shift of 0
    from_x = 0 ;
    to_x   = from_x + n ; // x[from_x,to_x)
    for(size_t i=0; i<shift; i++) // i is the shift of y
    {   from_y = i ;
        to_y   = from_y + n ; // y[from_y,to_y)

        // forward orientation -> x fw vs y fw
        {
            this->_sum_xy[from_x][from_y][Constants::FORWARD] = 0. ;
            for(size_t j=0; j<n; j++)
            {   this->_sum_xy[from_x][from_y][Constants::FORWARD] += x[from_x+j]*y[from_y+j] ; }


            double cor = (n*this->_sum_xy[from_x][from_y][Constants::FORWARD] - this->_sum_x[from_x]*this->_sum_y[from_y]) /
                         (sqrt(n*this->_sum_x2[from_x] - pow(this->_sum_x[from_x], 2))*sqrt(n*this->_sum_y2[from_y] - pow(this->_sum_y[from_y], 2))) ;
            CorDistanceComputer::correct_cor(cor) ;
            double d = 1. - cor ; roundToZero(d) ;
            // the distance
            this->_distances[from_x][from_y][Constants::FORWARD] = d ;
            // the distance score
            this->_scores[from_x][from_y] = abs(l_half - (n_half+static_cast<int>(from_x))) +  // the absolute distance of the middle of v1 slice to middle of v1
                                     abs(l_half - (n_half+static_cast<int>(from_y))) ;         // the absolute distance of the middle of v2 slice to middle of v2
            // keep track of the min distance
            if(d < min_distance)
            {   min_distance = d ;
                min_distances = std::vector<std::vector<size_t>>() ;
                min_distances.push_back(std::vector<size_t>({from_x, from_y, Constants::FORWARD})) ;
            }
            else if(d == min_distance)
            {   min_distances.push_back(std::vector<size_t>({from_x, from_y, Constants::FORWARD})) ; }

        }
        // reverse orientation -> x fw vs y rev
        if(flip)
        {   this->_sum_xy[from_x][from_y][Constants::REVERSE] = 0. ;
            for(size_t j=0; j<n; j++)
            {   this->_sum_xy[from_x][from_y][Constants::REVERSE] += x[from_x+j]*y[to_y-j-1] ; }

            double cor = (n*this->_sum_xy[from_x][from_y][Constants::REVERSE] - this->_sum_x[from_x]*this->_sum_y[from_y]) /
                         (sqrt(n*this->_sum_x2[from_x] - pow(this->_sum_x[from_x], 2))*sqrt(n*this->_sum_y2[from_y] - pow(this->_sum_y[from_y], 2))) ;
            CorDistanceComputer::correct_cor(cor) ;
            double d = 1. - cor ; roundToZero(d) ;
            // distance
            this->_distances[from_x][from_y][Constants::REVERSE] = d ;
            // the distance score
            // no need to compute, same as forward
            // keep track of the min distance
            if(d < min_distance)
            {   min_distance = d ;
                min_distances = std::vector<std::vector<size_t>>() ;
                min_distances.push_back(std::vector<size_t>({from_x, from_y, Constants::REVERSE})) ;
            }
            else if(d == min_distance)
            {   min_distances.push_back(std::vector<size_t>({from_x, from_y, Constants::REVERSE})) ;}
        }
    }

    // fill 1st column of dist -> compute the 1-cor with y having a shift of 0 in fw and a shift of <shift-1> in rev
    from_y     = 0 ;
    to_y       = from_y + n ;         // y[from_y,to_y)
    from_y_rev = shift - from_y - 1 ; // correspond to the slice with an equal shift but starting at the other edge of the vector
    to_y_rev   = from_y_rev + n ;     // y[from_y_rev, to_y_rev)
    for(size_t i=0; i<shift; i++) // i is the shift of x
    {   from_x = i ;
        to_x   = from_x + n ; // x[from_x,to_x)

        // forward orientation -> x fw vs y fw
        {   this->_sum_xy[from_x][from_y][Constants::FORWARD] = 0. ;
            for(size_t j=0; j<n; j++)
            {   this->_sum_xy[from_x][from_y][Constants::FORWARD] += x[from_x+j]*y[from_y+j] ; }

            double cor = (n*this->_sum_xy[from_x][from_y][Constants::FORWARD] - this->_sum_x[from_x]*this->_sum_y[from_y]) /
                         (sqrt(n*this->_sum_x2[from_x] - pow(this->_sum_x[from_x], 2))*sqrt(n*this->_sum_y2[from_y] - pow(this->_sum_y[from_y], 2))) ;
            CorDistanceComputer::correct_cor(cor) ;
            double d = 1. - cor ; roundToZero(d) ;
            // distance
            this->_distances[from_x][from_y][Constants::FORWARD] = d ;
            // the distance score
            this->_scores[from_x][from_y] = abs(l_half - (n_half+static_cast<int>(from_x))) +  // the absolute distance of the middle of v1 slice to middle of v1
                                     abs(l_half - (n_half+static_cast<int>(from_y))) ;         // the absolute distance of the middle of v2 slice to middle of v2
            // keep track of the min distance
            if(d < min_distance)
            {   min_distance = d ;
                min_distances = std::vector<std::vector<size_t>>() ;
                min_distances.push_back(std::vector<size_t>({from_x, from_y, Constants::FORWARD})) ;
            }
            else if(d == min_distance)
            {   min_distances.push_back(std::vector<size_t>({from_x, from_y, Constants::FORWARD})) ;}

        }
        // reverse orientation -> x fw vs y rev
        if(flip)
        {   this->_sum_xy[from_x][shift-1][Constants::REVERSE] = 0. ;
            for(size_t j=0; j<n; j++)
            {   this->_sum_xy[from_x][shift-1][Constants::REVERSE] += x[from_x+j]*y[to_y_rev-1-j] ; }

            double cor = (n*this->_sum_xy[from_x][shift-1][Constants::REVERSE] -this->_sum_x[from_x]*this->_sum_y[shift-1]) /
                         (sqrt(n*this->_sum_x2[from_x] - pow(this->_sum_x[from_x], 2))*sqrt(n*this->_sum_y2[shift-1] - pow(this->_sum_y[shift-1], 2))) ;
            CorDistanceComputer::correct_cor(cor) ;
            double d = 1. - cor ; roundToZero(d) ;
            // distance
            this->_distances[from_x][shift-1][Constants::REVERSE] = d ;
            // the distance score
            // no need to compute, same as forward
            // keep track of the min distance
            if(d < min_distance)
            {   min_distance = d ;
                min_distances = std::vector<std::vector<size_t>>() ;
                min_distances.push_back(std::vector<size_t>({from_x, shift-1, Constants::REVERSE})) ;
            }
            else if(d == min_distance)
            {   min_distances.push_back(std::vector<size_t>({from_x, shift-1, Constants::REVERSE})) ; }
        }
    }

    // compute the remaining distances
    for(size_t i=1; i<shift; i++)
    {   from_x = i ;
        to_x   = from_x + n ; // x[from_x,to_x)
        for(size_t j=1; j<shift; j++)
        {   from_y = j ;
            to_y   = from_y + n ; // y[from_y,to_y)

            // forward orientation -> x fw vs y fw
            {   // compute the current some of pairwise product from a previous one
                this->_sum_xy[from_x][from_y][Constants::FORWARD] = this->_sum_xy[from_x-1][from_y-1][Constants::FORWARD] - x[from_x-1]*y[from_y-1] + x[to_x-1]*y[to_y-1] ;
                double cor                                 = (n*this->_sum_xy[from_x][from_y][Constants::FORWARD] - this->_sum_x[from_x]*this->_sum_y[from_y]) /
                                                             (sqrt(n*this->_sum_x2[from_x] - pow(this->_sum_x[from_x], 2))*sqrt(n*this->_sum_y2[from_y] - pow(this->_sum_y[from_y], 2))) ;
                CorDistanceComputer::correct_cor(cor) ;
                double d = 1. - cor ; roundToZero(d) ;
                // distance
                this->_distances[from_x][from_y][Constants::FORWARD] = d ;
                // the distance score
                this->_scores[from_x][from_y] = abs(l_half - (n_half+static_cast<int>(from_x))) +  // the absolute distance of the middle of v1 slice to middle of v1
                                         abs(l_half - (n_half+static_cast<int>(from_y))) ;         // the absolute distance of the middle of v2 slice to middle of v2
                // keep track of the min distance
                if(d < min_distance)
                {   min_distance = d ;
                    min_distances = std::vector<std::vector<size_t>>() ;
                    min_distances.push_back(std::vector<size_t>({from_x, from_y, Constants::FORWARD})) ;

                }
                else if(d == min_distance)
                {   min_distances.push_back(std::vector<size_t>({from_x, from_y, Constants::FORWARD})) ; }
            }
            // reverse orientation -> x fw vs y rev
            // in forward order the values for <shift+1> can be obtained from the values for <shift>, it it the same
            // in reverse, except that the shift value is assigned at the end of the vector instead of the beginning
            // | shift -->
            // |from        to|
            //  1  2  3  4  5  6  7  8  9  10  11  12
            //                    |from_rev   to_rev|
            //                 <--|shift_rev
            // So to compute the values at <shift+1> for reverse (forward numbering), one can compute them from <shift+2>
            // (in forward numbering) (so shift in forward numbering decreases). However, if considering the shift numbering
            // in reverse (cf above), then the values at <shift+1> can be computed from those at <shift>.
            if(flip)
            {   size_t j_rev = shift - j - 1 ;    // gives the complementary shift state -> if shift=3 : 0<->2, 1<->1, 2<->0
                from_y_rev = j_rev  ;             // the 1st element of the slice if considering <shift> assigned at the end
                to_y_rev   = from_y_rev + n     ; // y[from_rev,to_rev)

                // compute the current some of pairwise product from a previous one
                this->_sum_xy[from_x][from_y_rev][Constants::REVERSE] = this->_sum_xy[from_x-1][from_y_rev+1][Constants::REVERSE] - x[from_x-1]*y[to_y_rev] + x[to_x-1]*y[from_y_rev] ;
                double cor  = (n*this->_sum_xy[from_x][from_y_rev][Constants::REVERSE] - this->_sum_x[from_x]*this->_sum_y[from_y_rev]) /
                              (sqrt(n*this->_sum_x2[from_x] - pow(this->_sum_x[from_x], 2))*sqrt(n*this->_sum_y2[from_y_rev] - pow(this->_sum_y[from_y_rev], 2))) ;
                CorDistanceComputer::correct_cor(cor) ;
                double d = 1. - cor ; roundToZero(d) ;
                // distance
                this->_distances[from_x][from_y_rev][Constants::REVERSE] = d ;
                // the distance score
                // no need to compute, same as forward
                // keep track of the min distance
                if(d < min_distance)
                {   min_distance = d ;
                    min_distances = std::vector<std::vector<size_t>>() ;
                    min_distances.push_back({from_x, from_y_rev, Constants::REVERSE}) ;
                }
                else if(d == min_distance)
                {   min_distances.push_back({from_x, from_y_rev, Constants::REVERSE}) ; }
            }
        }
    }

    // if several best distances, take the one with the min score
    size_t best_shift_x = 0 ;
    size_t best_shift_y = 0 ;
    size_t best_flip_y  = 0 ;
    int min_score = std::numeric_limits<int>::max() ;
    for(const auto& triplet : min_distances)
    {   size_t s1 = triplet[0] ;
        size_t s2 = triplet[1] ;
        size_t f  = triplet[2] ;
        if(this->_scores[s1][s2] < min_score)
        {   min_score    = this->_scores[s1][s2] ;
            best_shift_x = s1 ;
            best_shift_y = s2 ;
            best_flip_y  = f ;
        }
    }

    shift_x = best_shift_x ; shift_y = best_shift_y ; flip = best_flip_y ;
    return min_distance ;
}


bool CorDistanceComputer::reallocate_buffers(size_t shift) throw (std::invalid_argument)
{   if(shift == 0)
    {   throw std::invalid_argument("error! shift should be > 0!") ; }

    // need to reallocate
    bool realloc = false ;
    if(shift != this->_shift)
    {   this->_shift     = shift ;
        this->_distances = matrix3d_double(this->_shift, matrix2d_double(this->_shift, vector_double(2, 3.))) ;
        this->_sum_xy    = matrix3d_double(this->_shift, matrix2d_double(this->_shift, vector_double(2, 0.))) ;
        this->_sum_x     = vector_double(this->_shift, 0.) ;
        this->_sum_x2    = vector_double(this->_shift, 0.) ;
        this->_sum_y     = vector_double(this->_shift, 0.) ;
        this->_sum_y2    = vector_double(this->_shift, 0.) ;
        this->_scores    = matrix2d_int(shift, vector_int(shift)) ;
        realloc = true ;
    }

    return realloc ;
}


bool CorDistanceComputer::correct_cor(double& x)
{   if(isNaN(x))
    {   x = Constants::corr_fail ;
        return true ;
    }
    else if(roundToZero(x))
    {   return true ; }
    else if(roundToOne(x))
    {   return true ; }
    return false ;
}

