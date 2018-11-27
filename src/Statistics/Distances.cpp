#include "Distances.hpp"

#include <vector>
#include <cmath>
#include <limits>                    // numeric_limits<double>::epsilon() numeric_limits<size_t>::max()
#include <numeric>                   // std::accumulate()
#include <Statistics/Statistics.hpp>
#include <Utility/Utility.hpp>
#include <Utility/Constants.hpp>
#include <Utility/Vector_utility.hpp>


using matrix2d     = std::vector<std::vector<double>> ;
using matrix3d     = std::vector<std::vector<std::vector<double>>> ;
using vector       = std::vector<double> ;
using matrix2d_int = std::vector<std::vector<int>> ;
using vector_int   = std::vector<int> ;


double cor_distance(const std::vector<double>& x,
                    const std::vector<double>& y)
{
    assert(x.size() == y.size()) ;

    double cor = cor_pearson(x, y) ;
    correctCor(cor) ;
    double d = 1. - cor ; roundToZero(d) ;
    return d ;
}


double cor_distance(const std::vector<double>& x,
                    const std::vector<double>& y,
                    size_t shift,
                    size_t& shift1,
                    size_t& shift2)
{
    size_t l1 = x.size() ;
    size_t l2 = y.size() ;

    matrix2d     distances(shift, vector(shift, 0.)) ;
    matrix2d_int scores(shift, vector_int(shift, 0)) ;

    int l1_half = l1 / 2 ; // half x size
    int l2_half = l2 / 2 ; // half y size
    int n_half  = (l1 - shift + 1) / 2; // l1 - shift + 1 is the lenght of a vector slice taking the shift into account

    // to keep track of the best distance
    double min_distance  = 3. ;
    std::vector<std::vector<size_t>> min_distances ;


    // 1) compute the distances
    // 2) compute a score for each alignment -> sum of absolute difference of distance
    //    of middle sub-vector to middle vector, if several best distances, this will
    //    help choosing the one originating from the most central vector slices.
    // 3) keep track of the best distances and the shifts values (all the pairs leading to a dist
    //    equal to the minimal dist)
    // Do this for 1st row, 1st col and all the remaining others to ensure the same order as in cor_distance_fast().
    // Doing it this way ensure that the same optimal distance and shifts will be reported (the approach is greedy).

    // first row
    size_t s1     = 0 ;
    size_t x_from = s1 ;
    size_t x_to   = l1 - shift + x_from + 1 ;
    for(size_t s2=0; s2<shift; s2++)
    {   size_t y_from = s2 ;
        size_t y_to   = l2 - shift + y_from + 1 ;
        // distance
        double cor = cor_pearson(x ,y, x_from, x_to, y_from, y_to) ;
        correctCor(cor) ;
        double d = 1 - cor ; roundToZero(d) ;
        distances[s1][s2] = d ;
        // score
        scores[s1][s2] = abs(l1_half - (n_half+s1)) +  // the absolute distance of the middle of x slice to middle of x
                         abs(l2_half - (n_half+s2)) ;  // the absolute distance of the middle of y slice to middle of y
        // found new minimum
        if(distances[s1][s2] < min_distance)
        {   min_distance = distances[s1][s2] ;
            min_distances = std::vector<std::vector<size_t>>() ;
            min_distances.push_back(std::vector<size_t> {s1, s2}) ;
        }
        // found value equal to current minimum
        else if(distances[s1][s2] == min_distance)
        {   min_distances.push_back(std::vector<size_t> {s1, s2}) ; }
    }

    // first col
    size_t s2     = 0 ;
    size_t y_from = s2 ;
    size_t y_to   = l2 - shift + y_from + 1 ;
    for(size_t s1=1; s1<shift; s1++)
    {   size_t x_from = s1 ;
        size_t x_to   = l1 - shift + x_from + 1 ;
        // distance
        double cor = cor_pearson(x ,y, x_from, x_to, y_from, y_to) ;
        correctCor(cor) ;
        double d = 1 - cor ; roundToZero(d) ;
        distances[s1][s2] = d ;
        // score
        scores[s1][s2] = abs(l1_half - (n_half+s1)) +  // the absolute distance of the middle of x slice to middle of x
                         abs(l2_half - (n_half+s2)) ;  // the absolute distance of the middle of y slice to middle of y
        // found new minimum
        if(distances[s1][s2] < min_distance)
        {   min_distance = distances[s1][s2] ;
            min_distances = std::vector<std::vector<size_t>>() ;
            min_distances.push_back(std::vector<size_t> {s1, s2}) ;
        }
        // found value equal to current minimum
        else if(distances[s1][s2] == min_distance)
        {   min_distances.push_back(std::vector<size_t> {s1, s2}) ; }
    }


    // other
    for(size_t s1=1; s1<shift; s1++)
    {   size_t x_from = s1 ;
        size_t x_to   = l1 - shift + x_from + 1 ;
        for(size_t s2=1; s2<shift; s2++)
        {   size_t y_from = s2 ;
            size_t y_to   = l2 - shift + y_from + 1 ;
            // distance
            double cor = cor_pearson(x ,y, x_from, x_to, y_from, y_to) ;
            correctCor(cor) ;
            double d = 1 - cor ; roundToZero(d) ;
            distances[s1][s2] = d ;
            // score
            scores[s1][s2] = abs(l1_half - (n_half+s1)) +  // the absolute distance of the middle of x slice to middle of x
                             abs(l2_half - (n_half+s2)) ;  // the absolute distance of the middle of y slice to middle of y
            // found new minimum
            if(distances[s1][s2] < min_distance)
            {   min_distance = distances[s1][s2] ;
                min_distances = std::vector<std::vector<size_t>>() ;
                min_distances.push_back(std::vector<size_t> {s1, s2}) ;
            }
            // found value equal to current minimum
            else if(distances[s1][s2] == min_distance)
            {   min_distances.push_back(std::vector<size_t> {s1, s2}) ; }
        }
    }

    // search which of the minimum distances has the lowest score -> this will be distance to return
    int min_score = std::numeric_limits<double>::max() ;
    size_t best_s1 = 0 ; size_t best_s2 = 0 ;
    for(const auto& doublet : min_distances)
    {   size_t s1   = doublet[0] ;
        size_t s2   = doublet[1] ;
        if(scores[s1][s2] < min_score)
        {   min_score = scores[s1][s2] ;
            best_s1 = s1 ; best_s2 = s2 ;
        }
    }

    shift1 = best_s1 ; shift2 = best_s2 ;
    return distances[best_s1][best_s2] ;
}

double cor_distance(const std::vector<double>& x,
                    const std::vector<double>& y,
                    size_t shift,
                    size_t& shift1,
                    size_t& shift2,
                    bool& flip)
{
    size_t l1 = x.size() ;
    size_t l2 = y.size() ;

    matrix3d     distances(shift, matrix2d(shift, vector(2, 3))) ;
    matrix2d_int scores(shift, vector_int(shift, 0)) ;

    int l1_half = l1 / 2 ; // half v1 size
    int l2_half = l2 / 2 ; // half v2 size
    int n_half  = (l1 - shift + 1) / 2; // l1 - shift + 1 is the lenght of a vector slice taking the shift into account

    // to keep track of the best distance
    double min_distance  = 3. ;
    std::vector<std::vector<size_t>> min_distances ;

    // 1) compute the distances
    // 2) compute a score for each alignment -> sum of absolute difference of distance
    //    of middle sub-vector to middle vector, if several best distances, this will
    //    help choosing the one originating from the most central vector slices.
    // 3) keep track of the best distances and the shifts values (all the pairs leading to a dist
    //    equal to the minimal dist)
    // Do this for 1st row, 1st col and all the remaining others to ensure the same order as in cor_distance_fast().
    // Doing it this way ensure that the same optimal distance and shifts will be reported (the approach is greedy).

    // first row
    size_t s1      = 0 ;
    size_t v1_from = s1 ;
    size_t v1_to   = l1 - shift + v1_from + 1 ;
    for(size_t s2=0; s2<shift; s2++)
    {   size_t v2_from = s2 ;
        size_t v2_to   = l2 - shift + v2_from + 1 ;
        // distance
        // fw vs fw
        double cor = cor_pearson(x ,y, v1_from, v1_to, v2_from, v2_to) ;
        correctCor(cor) ;
        double d = 1 - cor ; roundToZero(d) ;
        distances[s1][s2][Constants::FORWARD] = d ;
        // fw vs rev
        cor = cor_pearson_rev(x, y, v1_from, v1_to, v2_to-1, v2_from-1) ;
        correctCor(cor) ;
        // in some cases dist should be zero but is rounded to the +/- epsilon machine, fix this issue
        d = 1 - cor ; roundToZero(d) ;
        distances[s1][s2][Constants::REVERSE] = d ;
        // score
        scores[s1][s2] = abs(l1_half - (n_half+static_cast<int>(s1))) +  // the absolute distance of the middle of v1 slice to middle of v1
                         abs(l2_half - (n_half+static_cast<int>(s2))) ;  // the absolute distance of the middle of v2 slice to middle of v2
        // found new minimum
        if(distances[s1][s2][Constants::FORWARD] < min_distance)
        {   min_distance = distances[s1][s2][Constants::FORWARD] ;
            min_distances = std::vector<std::vector<size_t>>() ;
            min_distances.push_back(std::vector<size_t> {s1, s2, Constants::FORWARD}) ;
        }
        // found value equal to current minimum
        else if(distances[s1][s2][Constants::FORWARD] == min_distance)
        {   min_distances.push_back(std::vector<size_t> {s1, s2, Constants::FORWARD}) ; }

        // found new minimum
        if(distances[s1][s2][Constants::REVERSE] < min_distance)
        {   min_distance = distances[s1][s2][Constants::REVERSE] ;
            min_distances = std::vector<std::vector<size_t>>() ;
            min_distances.push_back(std::vector<size_t> {s1, s2, Constants::REVERSE}) ;
        }
        // found value equal to current minimum
        else if(distances[s1][s2][Constants::REVERSE] == min_distance)
        {   min_distances.push_back(std::vector<size_t> {s1, s2, Constants::REVERSE}) ; }

    }

    // first column
    size_t s2      = 0 ;
    size_t v2_from = s2 ;
    size_t v2_to   = l2 - shift + v2_from + 1 ;
    for(size_t s1=1; s1<shift; s1++)
    {   size_t v1_from = s1 ;
        size_t v1_to   = l1 - shift + v1_from + 1 ;

        // distance
        // fw vs fw
        double cor = cor_pearson(x ,y, v1_from, v1_to, v2_from, v2_to) ;
        correctCor(cor) ;
        double d = 1 - cor ; roundToZero(d) ;
        distances[s1][s2][Constants::FORWARD] = d ;
        // fw vs rev
        cor = cor_pearson_rev(x, y, v1_from, v1_to, v2_to-1, v2_from-1) ;
        correctCor(cor) ;
        // in some cases dist should be zero but is rounded to the +/- epsilon machine, fix this issue
        d = 1 - cor ; roundToZero(d) ;
        distances[s1][s2][Constants::REVERSE] = d ;
        // score
        scores[s1][s2] = abs(l1_half - (n_half+static_cast<int>(s1))) +  // the absolute distance of the middle of v1 slice to middle of v1
                         abs(l2_half - (n_half+static_cast<int>(s2))) ;  // the absolute distance of the middle of v2 slice to middle of v2
        // found new minimum
        if(distances[s1][s2][Constants::FORWARD] < min_distance)
        {   min_distance = distances[s1][s2][Constants::FORWARD] ;
            min_distances = std::vector<std::vector<size_t>>() ;
            min_distances.push_back(std::vector<size_t> {s1, s2, Constants::FORWARD}) ;
        }
        // found value equal to current minimum
        else if(distances[s1][s2][Constants::FORWARD] == min_distance)
        {   min_distances.push_back(std::vector<size_t> {s1, s2, Constants::FORWARD}) ; }

        // found new minimum
        if(distances[s1][s2][Constants::REVERSE] < min_distance)
        {   min_distance = distances[s1][s2][Constants::REVERSE] ;
            min_distances = std::vector<std::vector<size_t>>() ;
            min_distances.push_back(std::vector<size_t> {s1, s2, Constants::REVERSE}) ;
        }
        // found value equal to current minimum
        else if(distances[s1][s2][Constants::REVERSE] == min_distance)
        {   min_distances.push_back(std::vector<size_t> {s1, s2, Constants::REVERSE}) ; }
    }

    // other
    for(size_t s1=1; s1<shift; s1++)
    {   size_t v1_from = s1 ;
        size_t v1_to   = l1 - shift + v1_from + 1 ;
        for(size_t s2=1; s2<shift; s2++)
        {   size_t v2_from = s2 ;
            size_t v2_to   = l2 - shift + v2_from + 1 ;
            // distance
            // fw vs fw
            double cor = cor_pearson(x ,y, v1_from, v1_to, v2_from, v2_to) ;
            correctCor(cor) ;
            double d = 1 - cor ; roundToZero(d) ;
            distances[s1][s2][Constants::FORWARD] = d ;
            // fw vs rev
            cor = cor_pearson_rev(x, y, v1_from, v1_to, v2_to-1, v2_from-1) ;
            correctCor(cor) ;
            // in some cases dist should be zero but is rounded to the +/- epsilon machine, fix this issue
            d = 1 - cor ; roundToZero(d) ;
            distances[s1][s2][Constants::REVERSE] = d ;
            // score
            scores[s1][s2] = abs(l1_half - (n_half+static_cast<int>(s1))) +  // the absolute distance of the middle of v1 slice to middle of v1
                             abs(l2_half - (n_half+static_cast<int>(s2))) ;  // the absolute distance of the middle of v2 slice to middle of v2
            // found new minimum
            if(distances[s1][s2][Constants::FORWARD] < min_distance)
            {   min_distance = distances[s1][s2][Constants::FORWARD] ;
                min_distances = std::vector<std::vector<size_t>>() ;
                min_distances.push_back(std::vector<size_t> {s1, s2, Constants::FORWARD}) ;
            }
            // found value equal to current minimum
            else if(distances[s1][s2][Constants::FORWARD] == min_distance)
            {   min_distances.push_back(std::vector<size_t> {s1, s2, Constants::FORWARD}) ; }

            // found new minimum
            if(distances[s1][s2][Constants::REVERSE] < min_distance)
            {   min_distance = distances[s1][s2][Constants::REVERSE] ;
                min_distances = std::vector<std::vector<size_t>>() ;
                min_distances.push_back(std::vector<size_t> {s1, s2, Constants::REVERSE}) ;
            }
            // found value equal to current minimum
            else if(distances[s1][s2][Constants::REVERSE] == min_distance)
            {   min_distances.push_back(std::vector<size_t> {s1, s2, Constants::REVERSE}) ; }
        }
    }

    // search which of the minimum distances has the lowest score -> this will be distance to return
    int min_score = std::numeric_limits<double>::max() ;
    size_t best_s1 = 0 ; size_t best_s2 = 0 ; size_t best_flip = 0 ;
    for(const auto& triplet : min_distances)
    {   size_t s1   = triplet[0] ;
        size_t s2   = triplet[1] ;
        size_t flip = triplet[2] ;
        if(scores[s1][s2] < min_score)
        {   min_score = scores[s1][s2] ;
            best_s1 = s1 ; best_s2 = s2 ; best_flip = flip ;
        }
    }

    shift1 = best_s1 ; shift2 = best_s2 ; flip = best_flip ;
    return distances[best_s1][best_s2][best_flip] ;
}


double cor_distance_fast(const std::vector<double>& x,
                         const std::vector<double>& y,
                         size_t shift,
                         size_t& shift1,
                         size_t& shift2)
{
    // the length of the vectors
    size_t l   = x.size() ;
    int l_half = l / 2 ;
    // the lenght of a vector slice taking the shift into account
    size_t n    = l - shift + 1 ;
    int n_half  = n / 2;

    assert(y.size() == l) ;
    assert(n >= 1) ;

    // data structures
    // for distance computations
    // 1-cor. of x for each shift (on the rows) with y for each shift (on the columns)
    matrix2d distances(shift, vector(shift, 3.)) ;
    // the sum of pairwise product of x elements for each shift (rows) and y elements for each shift (columns)
    matrix2d sum_xy(shift, vector(shift, 0.)) ;
    // the sum of x elements with a given shift
    vector sum_x(shift, 0.) ;
    // the sum of the square of the x elements with for each shift
    vector sum_x2(shift, 0.) ;
    // the sum of y elements with for each shift
    vector sum_y(shift, 0.) ;
    // the sum of the square of the y elements for each shift
    vector sum_y2(shift, 0.) ;

    // compute a score for each alignment -> sum of absolute difference of distance
    // of middle sub-vector to middle vector, if several best distances, this will
    // help choosing the one originating from the most central vector slices.
    matrix2d_int scores(shift, vector_int(shift)) ;

    // some iterators on x and y -> always x[from_x,to_x) and y[from_y,to_y)
    size_t from_x, to_x ;
    size_t from_y, to_y ;

    // to track the minimum distance
    double min_distance  = 3. ;                      // the current minimum distance value
    std::vector<std::vector<size_t>> min_distances ; // the coordinates of all distance equal to the current min distance

    // initialize the sums of x and y elements for each shift
    for(size_t i=0; i<n; i++)
    {   sum_x[0]  += x[i] ;
        sum_x2[0] += pow(x[i],2) ;
        sum_y[0]  += y[i] ;
        sum_y2[0] += pow(y[i],2) ;
    }
    for(size_t s=1, i=n; s<shift; s++, i++)
    {   sum_x[s]  += sum_x[s-1]  - x[i-n]        + x[i] ;
        sum_x2[s] += sum_x2[s-1] - pow(x[i-n],2) + pow(x[i],2) ;
        sum_y[s]  += sum_y[s-1]  - y[i-n]        + y[i] ;
        sum_y2[s] += sum_y2[s-1] - pow(y[i-n],2) + pow(y[i],2) ;
    }

    // fill 1st row of dist -> compute 1-cor with x having a shift of 0
    from_x  = 0 ;          // from_x is also the shift of x
    to_x    = from_x + n ; // x[from_x,to_x)
    for(size_t i=0; i<shift; i++) // i is the shift of y
    {   from_y = i ;
        to_y   = from_y + n ; // y[from_y,to_y)

        sum_xy[from_x][from_y] = 0. ;
        for(size_t j=0; j<n; j++)
        {   sum_xy[from_x][from_y] += x[from_x+j]*y[from_y+j] ; }

        double cor = (n*sum_xy[from_x][from_y] - sum_x[from_x]*sum_y[from_y]) /
                     (sqrt(n*sum_x2[from_x] - pow(sum_x[from_x], 2))*sqrt(n*sum_y2[from_y] - pow(sum_y[from_y], 2))) ;
        if(isNaN(cor))
        {   cor = 0. ; }
        double d = 1. - cor ;
        // in some cases d_tmp should be zero but is rounded to the +/- epsilon machine, fix this issue
        if(d <= std::numeric_limits<double>::epsilon() and d >= -std::numeric_limits<double>::epsilon())
        {   d = 0. ; }
        // the distance
        distances[from_x][from_y] = d ;
        // the distance score
        scores[from_x][from_y] = abs(l_half - (n_half+static_cast<int>(from_x))) +  // the absolute distance of the middle of v1 slice to middle of v1
                                 abs(l_half - (n_half+static_cast<int>(from_y))) ;  // the absolute distance of the middle of v2 slice to middle of v2
        // keep track of the min distance
        if(d < min_distance)
        {   min_distance = d ;
            min_distances = std::vector<std::vector<size_t>>() ;
            min_distances.push_back(std::vector<size_t>({from_x, from_y})) ;
        }
        else if(d == min_distance)
        {   min_distances.push_back(std::vector<size_t>({from_x, from_y})) ; }
    }

    // fill 1st column of dist -> compute the 1-cor with y having a shift of 0
    from_y = 0 ;            // from_y is also the shift of y
    to_y   = from_y + n ;   // y[from_y,to_y)
    for(size_t i=1; i<shift; i++)  // i is the shift of x
    {   from_x = i ;
        to_x   = from_x + n ; //x[from_x,to_x)

        sum_xy[from_x][from_y] = 0. ;
        for(size_t j=0; j<n; j++)
        {   sum_xy[from_x][from_y] += x[from_x+j]*y[from_y+j] ; }

        double cor = (n*sum_xy[from_x][from_y] - sum_x[from_x]*sum_y[from_y]) /
                     (sqrt(n*sum_x2[from_x] - pow(sum_x[from_x], 2))*sqrt(n*sum_y2[from_y] - pow(sum_y[from_y], 2))) ;
        if(isNaN(cor))
        {   cor = 0. ; }
        double d = 1. - cor ;
        // in some cases d_tmp should be zero but is rounded to the +/- epsilon machine, fix this issue
        if(d <= std::numeric_limits<double>::epsilon() and d >= -std::numeric_limits<double>::epsilon())
        {   d = 0. ; }
        distances[from_x][from_y] = d ;

        // the distance score
        scores[from_x][from_y] = abs(l_half - (n_half+static_cast<int>(from_x))) +  // the absolute distance of the middle of v1 slice to middle of v1
                                 abs(l_half - (n_half+static_cast<int>(from_y))) ;  // the absolute distance of the middle of v2 slice to middle of v2
        // keep track of the min distance
        if(d < min_distance)
        {   min_distance = d ;
            min_distances = std::vector<std::vector<size_t>>() ;
            min_distances.push_back(std::vector<size_t>({from_x, from_y})) ;
        }
        else if(d == min_distance)
        {   min_distances.push_back(std::vector<size_t>({from_x, from_y})) ; }
    }

    // compute the remaining distances
    for(size_t i=1; i<shift; i++)
    {   from_x = i ;
        to_x   = from_x + n ; // x[from_x,to_x)
        for(size_t j=1; j<shift; j++)
        {   from_y = j ;
            to_y   = from_y + n ; // y[from_y,to_y)
            // compute the current some of pairwise product from a previous one
            sum_xy[from_x][from_y] = sum_xy[from_x-1][from_y-1] - x[from_x-1]*y[from_y-1] + x[to_x-1]*y[to_y-1] ;
            double cor             = (n*sum_xy[from_x][from_y] - sum_x[from_x]*sum_y[from_y]) /
                                     (sqrt(n*sum_x2[from_x] - pow(sum_x[from_x], 2))*sqrt(n*sum_y2[from_y] - pow(sum_y[from_y], 2))) ;
            if(isNaN(cor))
            {   cor = 0. ; }
            double d = 1. - cor ;
            if(d <= std::numeric_limits<double>::epsilon() and d >= -std::numeric_limits<double>::epsilon())
            {   d = 0. ; }
            distances[from_x][from_y] = d ;
            // the distance score
            scores[from_x][from_y] = abs(l_half - (n_half+static_cast<int>(from_x))) +  // the absolute distance of the middle of v1 slice to middle of v1
                                     abs(l_half - (n_half+static_cast<int>(from_y))) ;  // the absolute distance of the middle of v2 slice to middle of v2
            // keep track of the min distance
            if(d < min_distance)
            {   min_distance = d ;
                min_distances = std::vector<std::vector<size_t>>() ;
                min_distances.push_back(std::vector<size_t>({from_x, from_y})) ;
            }
            else if(d == min_distance)
            {   min_distances.push_back(std::vector<size_t>({from_x, from_y})) ; }
        }
    }

    // search the best comparison, do this here to be sure to report the lowest possible shift1 and shift2
    min_distance = 3. ;
    min_distances = std::vector<std::vector<size_t>>() ;
    for(size_t i=0; i<shift; i++)
    {   for(size_t j=0; j<shift; j++)
        {   if(distances[i][j] < min_distance)
            {   min_distance = distances[i][j] ;
                min_distances = std::vector<std::vector<size_t>>() ;
                min_distances.push_back(std::vector<size_t>({i, j})) ;
            }
            else if(distances[i][j] == min_distance)
            {   min_distances.push_back(std::vector<size_t>({i, j})) ; }
        }
    }

    // if several best distances, take the one with the min score
    size_t best_shift1 = 0 ;
    size_t best_shift2 = 0 ;
    int min_score = std::numeric_limits<int>::max() ;
    for(const auto& doublet : min_distances)
    {   size_t s1 = doublet[0] ;
        size_t s2 = doublet[1] ;
        if(scores[s1][s2] < min_score)
        {   min_score   = scores[s1][s2] ;
            best_shift1 = s1 ;
            best_shift2 = s2 ;
        }
    }

    // return results
    shift1 = best_shift1 ; shift2 = best_shift2 ;
    return min_distance ;
}


double cor_distance_fast(const std::vector<double>& x,
                    const std::vector<double>& y,
                    size_t shift,
                    size_t& shift1,
                    size_t& shift2,
                    bool& flip)
{
    // the length of the vectors
    size_t l   = x.size() ;
    int l_half = l / 2 ;
    // the lenght of a vector slice taking the shift into account
    size_t n    = l - shift + 1 ;
    int n_half  = static_cast<int>(n) / 2;

    assert(y.size() == l) ;
    assert(n >= 1) ;

    // data structures
    // 1-cor. of x for each shift (on the rows) with y for each shift (on the columns) and on each flip
    matrix3d distances(shift, std::vector<std::vector<double>>(shift, std::vector<double>(2, 3.))) ;
    // the sum of pairwise product of x elements for each shift (rows) and y elements for each shift (columns)
    matrix3d sum_xy(shift, std::vector<std::vector<double>>(shift, std::vector<double>(2, 0.))) ;
    // the sum of x elements with a given shift
    vector sum_x(shift, 0.) ;
    // the sum of the square of the x elements with for each shift
    vector sum_x2(shift, 0.) ;
    // the sum of y elements with for each shift (the flip state does not influence these values)
    vector sum_y(shift, 0.) ;
    // the sum of the square of the y elements for each shift (the flip state does not influence these values)
    vector sum_y2(shift, 0.) ;

    // compute a score for each alignment -> sum of absolute difference of distance
    // of middle sub-vector to middle vector, if several best distances, this will
    // help choosing the one originating from the most central vector slices.
    matrix2d_int scores(shift, vector_int(shift)) ;

    // some iterators on x and y -> always x[from_x,to_x) and y[from_y,to_y)
    size_t from_x, to_x ;
    size_t from_y, to_y, from_y_rev, to_y_rev ;

    // to track the minimum distance
    double min_distance  = 3. ;                      // the current minimum distance value
    std::vector<std::vector<size_t>> min_distances ; // the coordinates of all distance equal to the current min distance

    // initialize the sums of x and y elements for each shift
    for(size_t i=0; i<n; i++)
    {   sum_x[0]  += x[i]        ;
        sum_x2[0] += pow(x[i],2) ;
        sum_y[0]  += y[i]        ;
        sum_y2[0] += pow(y[i],2) ;
    }
    for(size_t s=1, i=n; s<shift; s++, i++)
    {   sum_x[s]  += sum_x[s-1]  - x[i-n]        + x[i] ;
        sum_x2[s] += sum_x2[s-1] - pow(x[i-n],2) + pow(x[i],2) ;
        sum_y[s]  += sum_y[s-1]  - y[i-n]        + y[i] ;
        sum_y2[s] += sum_y2[s-1] - pow(y[i-n],2) + pow(y[i],2) ;
    }

    // fill 1st row of dist -> compute the 1-cor with x having a shift of 0
    from_x = 0 ;
    to_x   = from_x + n ; // x[from_x,to_x)
    for(size_t i=0; i<shift; i++) // i is the shift of y
    {   from_y = i ;
        to_y   = from_y + n ; // y[from_y,to_y)

        // forward orientation -> x fw vs y fw
        {
            sum_xy[from_x][from_y][Constants::FORWARD] = 0. ;
            for(size_t j=0; j<n; j++)
            {   sum_xy[from_x][from_y][Constants::FORWARD] += x[from_x+j]*y[from_y+j] ; }


            double cor = (n*sum_xy[from_x][from_y][Constants::FORWARD] - sum_x[from_x]*sum_y[from_y]) /
                         (sqrt(n*sum_x2[from_x] - pow(sum_x[from_x], 2))*sqrt(n*sum_y2[from_y] - pow(sum_y[from_y], 2))) ;
            correctCor(cor) ;
            double d = 1. - cor ; roundToZero(d) ;
            // the distance
            distances[from_x][from_y][Constants::FORWARD] = d ;
            // the distance score
            scores[from_x][from_y] = abs(l_half - (n_half+static_cast<int>(from_x))) +  // the absolute distance of the middle of v1 slice to middle of v1
                                     abs(l_half - (n_half+static_cast<int>(from_y))) ;  // the absolute distance of the middle of v2 slice to middle of v2
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
        {   sum_xy[from_x][from_y][Constants::REVERSE] = 0. ;
            for(size_t j=0; j<n; j++)
            {   sum_xy[from_x][from_y][Constants::REVERSE] += x[from_x+j]*y[to_y-j-1] ; }

            double cor = (n*sum_xy[from_x][from_y][Constants::REVERSE] - sum_x[from_x]*sum_y[from_y]) /
                         (sqrt(n*sum_x2[from_x] - pow(sum_x[from_x], 2))*sqrt(n*sum_y2[from_y] - pow(sum_y[from_y], 2))) ;
            correctCor(cor) ;
            double d = 1. - cor ; roundToZero(d) ;
            // distance
            distances[from_x][from_y][Constants::REVERSE] = d ;
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
    for(size_t i=0; i<shift; i++)     // i is the shift of x
    {   from_x = i ;
        to_x   = from_x + n ; // x[from_x,to_x)
        // forward orientation -> x fw vs y fw
        {   sum_xy[from_x][from_y][Constants::FORWARD] = 0. ;
            for(size_t j=0; j<n; j++)
            {   sum_xy[from_x][from_y][Constants::FORWARD] += x[from_x+j]*y[from_y+j] ; }

            double cor = (n*sum_xy[from_x][from_y][Constants::FORWARD] - sum_x[from_x]*sum_y[from_y]) /
                         (sqrt(n*sum_x2[from_x] - pow(sum_x[from_x], 2))*sqrt(n*sum_y2[from_y] - pow(sum_y[from_y], 2))) ;
            correctCor(cor) ;
            double d = 1. - cor ; roundToZero(d) ;
            // distance
            distances[from_x][from_y][Constants::FORWARD] = d ;
            // the distance score
            scores[from_x][from_y] = abs(l_half - (n_half+static_cast<int>(from_x))) +  // the absolute distance of the middle of v1 slice to middle of v1
                                     abs(l_half - (n_half+static_cast<int>(from_y))) ;  // the absolute distance of the middle of v2 slice to middle of v2
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
        {   sum_xy[from_x][shift-1][Constants::REVERSE] = 0. ;
            for(size_t j=0; j<n; j++)
            {   sum_xy[from_x][shift-1][Constants::REVERSE] += x[from_x+j]*y[to_y_rev-1-j] ; }

            double cor = (n*sum_xy[from_x][shift-1][Constants::REVERSE] - sum_x[from_x]*sum_y[shift-1]) /
                         (sqrt(n*sum_x2[from_x] - pow(sum_x[from_x], 2))*sqrt(n*sum_y2[shift-1] - pow(sum_y[shift-1], 2))) ;
            correctCor(cor) ;
            double d = 1. - cor ; roundToZero(d) ;
            // distance
            distances[from_x][shift-1][Constants::REVERSE] = d ;
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
                sum_xy[from_x][from_y][Constants::FORWARD] = sum_xy[from_x-1][from_y-1][Constants::FORWARD] - x[from_x-1]*y[from_y-1] + x[to_x-1]*y[to_y-1] ;
                double cor                                 = (n*sum_xy[from_x][from_y][Constants::FORWARD] - sum_x[from_x]*sum_y[from_y]) /
                                                             (sqrt(n*sum_x2[from_x] - pow(sum_x[from_x], 2))*sqrt(n*sum_y2[from_y] - pow(sum_y[from_y], 2))) ;
                correctCor(cor) ;
                double d = 1. - cor ; roundToZero(d) ;
                // distance
                distances[from_x][from_y][Constants::FORWARD] = d ;
                // the distance score
                scores[from_x][from_y] = abs(l_half - (n_half+static_cast<int>(from_x))) +  // the absolute distance of the middle of v1 slice to middle of v1
                                         abs(l_half - (n_half+static_cast<int>(from_y))) ;  // the absolute distance of the middle of v2 slice to middle of v2
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
            {   size_t j_rev = shift - j - 1 ;    // gives the complementary shift state -> if shift=3 : 0<->2, 1<->1, 2<->0
                from_y_rev = j_rev  ;             // the 1st element of the slice if considering <shift> assigned at the end
                to_y_rev   = from_y_rev + n     ; // y[from_rev,to_rev)

                // compute the current some of pairwise product from a previous one
                sum_xy[from_x][from_y_rev][Constants::REVERSE] = sum_xy[from_x-1][from_y_rev+1][Constants::REVERSE] - x[from_x-1]*y[to_y_rev] + x[to_x-1]*y[from_y_rev] ;
                double cor                                     = (n*sum_xy[from_x][from_y_rev][Constants::REVERSE] - sum_x[from_x]*sum_y[from_y_rev]) /
                                                                 (sqrt(n*sum_x2[from_x] - pow(sum_x[from_x], 2))*sqrt(n*sum_y2[from_y_rev] - pow(sum_y[from_y_rev], 2))) ;
                correctCor(cor) ;
                double d = 1. - cor ; roundToZero(d) ;
                // distance
                distances[from_x][from_y_rev][Constants::REVERSE] = d ;
                // the distance score
                // no need to compute, same as forward
                // keep track of the min distance
                if(d < min_distance)
                {   min_distance = d ;
                    min_distances = std::vector<std::vector<size_t>>() ;
                    min_distances.push_back({from_x, from_y_rev, Constants::REVERSE}) ;
                }
                else if(d == min_distance)
                // else if(isEqual(d, min_distance, 1e-6))
                {   min_distances.push_back({from_x, from_y_rev, Constants::REVERSE}) ; }
            }
        }
    }

    // if several best distances, take the one with the min score
    size_t best_shift1 = 0 ;
    size_t best_shift2 = 0 ;
    size_t best_flip   = 0 ;
    int min_score = std::numeric_limits<int>::max() ;
    for(const auto& triplet : min_distances)
    {   size_t s1 = triplet[0] ;
        size_t s2 = triplet[1] ;
        size_t f  = triplet[2] ;
        if(scores[s1][s2] < min_score)
        {   min_score   = scores[s1][s2] ;
            best_shift1 = s1 ;
            best_shift2 = s2 ;
            best_flip   = f ;
        }
    }

    shift1 = best_shift1 ; shift2 = best_shift2 ; flip = best_flip ;
    return min_distance ;
}


bool correctCor(double& x)
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
