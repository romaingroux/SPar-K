
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>                // reverse()
#include <cstdlib>
#include <stdexcept>
#include <UnitTest++/UnitTest++.h>

#include "Statistics/Statistics.hpp"
#include "Statistics/Distances.hpp"
#include "Utility/Utility.hpp"
#include "Utility/Constants.hpp"
#include "Matrix/Matrix2D.hpp"

using namespace std;


// Statistics.cpp unittests
SUITE(Statistics_testsuit)
{
    TEST(start_message)
    {   cout << "Running Statistics tests..." << endl ; }

    // tests the mean() function
    TEST(mean_test)
    {   // tolerated error for equality testing
        double error = 0.0001 ;

        vector<double> v = {-2.5, 4.3, 5.0, -10, 33.0} ;

        // mean of the of the all vector
        double results1  = mean(v) ;
        double results2  = mean(v,-1,-1) ;
        double results3  = mean(v, 0, 5) ;
        double expected1 = 5.96 ;
        CHECK_CLOSE(results1, expected1, error) ;
        CHECK_CLOSE(results2, expected1, error) ;
        CHECK_CLOSE(results3, expected1, error) ;

        // mean over a single value
        double results4  = mean(v,3,4) ;
        double expected4 = -10 ;
        CHECK_CLOSE(results4, expected4, error) ;
        // mean over a part of the vector
        double results5  = mean(v,0,4) ;
        double expected5 = -0.8 ;
        CHECK_CLOSE(results5, expected5, error) ;
    }

    TEST(mean_weigthed)
    {   // tolerated error for equality testing
        double error = 0.0001 ;

        vector<double> x  = {-2.5, 4.3, 5.0, -10,  33.0} ;
        vector<double> p1 = {0.2,  0.2, 0.2,  0.2, 0.2} ;
        vector<double> p2 = {0.1,  0.2, 0.3,  0.4, 0.0} ;
        vector<double> p3 = {1.0,  2.0, 3.0,  4.0, 5.0} ;
        vector<double> p4 = {0.0,  0.0, 0.0,  0.0, 0.0} ;

        double expected1   = 5.96 ;
        double expected2   = -1.89 ;
        double expected3   = 9.74 ;
        bool   expected4   = true ;

        CHECK_CLOSE(mean(x, p1), expected1, error) ;
        CHECK_CLOSE(mean(x, p2), expected2, error) ;
        CHECK_CLOSE(mean(x, p3), expected3, error) ;
        CHECK_EQUAL(isNaN(mean(x, p4)), expected4) ; // NaN, division by 0
    }

    TEST(sd_test)
    {   // tolerated error for equality testing
        double error = 0.0001 ;

        vector<double> v = {-2.5, 4.3, 5.0, -10, 33.0} ;

        // sd of the of the all vector
        double results1  = sd(v) ;
        double results2  = sd(v,-1,-1) ;
        double results3  = sd(v, 0, 5) ;
        double expected1 = 16.28751 ;
        CHECK_CLOSE(results1, expected1, error) ;
        CHECK_CLOSE(results2, expected1, error) ;
        CHECK_CLOSE(results3, expected1, error) ;

        // mean over a single value
        double results4  = sd(v,3,4) ;
        double expected4 = nan("") ;
        CHECK_EQUAL(isNaN(results4), isNaN(expected4)) ;

        // mean over a part of the vector
        double results5  = sd(v,0,4) ;
        double expected5 = 7.004284 ;
        CHECK_CLOSE(results5, expected5, error) ;

    }

    TEST(sd_weigthed)
    {   // tolerated error for equality testing
        double error = 0.0001 ;

        vector<double> x  = {-2.5, 4.3, 5.0, -10,  33.0} ;
        vector<double> p1 = {0.2,  0.2, 0.2,  0.2, 0.2} ;
        vector<double> p2 = {0.1,  0.2, 0.3,  0.4, 0.0} ;
        vector<double> p3 = {1.0,  2.0, 3.0,  4.0, 5.0} ;
        vector<double> p4 = {0.0,  0.0, 0.0,  0.0, 0.0} ;

        double expected1   = 16.28751 ;
        double expected2   = 8.302401 ;
        double expected3   = 20.01518 ;
        bool   expected4   = true ;

        CHECK_CLOSE(sd(x, p1), expected1, error) ;
        CHECK_CLOSE(sd(x, p2), expected2, error) ;
        CHECK_CLOSE(sd(x, p3), expected3, error) ;
        CHECK_EQUAL(isNaN(sd(x, p4)), expected4) ; // NaN, division by 0
    }

    TEST(cor_pearson_test)
    {   // tolerated error for equality testing
        double error = 0.000001 ;

        vector<int> v1 = {0, 0, 1, 2, 3, 0, 0, 1} ;
        vector<int> v2 = {1, 2, 3, 0, 0, 0, 0, 1} ;

        // correlation over the whole vector
        double expected1 = -0.2394366 ;
        double results1  = cor_pearson(v1, v2) ;
        double results2  = cor_pearson(v1, v2, -1, -1, -1, -1) ;
        double results3  = cor_pearson(v1, v2,  0,  8,  0,  8) ;
        double results4  = cor_pearson(v2, v1) ;
        CHECK_CLOSE(expected1, results1, error) ;
        CHECK_EQUAL(results1, results2) ;
        CHECK_EQUAL(results2, results3) ;
        CHECK_EQUAL(results3, results4) ;

        // correlation over part of the vector
        double expected5 = 1.0 ;
        double results5  = cor_pearson(v1, v2, 2, 5, 0, 3) ;
        CHECK_CLOSE(expected5, results5, error) ;

        // correlation of non-variable vectors (sd1 == sd2 == 0)
        double expected6 = true ;
        double results6  = cor_pearson(v1, v2, 0, 2, 3, 5) ;
        CHECK_EQUAL(expected6, isNaN(results6)) ;

        // correlation using single values
        double expected7 = true ;
        double results7  = cor_pearson(v1, v2, 0, 1, 3, 4) ;
        CHECK_EQUAL(expected7, isNaN(results7)) ;
    }

    TEST(cor_pearson_rev_test)
    {   // tolerated error for equality testing
        double error = 0.000001 ;

        vector<int> v1 = {0, 0, 1, 2, 3, 0, 0, 1} ;
        vector<int> v2 = {1, 0, 0, 0, 0, 3, 2, 1} ;

        // correlation over the whole vector
        double expected1 = -0.2394366 ;
        bool   expected2 = true ;
        double results1  = cor_pearson_rev(v1, v2) ;
        double results2  = cor_pearson_rev(v1, v2, -1, -1, -1, -1) ;
        double results3  = cor_pearson_rev(v1, v2,  0,  8,  7, -1) ;
        double results4  = cor_pearson_rev(v2, v1) ;
        CHECK_CLOSE(expected1, results1, error) ;
        CHECK_EQUAL(results1,  results2) ;
        CHECK_EQUAL(results2,  results3) ;
        CHECK_EQUAL(expected2, results3 == results4) ;

        // correlation over part of the vector
        double expected5 = 1.0 ;
        double results5  = cor_pearson_rev(v1, v2, 2, 5, 7, 4) ;
        CHECK_CLOSE(expected5, results5, error) ;

        // correlation of non-variable vectors (sd1 == sd2 == 0)
        double expected6 = true ;
        double results6  = cor_pearson_rev(v1, v2, 0, 2, 4, 2) ;
        CHECK_EQUAL(expected6, isNaN(results6)) ;

        // correlation using single values
        double expected7 = true ;
        double results7  = cor_pearson(v1, v2, 0, 1, 3, 4) ;
        CHECK_EQUAL(expected7, isNaN(results7)) ;
    }
}


SUITE(Distances_testsuit)
{   /* this test suite checks that
       1) the naive way to compute correlation distances is correct and
       2) the fast way (using the dynamic programming algo) is totally
          equivalent to the naive approach
    */
    TEST(start_message)
    {   cout << "Running Distances tests..." << endl ; }

    TEST(cor_distance_test)
    {
        // tolerated error for equality testing, the expected values are
        // reported with 7 decimals only
        double error = 0.000001 ;
        // some variables for the distance functions
        double d, d_exp ;
        size_t shift1, shift2 ;
        bool flip2 ;

        // --------------------- vector are the same -> distance of 0 without even shifting ---------------------
        d_exp = 0. ;
        std::vector<double> v1 = {0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0} ;
        std::vector<double> v2 = {0, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 0} ;

        CHECK_CLOSE(cor_distance(v1, v2), d_exp, error) ;
        CHECK_CLOSE(cor_distance(v1, v2), cor_distance(v2, v1), error) ;

        // allow shifting (which is not required)
        d = cor_distance(v1, v2, 1, shift1, shift2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 0) ;

        d = cor_distance(v1, v2, 3, shift1, shift2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 1) ;
        CHECK_EQUAL(shift2, 1) ;

        d = cor_distance(v1, v2, 5, shift1, shift2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 2) ;
        CHECK_EQUAL(shift2, 2) ;

        // allow shift and flip (which are not required)
        d = cor_distance(v1, v2, 1, shift1, shift2, flip2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 0) ;
        CHECK_EQUAL(flip2, Constants::FORWARD) ;

        d = cor_distance(v1, v2, 3, shift1, shift2, flip2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 1) ;
        CHECK_EQUAL(shift2, 1) ;
        CHECK_EQUAL(flip2, Constants::FORWARD) ;

        d = cor_distance(v1, v2, 5, shift1, shift2, flip2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 2) ;
        CHECK_EQUAL(shift2, 2) ;
        CHECK_EQUAL(flip2, Constants::FORWARD) ;

        // --------------------- vector are the shifted -> distance of 0 without even shifting ---------------------
        std::vector<double> v3 = {1, 2, 3, 0, 4, 5, 5, 5, 0, 0, 1, 0, 0} ;
        std::vector<double> v4 = {0, 0, 1, 2, 3, 0, 4, 5, 5, 5, 0, 0, 1} ;

        d_exp = 0.8518519 ;
        CHECK_CLOSE(cor_distance(v3, v4), d_exp, error) ;
        CHECK_CLOSE(cor_distance(v3, v4), cor_distance(v4, v3), error) ;

        // allow shifting (shift=1 is equal to no shift, results should be as above)
        d = cor_distance(v3, v4, 1, shift1, shift2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 0) ;

        // allow the correct shifting
        d_exp = 0. ;
        d = cor_distance(v3, v4, 5, shift1, shift2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 2) ;

        // allow a bit more than the correct shifting
        d_exp = 0. ;
        d = cor_distance(v3, v4, 7, shift1, shift2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 1) ;
        CHECK_EQUAL(shift2, 3) ;

        // allow flip and shift (flip not required but shifting is not enough here)
        d_exp = 0.722222 ;
        d = cor_distance(v3, v4, 1, shift1, shift2, flip2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 0) ;
        CHECK_EQUAL(flip2, Constants::REVERSE) ;

        // allow the correct shifting
        d_exp = 0. ;
        d = cor_distance(v3, v4, 5, shift1, shift2, flip2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 2) ;
        CHECK_EQUAL(flip2, Constants::FORWARD) ;

        // allow a bit more than the correct shifting
        d_exp = 0. ;
        d = cor_distance(v3, v4, 7, shift1, shift2, flip2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 1) ;
        CHECK_EQUAL(shift2, 3) ;
        CHECK_EQUAL(flip2, Constants::FORWARD) ;

        // --------------------- vectors are flipped ---------------------
        std::vector<double> v5 = {1, 2, 3, 0, 4, 5, 5, 5, 0, 0, 1, 0, 0} ;
        std::vector<double> v6(v5) ; std::reverse(v6.begin(), v6.end()) ;

        d_exp = 0.462963 ;
        CHECK_CLOSE(cor_distance(v5, v6), d_exp, error) ;
        CHECK_CLOSE(cor_distance(v5, v6), cor_distance(v6, v5), error) ;

        // allow shifting (shift=1 is equal to no shift)
        d_exp = 0.462963 ;
        d = cor_distance(v5, v6, 1, shift1, shift2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 0) ;

        // allow more shifting
        d_exp = 0.113708 ;
        d = cor_distance(v5, v6, 5, shift1, shift2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 1) ;

        d_exp = 0.124405 ;
        d = cor_distance(v5, v6, 7, shift1, shift2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 1) ;

        // allow flip and shift (flip required but not shifting)
        d_exp = 0. ;
        d = cor_distance(v5, v6, 1, shift1, shift2, flip2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 0) ;
        CHECK_EQUAL(flip2, Constants::REVERSE) ;

        // allow more shifting
        d_exp = 0. ;
        d = cor_distance(v5, v6, 5, shift1, shift2, flip2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 2) ;
        CHECK_EQUAL(shift2, 2) ;
        CHECK_EQUAL(flip2, Constants::REVERSE) ;

        // allow even more shifting
        d_exp = 0. ;
        d = cor_distance(v5, v6, 7, shift1, shift2, flip2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 3) ;
        CHECK_EQUAL(shift2, 3) ;
        CHECK_EQUAL(flip2, Constants::REVERSE) ;

        // --------------------- vectors are shifted and flipped ---------------------
        std::vector<double> v7 = {1, 2, 3, 0, 4, 5, 5, 5, 0, 0, 1, 0, 0} ;
        std::vector<double> v8 = {0, 0, 1, 2, 3, 0, 4, 5, 5, 5, 0, 0, 1} ; std::reverse(v8.begin(), v8.end()) ;

        d_exp = 0.7222222 ;
        CHECK_CLOSE(cor_distance(v7, v8), d_exp, error) ;
        CHECK_CLOSE(cor_distance(v7, v8), cor_distance(v8, v7), error) ;

        // allow shifting (shift=1 is equal to no shift)
        d_exp = 0.7222222 ;
        d = cor_distance(v7, v8, 1, shift1, shift2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 0) ;

        // allow more shifting
        d_exp = 0.113708 ;
        d = cor_distance(v7, v8, 5, shift1, shift2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 3) ;
        CHECK_EQUAL(shift2, 2) ;

        // allow even more shifting
        d_exp = 0.124405 ;
        d = cor_distance(v7, v8, 7, shift1, shift2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 5) ;
        CHECK_EQUAL(shift2, 4) ;

        // allow flip and shift (flip required and shifting required)
        d_exp = 0.7222222 ;
        d = cor_distance(v7, v8, 1, shift1, shift2, flip2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 0) ;
        CHECK_EQUAL(flip2, Constants::FORWARD) ;

        // allow the good shifting
        // this is a bit counter intuitive that both shifts are 0 : it is because the vectors are
        // first sliced (to get the shift) and then reversed
        d_exp = 0. ;
        d = cor_distance(v7, v8, 3, shift1, shift2, flip2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 0) ;
        CHECK_EQUAL(flip2, Constants::REVERSE) ;

        // allow more shifting
        d_exp = 0. ;
        d = cor_distance(v7, v8, 5, shift1, shift2, flip2) ;
        CHECK_CLOSE(d, d_exp, error) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 2) ;
        CHECK_EQUAL(flip2, Constants::REVERSE) ;

        // --------------------- one vector as sd of 0 ---------------------

        std::vector<double> v9  = {1, 2, 3, 0, 4, 5, 5, 5, 0, 0, 1, 0, 0} ;
        std::vector<double> v10 = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1} ;

        // cor(v9,v10) should return Constants::corr_fail, so the distance should be
        // 1 - Constants::corr_fail
        d_exp = 1 - Constants::corr_fail ;
        CHECK_EQUAL(cor_distance(v9, v10), d_exp) ;
        CHECK_EQUAL(cor_distance(v9, v10), cor_distance(v10, v9)) ;

        // should also be the case with shifting (whatever the freedom)
        d = cor_distance(v9, v10, 1, shift1, shift2) ;
        CHECK_EQUAL(d, d_exp) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 0) ;

        d = cor_distance(v9, v10, 5, shift1, shift2) ;
        CHECK_EQUAL(d, d_exp) ;
        CHECK_EQUAL(shift1, 2) ;
        CHECK_EQUAL(shift2, 2) ;

        d = cor_distance(v9, v10, 7, shift1, shift2) ;
        CHECK_EQUAL(d, d_exp) ;
        CHECK_EQUAL(shift1, 3) ;
        CHECK_EQUAL(shift2, 3) ;

        // should also be the case with shifting (whatever the freedom) and flipping
        d = cor_distance(v9, v10, 1, shift1, shift2, flip2) ;
        CHECK_EQUAL(d, d_exp) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 0) ;
        CHECK_EQUAL(flip2, Constants::FORWARD) ;

        d = cor_distance(v9, v10, 5, shift1, shift2, flip2) ;
        CHECK_EQUAL(d, d_exp) ;
        CHECK_EQUAL(shift1, 2) ;
        CHECK_EQUAL(shift2, 2) ;
        CHECK_EQUAL(flip2, Constants::FORWARD) ;

        d = cor_distance(v9, v10, 7, shift1, shift2, flip2) ;
        CHECK_EQUAL(d, d_exp) ;
        CHECK_EQUAL(shift1, 3) ;
        CHECK_EQUAL(shift2, 3) ;
        CHECK_EQUAL(flip2, Constants::FORWARD) ;


        // compare the naive method (tested above) and the fast method
        size_t shift1_naive, shift2_naive ; bool flip_naive ; double d_naive ;
        size_t shift1_fast, shift2_fast ;   bool flip_fast ;  double d_fast ;

        // 2 vectors which need shifting
        for(size_t shift=1; shift<7; shift++)
        {   d_naive = cor_distance(    v3, v4, shift, shift1_naive, shift2_naive) ;
            d_fast = cor_distance_fast(v3, v4, shift, shift1_fast,  shift2_fast) ;
            CHECK_CLOSE(d_naive, d_fast, error) ;
            CHECK_EQUAL(shift1_naive, shift1_fast) ;
            CHECK_EQUAL(shift2_naive, shift2_fast) ;

            d_naive = cor_distance(    v3, v4, shift, shift1_naive, shift2_naive, flip_naive) ;
            d_fast = cor_distance_fast(v3, v4, shift, shift1_fast,  shift2_fast, flip_fast) ;

            CHECK_CLOSE(d_naive, d_fast, error) ;
            CHECK_EQUAL(shift1_naive, shift1_fast) ;
            CHECK_EQUAL(shift2_naive, shift2_fast) ;
            CHECK_EQUAL(flip_naive, flip_fast) ;
        }

        // 2 vectors out of which one has a sd==0
        for(size_t shift=1; shift<7; shift++)
        {   d_naive = cor_distance(    v9, v10, shift, shift1_naive, shift2_naive) ;
            d_fast = cor_distance_fast(v9, v10, shift, shift1_fast,  shift2_fast) ;
            CHECK_CLOSE(d_naive, d_fast, error) ;
            CHECK_EQUAL(shift1_naive, shift1_fast) ;
            CHECK_EQUAL(shift2_naive, shift2_fast) ;

            d_naive = cor_distance(    v9, v10, shift, shift1_naive, shift2_naive, flip_naive) ;
            d_fast = cor_distance_fast(v9, v10, shift, shift1_fast,  shift2_fast, flip_fast) ;
            CHECK_CLOSE(d_naive, d_fast, error) ;
            CHECK_EQUAL(shift1_naive, shift1_fast) ;
            CHECK_EQUAL(shift2_naive, shift2_fast) ;
            CHECK_EQUAL(flip_naive, flip_fast) ;
        }

        // the smallest possible vector for which it is possible to compute a correlation
        std::vector<double> v11 = {1, 2} ;
        std::vector<double> v12 = {2, 1} ;
        for(size_t shift=1; shift<1; shift++)
        {   d_naive = cor_distance(    v11, v12, shift, shift1_naive, shift2_naive) ;
            d_fast = cor_distance_fast(v11, v12, shift, shift1_fast,  shift2_fast) ;
            CHECK_CLOSE(d_naive, d_fast, error) ;
            CHECK_EQUAL(shift1_naive, shift1_fast) ;
            CHECK_EQUAL(shift2_naive, shift2_fast) ;

            d_naive = cor_distance(    v11, v12, shift, shift1_naive, shift2_naive, flip_naive) ;
            d_fast = cor_distance_fast(v11, v12, shift, shift1_fast,  shift2_fast, flip_fast) ;
            CHECK_CLOSE(d_naive, d_fast, error) ;
            CHECK_EQUAL(shift1_naive, shift1_fast) ;
            CHECK_EQUAL(shift2_naive, shift2_fast) ;
            CHECK_EQUAL(flip_naive, flip_fast) ;
        }

        // vectors with 1 element -> no correlation should be possible, distance should be
        // equal to 1 - Constants::corr_fail
        std::vector<double> v13 = {1} ;
        std::vector<double> v14 = {2} ;
        d_naive = cor_distance(    v13, v14, 1, shift1_naive, shift2_naive) ;
        d_fast = cor_distance_fast(v13, v14, 1, shift1_fast,  shift2_fast) ;
        CHECK_EQUAL(d_naive, 1. - Constants::corr_fail) ;
        CHECK_EQUAL(d_naive, d_fast) ;
        CHECK_EQUAL(shift1_naive, shift1_fast) ;
        CHECK_EQUAL(shift2_naive, shift2_fast) ;

        // stress test
        std::cout << "Running Distances stress tests..." << std::endl ;
        error = 10e-8 ;
        size_t n_tests = 1000000 ;
        for(size_t i=0; i<n_tests; i++)
        {   // vector size
            size_t l = (rand() % 200) + 2 ;
            // shift
            // l can be 1 -> l/2 becomes 0 -> modulo 0 is an error, make it 1
            size_t shift = (rand() % (l/2 == 0 ? 1 : l/2)) + 1 ;
            // vectors
            std::vector<double> v1(l) ; for(size_t j=0; j<l; j++) { v1[j] = rand() % 10000 ; }
            std::vector<double> v2(l) ; for(size_t j=0; j<l; j++) { v2[j] = rand() % 10000 ; }

            // no flip
            double d_naive = cor_distance(v1, v2, shift, shift1_naive, shift2_naive) ;
            double d_fast  = cor_distance_fast(v1, v2, shift, shift1_fast,  shift2_fast) ;
            if(not isEqual(d_naive, d_fast, error) or
               (shift1_naive != shift1_fast) or
               (shift2_naive != shift2_fast))
            {   std::cerr << "ERROR" << std::endl ;
                std::cerr << "length      " << l << std::endl ;
                std::cerr << "shift       " << shift << std::endl ;
                std::cerr << "dist naive  " << d_naive << std::endl ;
                std::cerr << "dist fast   " << d_fast << std::endl ;
                std::cerr << "shift naive " << shift1_naive << " " << shift2_naive << std::endl ;
                std::cerr << "shift fast  " << shift1_fast  << " " << shift2_fast  << std::endl ;
                std::cerr << v1 << std::endl ;
                std::cerr << v2 << std::endl ;
                throw std::runtime_error("failed assertion in stress test!") ;
            }

            // with flip
            d_naive = cor_distance(v1, v2, shift, shift1_naive, shift2_naive, flip_naive) ;
            d_fast  = cor_distance_fast(v1, v2, shift, shift1_fast,  shift2_fast, flip_fast) ;
            if(not isEqual(d_naive, d_fast, error) or
               (shift1_naive != shift1_fast) or
               (shift2_naive != shift2_fast) or
               (flip_naive   != flip_fast))
            {   std::cerr << "ERROR" << std::endl ;
                std::cerr << "length      " << l << std::endl ;
                std::cerr << "shift       " << shift << std::endl ;
                std::cerr << "dist naive  " << d_naive << std::endl ;
                std::cerr << "dist fast   " << d_fast << std::endl ;
                std::cerr << "shift naive " << shift1_naive << " " << shift2_naive << std::endl ;
                std::cerr << "shift fast  " << shift1_fast  << " " << shift2_fast  << std::endl ;
                std::cerr << "flip naive  " << flip_naive << std::endl ;
                std::cerr << "flip fast   " << flip_fast << std::endl ;
                std::cerr << v1 << std::endl ;
                std::cerr << v2 << std::endl ;
                throw std::runtime_error("failed assertion in stress test!") ;
            }
            if((i+1) % 100000 == 0)
            {   std::cout << i+1 << " tests run" << std::endl ; }
        }

        std::cout << std::endl ;
        std::cout << "Distances stress test done! " << std::endl ;
        std::cout << "tests run successfully : " << n_tests << std::endl << std::endl ;
    }
}



































