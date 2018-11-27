#include <cmath>
#include <UnitTest++/UnitTest++.h>

#include "Utility/Utility.hpp"
#include "Matrix/Matrix2D.hpp"
#include "Utility/Constants.hpp"

using namespace std ;


SUITE(Utility_testsuit)
{
    // a message indicating the tests started
    TEST(start_message)
    {   cout << "Running Utility unittests..." << endl ; }

    //test seq() function
    TEST(seq_test)
    {   vector<int> expected1 = {-5,-4,-3,-2,-1,0,1,2,3,4,5} ;
        vector<int> results1  = seq(-5,5,1) ;
        CHECK_ARRAY_EQUAL(results1, expected1, expected1.size()) ;

        vector<int> expected2 = {-5,-4,-3,-2,-1,0} ;
        vector<int> results2  = seq(-5,0,1) ;
        CHECK_ARRAY_EQUAL(results2, expected2, expected2.size()) ;

        vector<int> expected3 = {-57,-47,-37,-27,-17,-7,3} ;
        vector<int> results3  = seq(-57,7,10) ;
        CHECK_ARRAY_EQUAL(results3, expected3, expected3.size()) ;
    }

    // tests pad_vector() function
    TEST(pad_vector_test)
    {   vector<int> expected1 = {1,2,3,0,0,0,0,0,0} ;
        vector<int> results1  = {0,0,0,1,2,3,0,0,0} ;
        pad_vector(results1, -3) ;
        CHECK_ARRAY_EQUAL(results1, expected1, expected1.size()) ;

        vector<int> expected2 = {0,0,0,0,0,0,0,0,0} ;
        pad_vector(results1, -3) ;
        CHECK_ARRAY_EQUAL(results1, expected2, expected1.size()) ;

        vector<int> expected3 = {0,0,0,0,0,1,2,3,0} ;
        vector<int> results3  = {0,0,0,1,2,3,0,0,0} ;
        pad_vector(results3, 2) ;
        CHECK_ARRAY_EQUAL(results3, expected3, expected3.size()) ;

        vector<int> expected4 = {0,0,0,0,0,0,0,1,2} ;
        pad_vector(results3, 2) ;
        CHECK_ARRAY_EQUAL(results3, expected4, expected4.size()) ;
    }

    // tests find_min() function for 2D matrix
    TEST(find_min_2d)
    {   Matrix2D<double> m(3,4,0) ;
        m(0,2) = 4. ;
        m(2,0) = -9. ;
        m(2,2) = 7. ;

        vector<size_t> expected1 = {2,0} ;
        vector<size_t> results1   = find_min(m) ;
        CHECK_ARRAY_EQUAL(results1, expected1, expected1.size()) ;

        // if there are two equally low values in the matrix
        // the function should keep the 1st one found
        m(2,3) = -9. ;
        vector<size_t> results2   = find_min(m) ;
        CHECK_ARRAY_EQUAL(results2, expected1, expected1.size()) ;

        m(1,1) = -18. ;
        vector<size_t> expected2 = {1,1} ;
        vector<size_t> results3   = find_min(m) ;
        CHECK_ARRAY_EQUAL(results3, expected2, expected2.size()) ;
    }

    // tests find_max() function for 2D matrix
    TEST(find_max_2d)
    {   Matrix2D<double> m(3,4,0) ;
        m(0,2) = 4. ;
        m(2,0) = 9. ;
        m(2,2) = 7. ;

        vector<size_t> expected1 = {2,0} ;
        vector<size_t> results1   = find_max(m) ;
        CHECK_ARRAY_EQUAL(results1, expected1, expected1.size()) ;

        // if there are two equally high values in the matrix
        // the function should keep the 1st one found
        m(2,3) = 9. ;
        vector<size_t> results2   = find_max(m) ;
        CHECK_ARRAY_EQUAL(results2, expected1, expected1.size()) ;

        m(1,1) = 18. ;
        vector<size_t> expected2 = {1,1} ;
        vector<size_t> results3   = find_max(m) ;
        CHECK_ARRAY_EQUAL(results3, expected2, expected2.size()) ;

    }

    // tests isNaN() function
    TEST(nan)
    {   double x = nan("") ;
        bool results = isNaN(x) ;
        bool expected = true ;
        CHECK_EQUAL(results, expected) ;
    }
}
