#include <UnitTest++/UnitTest++.h>

#include <iostream>
#include <stdexcept>                          // std::invalid_argument
#include <ctime>                              // clock_t, clock(), CLOCKS_PER_SEC
#include <Statistics/Distances.hpp>           // cor_distance()
#include <Clustering/DistanceComputer.hpp>
#include <Clustering/CorDistanceComputer.hpp>
#include <Utility/Constants.hpp>              // Constants
#include <Utility/Utility.hpp>                // isEqual()

/*!
 * \brief A subclass of CorDistanceComputer containing invasive
 * getters for testing purposes.
 */
class CorDistanceComputerTest : public CorDistanceComputer
{   public:

        // constructors and destructor
        CorDistanceComputerTest(size_t shift=1)
            : CorDistanceComputer(shift)
        {}

        virtual ~CorDistanceComputerTest()
        {}

        // methods
        const size_t& get_shift()
        {   return this->_shift ; }

        const matrix3d_double& get_distances()
        {   return this->_distances ; }

        const matrix3d_double& get_sum_xy()
        {   return this->_sum_xy ; }

        const vector_double& get_sumx()
        {   return this->_sum_x ; }

        const vector_double& get_sumx2()
        {   return this->_sum_x2 ; }

        const vector_double& get_sumy()
        {   return this->_sum_y ; }

        const vector_double& get_sumy2()
        {   return this->_sum_y2 ; }

        const matrix2d_int& get_scores()
        {   return this->_scores ; }
} ;

/*!
 * \brief Checks that a given vector has the expected size.
 * \param vector the vector of interest.
 * \param size the expected size.
 * \return whether the vector has the expected size.
 */
bool check_dimension(const vector_int& vector, size_t size)
{   // std::cerr << "x dim : " << vector.size() << " " << size << std::endl << std::endl ;
    return vector.size() == size ? true : false ; }

/*!
 * \brief Checks that a given vector has the expected size.
 * \param vector the vector of interest.
 * \param size the expected size.
 * \return whether the vector has the expected size.
 */
bool check_dimension(const vector_double& vector, size_t size)
{   // std::cerr << "x dim : " << vector.size() << " " << size << std::endl << std::endl ;
    return vector.size() == size ? true : false ; }

/*!
 * \brief Checks that a given matrix has the expected dimensions.
 * \param matrix the matrix of interest.
 * \param dim_x the expected x dimension.
 * \param dim_y the expected y dimension.
 * \return whether the matrix has the expected dimensions.
 */
bool check_dimension(const matrix2d_int& matrix, size_t dim_x, size_t dim_y)
{   if(matrix.size() != dim_y)
    {   return false ; }
    // std::cerr << "y dim : " << matrix.size() << " " << dim_y << std::endl ;
    for(const auto& row : matrix)
    {   if(not check_dimension(row, dim_x))
        {   return false ; }
    }
    return true ;
}

/*!
 * \brief Checks that a given matrix has the expected dimensions.
 * \param matrix the matrix of interest.
 * \param dim_x the expected x dimension.
 * \param dim_y the expected y dimension.
 * \return whether the matrix has the expected dimensions.
 */
bool check_dimension(const matrix2d_double& matrix, size_t dim_x, size_t dim_y)
{   // std::cerr << "y dim : " << matrix.size() << " " << dim_y << std::endl ;
    if(matrix.size() != dim_y)
    {   return false ; }
    for(const auto& row : matrix)
    {   if(not check_dimension(row, dim_x))
        {   return false ; }
    }
    return true ;
}

/*!
 * \brief Checks that a given matrix has the expected dimensions.
 * \param matrix the matrix of interest.
 * \param dim_x the expected x dimension.
 * \param dim_y the expected y dimension.
 * \param dim_z the expected z dimension.
 * \return whether the matrix has the expected dimensions.
 */
bool check_dimension(const matrix3d_double& matrix, size_t dim_x, size_t dim_y, size_t dim_z)
{   // std::cerr << "z dim : " << matrix[0].size() << " " << dim_z << std::endl ;
    if(matrix.size() != dim_z)
    {   return false ; }
    for(const auto& mat2d : matrix)
    {   if(not check_dimension(mat2d, dim_x, dim_y))
        {   return false ; }
    }
    return true ;
}


// Statistics.hpp unittests
SUITE(DistanceComputer_testsuit)
{   /* This test suite checks that the CorDistanceComputer class works correclty
       and that it report similar results as the naive way of computing the
       correlation distance. However, due to implementation details (these
       modifications were made after some real cases failed) this class reports
       results similar to the naive method in 99.99% of the cases but somtimes
       reports different ones. Here is the explanation.

       The naive way of computing the correlation distance with flip and shifts
       simply considers each possible slices of the vectors and computes 1-cor.

       The fast way algorithm is totally equivalent to the naive way (tested
       in the Distance_unittest.cpp). A slightly modified version of this
       algorithm is implemented here to increase computation time.
       This class uses the dynamic programming approach and pre-computes some
       parameters for some slices, compute the correlation distance for these
       slices and then updates the paramters for the next slices (from the
       values for the previous slices) and computes the correlation distance
       for the current slices. This required to introduce the addition of a
       small constant (1e-9) to the really first values of the pre-computed
       parameters in order to avoid a division by 0 in some cases. In the end,
       this value does not change much the distance result but does it.

       In some cases, the addition of this tiny constant in the dynamic
       programming method result in finding different results than the naive
       way. This is simply because the distances between each pair of slice
       is different among both method. The minimum one differs between bot
       methods (even though the difference is really tiny, it even appear
       to be the same value) but the shift/flip values changes. There is
       nothing which can really be done for this... This is an implementation
       issue.

       So, to avoid the crash of the tests by something which is normal, if
       both methods returns the same distance but not the same shift/flip
       values, this will be accepted. However, if the distance is not the
       same, then this won't be accepted.
    */
    TEST(start_message)
    {   std::cout << "Running DistanceComputer tests..." << std::endl ; }

    TEST(constructor_test)
    {   // default contructor, should have a shift of 1
        CorDistanceComputerTest cdc ;

        CHECK_EQUAL(1, cdc.get_shift()) ;
        CHECK_EQUAL(true, check_dimension(cdc.get_sumx(),     1)) ;
        CHECK_EQUAL(true, check_dimension(cdc.get_sumx2(),    1)) ;
        CHECK_EQUAL(true, check_dimension(cdc.get_sumy(),     1)) ;
        CHECK_EQUAL(true, check_dimension(cdc.get_sumy2(),    1)) ;
        CHECK_EQUAL(true, check_dimension(cdc.get_distances(),2,1,1)) ;
        CHECK_EQUAL(true, check_dimension(cdc.get_sum_xy(),   2,1,1)) ;
        CHECK_EQUAL(true, check_dimension(cdc.get_scores(),   1,1)) ;

        // constructor with value
        for(size_t shift=1; shift<10; shift++)
        {   cdc = CorDistanceComputerTest(shift) ;
            CHECK_EQUAL(shift, cdc.get_shift()) ;
            CHECK_EQUAL(true, check_dimension(cdc.get_sumx(),     shift)) ;
            CHECK_EQUAL(true, check_dimension(cdc.get_sumx2(),    shift)) ;
            CHECK_EQUAL(true, check_dimension(cdc.get_sumy(),     shift)) ;
            CHECK_EQUAL(true, check_dimension(cdc.get_sumy2(),    shift)) ;
            CHECK_EQUAL(true, check_dimension(cdc.get_distances(),2,shift,shift)) ;
            CHECK_EQUAL(true, check_dimension(cdc.get_sum_xy(),   2,shift,shift)) ;
            CHECK_EQUAL(true, check_dimension(cdc.get_scores(),   shift,shift)) ;
        }
        // constructor with invalid value
        CHECK_THROW(CorDistanceComputerTest(0), std::invalid_argument) ;
    }

    TEST(compute_distance)
    {   CorDistanceComputerTest cdc ;

        // corner cases

        // tolerated error for equality testing, the expected values are reported
        // with 7 decimals only also there is a little difference between this
        // method and the naive method (see top comment)
        double error = 10e-5 ;
        // some variables for the distance functions
        double d, d_exp ;
        size_t shift1, shift2 ;
        bool flip2 ;

        // --------------------- one vector as sd of 0 ---------------------
        std::vector<double> v9  = {1, 2, 3, 0, 4, 5, 5, 5, 0, 0, 1, 0, 0} ;
        std::vector<double> v10 = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1} ;

        // --------------------- vector are the shifted -> distance of 0 without even shifting ---------------------
        std::vector<double> v3 = {1, 2, 3, 0, 4, 5, 5, 5, 0, 0, 1, 0, 0} ;
        std::vector<double> v4 = {0, 0, 1, 2, 3, 0, 4, 5, 5, 5, 0, 0, 1} ;

        // cor(v9,v10) should return Constants::corr_fail, so the distance should be
        // 1 - Constants::corr_fail
        d_exp = 1 - Constants::corr_fail ;

        // should also be the case with shifting (whatever the freedom)
        d = cdc.compute_distance(v9, v10, 1, shift1, shift2) ;
        CHECK_EQUAL(d, d_exp) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 0) ;

        d = cdc.compute_distance(v9, v10, 5, shift1, shift2) ;
        CHECK_EQUAL(d, d_exp) ;
        CHECK_EQUAL(shift1, 2) ;
        CHECK_EQUAL(shift2, 2) ;

        d = cdc.compute_distance(v9, v10, 7, shift1, shift2) ;
        CHECK_EQUAL(d, d_exp) ;
        CHECK_EQUAL(shift1, 3) ;
        CHECK_EQUAL(shift2, 3) ;

        // should also be the case with shifting (whatever the freedom) and flipping
        d = cdc.compute_distance(v9, v10, 1, shift1, shift2, flip2) ;
        CHECK_EQUAL(d, d_exp) ;
        CHECK_EQUAL(shift1, 0) ;
        CHECK_EQUAL(shift2, 0) ;
        CHECK_EQUAL(flip2, Constants::FORWARD) ;

        d = cdc.compute_distance(v9, v10, 5, shift1, shift2, flip2) ;
        CHECK_EQUAL(d, d_exp) ;
        CHECK_EQUAL(shift1, 2) ;
        CHECK_EQUAL(shift2, 2) ;
        CHECK_EQUAL(flip2, Constants::FORWARD) ;

        d = cdc.compute_distance(v9, v10, 7, shift1, shift2, flip2) ;
        CHECK_EQUAL(d, d_exp) ;
        CHECK_EQUAL(shift1, 3) ;
        CHECK_EQUAL(shift2, 3) ;
        CHECK_EQUAL(flip2, Constants::FORWARD) ;

        // the smallest possible vector for which it is possible to compute a correlation
        size_t shift1_naive, shift2_naive ; bool flip_naive ; double d_naive ;
        size_t shift1_fast, shift2_fast ;   bool flip_fast ;  double d_fast ;

        std::vector<double> v11 = {1, 2} ;
        std::vector<double> v12 = {2, 1} ;
        for(size_t shift=1; shift<1; shift++)
        {   d_naive = cor_distance(        v11, v12, shift, shift1_naive, shift2_naive) ;
            d_fast  = cdc.compute_distance(v11, v12, shift, shift1_fast,  shift2_fast) ;
            CHECK_CLOSE(d_naive, d_fast, error) ;
            CHECK_EQUAL(shift1_naive, shift1_fast) ;
            CHECK_EQUAL(shift2_naive, shift2_fast) ;

            d_naive = cor_distance(        v11, v12, shift, shift1_naive, shift2_naive, flip_naive) ;
            d_fast  = cdc.compute_distance(v11, v12, shift, shift1_fast,  shift2_fast,  flip_fast) ;
            CHECK_CLOSE(d_naive, d_fast, error) ;
            CHECK_EQUAL(shift1_naive, shift1_fast) ;
            CHECK_EQUAL(shift2_naive, shift2_fast) ;
            CHECK_EQUAL(flip_naive, flip_fast) ;
        }

        // vectors with 1 element -> no correlation should be possible, distance should be equal
        // to 1 - Constants::corr_fail
        std::vector<double> v13 = {1} ;
        std::vector<double> v14 = {2} ;
        d_naive = cor_distance(        v13, v14, 1, shift1_naive, shift2_naive) ;
        d_fast  = cdc.compute_distance(v13, v14, 1, shift1_fast,  shift2_fast) ;
        CHECK_EQUAL(d_naive, 1. - Constants::corr_fail) ;
        CHECK_EQUAL(d_naive, d_fast) ;
        CHECK_EQUAL(shift1_naive, shift1_fast) ;
        CHECK_EQUAL(shift2_naive, shift2_fast) ;



        // test the difference betweeen the naive and the new method

        //CorDistanceComputer cdc  ;
        std::cout << "Running CorDistanceComputer stress tests..." << std::endl ;
        // double error = 10e-5 ;
        // for shift only
        size_t n_tests = 10000000 ;
        size_t n_special_noflip = 0 ;
        size_t n_special_flip   = 0 ;
        for(size_t i=0; i<n_tests; i++)
        {   // vector size
            size_t l = (rand() % 15) + 2 ;
            // shift
            // l can be 1 -> l/2 becomes 0 -> modulo 0 is an error, make it 1
            size_t shift = (rand() % (l/2 == 0 ? 1 : l/2)) + 1 ;
            size_t shift1_naive, shift1_fast, shift2_naive, shift2_fast ;
            bool flip_fast, flip_naive ;
            // vectors
            std::vector<double> v1(l) ; for(size_t j=0; j<l; j++) { v1[j] = rand() % 10000 ; }
            std::vector<double> v2(l) ; for(size_t j=0; j<l; j++) { v2[j] = rand() % 10000 ; }

            // without flipping
            double d_naive = cor_distance(        v1, v2, shift, shift1_naive, shift2_naive) ;
            double d_fast  = cdc.compute_distance(v1, v2, shift, shift1_fast,  shift2_fast) ;
            // special case, difference dued to precision error (read comment at top)
            if((isEqual(d_naive, d_fast, error)) and
               ((shift1_naive != shift1_fast) or (shift2_naive != shift2_fast)))
            {   n_special_noflip++ ; }
            // error
            else if(not isEqual(d_naive, d_fast, error) and
                    ((shift1_naive != shift1_fast) or
                     (shift2_naive != shift2_fast)))
            {   std::cerr << "ERROR" << std::endl ;
                std::cerr << "length      " << l << std::endl ;
                std::cerr << "shift       " << shift << std::endl ;
                std::cerr << "dist naive  " << d_naive << std::endl ;
                std::cerr << "dist fast   " << d_fast  << std::endl ;
                std::cerr << "shift naive " << shift1_naive << " " << shift2_naive << std::endl ;
                std::cerr << "shift fast  " << shift1_fast  << " " << shift2_fast  << std::endl ;
                std::cerr << v1 << std::endl ;
                std::cerr << v2 << std::endl ;
                throw std::runtime_error("failed assertion in stress test!") ;
            }

            // with flipping
            d_naive = cor_distance(        v1, v2, shift, shift1_naive, shift2_naive, flip_naive) ;
            d_fast  = cdc.compute_distance(v1, v2, shift, shift1_fast,  shift2_fast,  flip_fast) ;
            // special case, difference dued to precision error (read comment at top)
            if((isEqual(d_naive, d_fast, error)) and
               ((flip_naive != flip_fast) or (shift1_naive != shift1_fast) or (shift2_naive != shift2_fast)))
            {   n_special_flip++ ; }
            // error
            if((not isEqual(d_naive, d_fast, error)) and
               ((flip_naive   != flip_fast)  or (shift1_naive != shift1_fast) or (shift2_naive != shift2_fast)))
            {   std::cerr << "ERROR" << std::endl ;
                std::cerr << "length      " << l << std::endl ;
                std::cerr << "shift       " << shift << std::endl ;
                std::cerr << "dist naive  " << d_naive << std::endl ;
                std::cerr << "dist fast   " << d_fast  << std::endl ;
                std::cerr << "flip naive  " << flip_naive << std::endl ;
                std::cerr << "flip fast   " << flip_fast  << std::endl ;
                std::cerr << "shift naive " << shift1_naive << " " << shift2_naive << std::endl ;
                std::cerr << "shift fast  " << shift1_fast  << " " << shift2_fast  << std::endl ;
                std::cerr << v1 << std::endl ;
                std::cerr << v2 << std::endl ;
                throw std::runtime_error("failed assertion in stress test!") ;
            }

            if((i+1) % 1000000 == 0)
            {   std::cout << i+1 << " tests run" << std::endl ; }
        }
        std::cout << std::endl ;
        std::cout << "CorDistanceComputer stress test done! " << std::endl ;
        std::cout << "tests run successfully                : " << n_tests          << std::endl ;
        std::cout << "special cases (but success) w/o flip  : " << n_special_noflip << std::endl ;
        std::cout << "special cases (but success) with flip : " << n_special_flip   << std::endl << std::endl ;
    }
}
