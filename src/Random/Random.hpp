#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <random>
#include <vector>
#include <assert.h>

#include "RandomNumberGenerator.hpp"


/*!
 * \brief Generates a random number from a
 * Bernouilli distribution of parameter p.
 * \param p the probability of success.
 * \return a random number.
 */
bool rand_bernoulli(double p) ;


/*!
 * \brief Generates n random number from a
 * Bernouilli distribution of parameter p.
 * Not faster than rand_bernoulli(double p)
 * \param p the probability of success.
 * \param n the number of values to sample.
 * \return a vector of n random numbers.
 */
std::vector<bool> rand_bernoulli(double p, size_t n) ;


/*!
 * \brief Generates a random number from a
 * Normal distribution of mean m and standard
 * deviation sd.
 * \param m the mean.
 * \param sd the standard deviation.
 * \return a random number.
 */
double rand_normal(double m, double sd) ;


/*!
 * \brief Generates n random numbers from a
 * Normal distribution of mean m and standard
 * deviation sd.
 * More efficient for sampling than
 * rand_normal(double m, double sd).
 * \param m the mean.
 * \param sd the standard deviation.
 * \param n the number of values to sample.
 * \return a vector of n random numbers.
 */
std::vector<double> rand_normal(double m, double sd, size_t n) ;


/*! Generates a real random number from a uniform
 * distribution comprised between min and max.
 * \param min the lower limit of the distribution.
 * \param max the upper limit of the distribution.
 * \return a random number.
 */
template<typename T>
T rand_real_uniform(T min, T max)
{   std::uniform_real_distribution<T> dist(min, max) ;
    return dist(getRandomGenerator()) ;
}


/*! Generates n real random numbers from a uniform
 * distribution comprised between min and max.
 * \param min the lower limit of the distribution.
 * \param max the upper limit of the distribution.
 * \param n the number of value to sample.
 * \return a vector of n random number.
 */
template<typename T>
std::vector<T> rand_real_uniform(T min, T max, size_t n)
{
    assert(n > 0) ;

    std::vector<T> vector(n) ;
    std::uniform_real_distribution<T> dist(min, max) ;

    for(size_t i=0; i<n; i++)
    {   vector[i] = dist(getRandomGenerator()) ; }
    return vector ;
}

/*! Generates a random integer from a uniform
 * distribution comprised between min and max.
 * \param min the lower limit of the distribution.
 * \param max the upper limit of the distribution.
 * \return a random number.
 */
template<typename T>
T rand_int_uniform(T min, T max)
{   std::uniform_int_distribution<T> dist(min, max) ;
    return dist(getRandomGenerator()) ;
}


/*! Generates n random integers from a uniform
 * distribution comprised between min and max.
 * \param min the lower limit of the distribution.
 * \param max the upper limit of the distribution.
 * \param n the number of value to sample.
 * \return a vector of n random number.
 */
template<typename T>
std::vector<T> rand_int_uniform(T min, T max, size_t n)
{
    assert(n > 0) ;

    std::vector<T> vector(n) ;
    std::uniform_int_distribution<T> dist(min, max) ;

    for(size_t i=0; i<n; i++)
    {   vector[i] = dist(getRandomGenerator()) ; }
    return vector ;
}


#endif //RANDOM_HPP
