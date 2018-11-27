#ifndef RANDOMNUMBERGENERATOR_HPP
#define RANDOMNUMBERGENERATOR_HPP

#include <string>
#include <random>


/*!
 * \brief Initialise and randomly seeds a 32-bit Mesrenne Twister random
 * number generator and returns it.
 * Initialisation and seeding are static thus from the first time this
 * function is called, it will ALWAYS return the same generator.
 * \param seed a value to seed the random generator with. If seed is 0,
 * then the generator is seeded using a random generator.
 * \return ALWAYS THE SAME random number generator.
 */
std::mt19937& getRandomGenerator(std::string seed="") ;


#endif
