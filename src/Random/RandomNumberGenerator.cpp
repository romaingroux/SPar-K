#include "RandomNumberGenerator.hpp"


std::mt19937& getRandomGenerator(std::string seed)
{
    // initialise the generator once and for all
    static bool initialised = false ;
    static std::mt19937 generator ;

    // seeds the generator once and for all
    if(not initialised)
    {   if(seed == std::string(""))
        {   std::random_device rd ;
            generator.seed(rd()) ;
        }
        else
        {   std::seed_seq seed_sequence(seed.begin(),seed.end()) ;
            generator.seed(seed_sequence) ;
        }
        initialised = true ;
    }
    return generator ;
}
