
#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP


/*!
 * \brief The Constants class contains miscellaneous constants
 * for the program.
 */
class Constants
{
   public:
    Constants() = delete ;
    Constants(Constants&) = delete ;

    // numerical constants
    static const double delta ;     // a small value to add to avoid zeroes
    static const double epsilon ;   // an epsilon value for double comparisons
    static const double corr_fail ; // a special correlation value to be used in case
                                    // the correlation computation failed.
                                    // if cor(x,y) cannot be computed, assume that x and
                                    // y are as correlated
                                    // correlated with each other as they are
                                    // anti-correlated

    // collections of values
    enum flip {FORWARD=0, REVERSE, N_FLIP_STATES=2} ;
    enum dist_type {COR_DIST=0, NORM_COR_DIST, N_DIST=2} ;
    enum clustering_codes {CONVERGENCE=0, SUCCESS, FAILURE, N_CODES=3} ;
} ;


#endif // CONSTANTS_HPP
