#ifndef UPDATABLE_HPP
#define UPDATABLE_HPP

/*!
 * \brief The Updatable class is an abstract class for classes meant
 * to be updated over time.
 */
class Updatable
{
    public:
         /*!
         * \brief Destructor.
         */
        virtual ~Updatable() ;

         /*!
         * \brief Should trigger the update of the instance.
         */
        virtual void update() = 0;
} ;


#endif // UPDATABLE_HPP
