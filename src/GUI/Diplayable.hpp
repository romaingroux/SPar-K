#ifndef DISPLAYABLE_HPP
#define DISPLAYABLE_HPP

/*!
 * \brief The Displayable class is an abstract class for classes meant
 * to be displayed on the screen.
 */
class Displayable
{
    public:
        /*!
         * \brief Destructor.
         */
        virtual ~Displayable() ;

         /*!
         * \brief Should trigger the display of the instance
         * on a given widget/window/stream.
         */
        virtual void display() const = 0;
} ;


#endif // DISPLAYABLE_HPP
