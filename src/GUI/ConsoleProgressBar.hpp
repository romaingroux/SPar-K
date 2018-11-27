#ifndef CONSOLEPROGRESSBAR_HPP
#define CONSOLEPROGRESSBAR_HPP

#include <iostream>
#include <string>

#include "Diplayable.hpp"
#include "Updatable.hpp"

/*!
 * \brief The ConsoleProgressBar class displays a progress bar on a given
 * stream to illustrate the progress of a process.
 * Exemple : ConsoleProgrssBar(cout, 100, 10, "sending")
 *           after calling update() 30x should display
 *           on std::cout :
 * sending : progress [===.......] 30.00 %
 *
 * 30% because 30x/100
 *
 */
class ConsoleProgressBar : public Displayable, public Updatable
{
    public :
        ConsoleProgressBar() = delete ;
        /*!
         * \brief Constructor. Constructs an empty bar.
         * \param stream the stream to display the bar on.
         * \param repeats the maximal number of times a piece of the process needs to be
         * repeated for the whole process to be done (and the progress bar to reach 100%).
         * \param size the length of the progress bar, in characters on the stream (not
         * accounting for some text on both sides.
         * \param prefix a short piece of text describing what the progress bar is
         * tracking, for information purposes only.
         */
        ConsoleProgressBar(std::ostream& stream, size_t repeats, size_t size, const std::string& prefix) ;

        /*!
         * \brief produces a display representation of the bar and sends it to stream.
         */
        virtual void display() const override ;
        /*!
         * \brief calls display() and updates the state of the bar by 1 (this method
         * can now be called one time less than before the call before is reached
         * 100%).
         */
        virtual void update() override ;

    private:
        /*!
         * \brief the number of times update() should be called before showing 100%.
         */
        size_t _repeats ;
        /*!
         * \brief the size of the bar in number of characters when reaching 100%.
         */
        size_t _size ;
        /*!
         * \brief the number of times update() was called so far.
         */
        size_t _current ;
        /*!
         * \brief a brief description of what the bar tracks, to be displayed.
         */
        std::string _prefix ;
        /*!
         * \brief a reference to the stream to display the bar representation to.
         */
        std::ostream& _stream ;
} ;


#endif // CONSOLEPROGRESSBAR_HPP
