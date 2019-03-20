#ifndef APPLICATION_HPP
#define APPLICATION_HPP

#include <string>
#include <vector>

struct options
{   std::string file_data ;
    std::string file_references ;
    size_t thread_n ;
    size_t iteration_n ;
    size_t cluster_n ;
    size_t shift_n ;
    size_t width_n ;
    bool flip ;
    bool no_outlier ;
    bool debug ;
    std::string seed ;
    std::string seeding ;
    std::string distance ;
} ;


/*!
 * \brief The Application class is the interface to all
 * the program functionalities.
 * Mainly contains the run() method which runs the application.
 */
class Application
{
    public:

        Application() = delete ;
        Application(const Application& other) = delete ;

        /*!
         * \brief Constructor, initialize the object
         * running options from the command line given
         * arguments.
         * \param argv the vector of argument as fetch
         * from the command line.
         * \param argc the number of argument fetch from
         * the command line.
         */
        Application(int argn, char** argv) ;

        /*!
         * Destructor.
         */
        ~Application() ;

        /*!
         * \brief Runs the application and perform clustering.
         * \return the exit code, an error occured if non-0.
         */
        int run() ;

    private:
        /*!
         * \brief Sets this->options according to the options
         * given from the command line.
         * \param argn the argument corresponding to the argn argument of the main() function.
         * \param argv the argument corresponding to the argv argument of the main() function.
         */
        void setOptions(int argn, char** argv) ;

        /*!
         * \brief Checks whether debugging verbosity is on
         * (this->options.debug).
         * \return  whether debugging verbosity is on.
         */
        bool isDebugOn() const ;

        /*!
         * \brief Prints debug informations if the debug flag is on.
         * \param stream the stream to display on.
         */
        void printDebugInfo(std::ostream& stream) const ;

        // fields
        /*!
         * \brief stores the running options given from the command line.
         */
        struct options options ;

        /*!
         * \brief the application exit code. This value is returned by the
         * main() function to the operating system.
         * The values are :
         * -1 if an error occure while parsing the command line options
         * EXIT_SUCCESS if the clustering could be run (until convergence or not).
         * EXIT_FAILURE if the clustering failed.
         */
        int exit_code ;

} ;



#endif // APPLICATION_HPP
