#ifndef THREADPOOL_HPP
#define THREADPOOL_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include <thread>
#include <queue>
#include <mutex>
#include <string>


typedef void (*function_ptr)(size_t, size_t) ;


/*!
 * \brief The ThreadPool class implements a simple pool of working threads.
 * At construction, <n> threads are spawned and tries to get tasks to execute
 * from a queue. Any access to the queue is synchonized through the use of
 * a mutex.
 * Jobs are added to the queue by calling addJob(task) where task is the
 * result of std::bind(). The queue can be open - in which case the jobs are
 * effectively added to the queues - or closed - in which case nothing can
 * be push into the queue anymore.
 * Stopping the pool is done through calling join() which will close the queue
 * - and eventually let the threads empty it - and join all the threads.
 * Any access to the queue is synchonized through the use of a mutex.
 */
class ThreadPool
{
    public:
        /*!
         * \brief Default constructor.
         * Constructs a thread pool containing the
         * given number of threads, by default 1.
         * \param n_threads the number of threads, by default 1.
         * \param debug enables debugging verbosity.
         */
        ThreadPool(size_t n_threads=1, bool debug=false) ;

        ThreadPool(const ThreadPool& other) = delete ;

        /*!
          * \brief class destructor
          */
        ~ThreadPool() ;

        /*!
         * \brief gets the number of threads in the pool.
         * \return the number of threads.
         */
        size_t getNThread() const ;

        /*!
         * \brief adds a task in the queue for the threads.
         * Once join() has been called, the job queues are closed
         * and calling this method remains effectless.
         * \param task a function bound to its arguments using std::bind()
         * WARNING : I'dont know whether this is portable :-/
         */
        void addJob(std::function<void()>&& task) ;

        /*!
         * \brief closes the job queues and join the threads.
         * When this methods returns, all the jobs have been
         * run and each threads has been joined.
         */
        void join() ;

    private:
        /*!
         * \brief runs an inifinite routine during which
         * the caller tries to get a task to execute from
         * task queue.
         * The routine is interrupted when the queue is
         * closed AND empty.
         */
        void thread_routine() ;

        /*!
         * \brief locks the queue accession mutex.
         */
        void lock_mutex_queue() ;

        /*!
         * \brief unlocks the queue accession mutex.
         */
        void unlock_mutex_queue() ;

        /*!
         * \brief opens the jobs queue. Later calls to
         * addJob() have an effect.
         */
        void open_queue() ;

        /*!
         * \brief closes the jobs queue. Later calls to
         * addJob() remain effectless.
         */
        void close_queue() ;

        /*!
         * \brief checks whether the queue is open.
         * \return whether the job queue is open.
         */
        bool isQueueOpen() ;

        /*!
         * \brief checks whether the debugging verbosity is on.
         * \return whether the debugging verbosity is on.
         */
        bool isDebugOn() const ;

        /*!
         * \brief formats and prints the given debug message
         * (with the thread number) to the given stream.
         * The message is effectively printed only if this->debug
         * is set to true (set at construction time).
         * \param msg the message to print.
         */
        void debug_print(const std::string& msg, std::ostream& out=std::cerr) const ;

        // fields
        /*!
         * \brief the threads.
         */
        std::vector<std::thread> threads ;

        /*!
         * \brief the task queue.
         */
        std::queue<std::function<void()>> queue_task ;

        /*!
         * \brief the synchronization mutex to access
         * the queues.
         */
        std::mutex queue_mutex ;

        /*!
         * \brief whether the queues are open for pushing
         * or not.
         */
        bool queue_open ;

        /*!
         * \brief whether debugging verbosity is enabled.
         */
        bool debug ;

} ;


#endif // THREADPOOL_HPP
