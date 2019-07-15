// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

#include <thread>
#include <vector>
#include <queue>
#include <future>
#include <condition_variable>
#include <atomic>

namespace vinecopulib {

namespace tools_thread {

//! Implementation of the thread pool pattern based on `std::thread`.
class ThreadPool {
public:

    ThreadPool() = default;
    ThreadPool(ThreadPool&&) = delete;
    ThreadPool(const ThreadPool&) = delete;

    //! constructs a thread pool with `nThreads` threads.
    //! @param nThreads number of threads to create; 0 means that all work
    //!   pushed to the pool will be done sequentially in the main thread.
    ThreadPool(size_t nThreads) : num_busy_(0), stopped_(false)
    {
        for (size_t t = 0; t < nThreads; t++) {
            pool_.emplace_back([this] {
                // observe thread pool as long there are jobs or pool has been
                // stopped
                while (!stopped_ | !jobs_.empty()) {
                    // thread must hold the lock while modifying shared
                    // variables
                    std::unique_lock<std::mutex> lk(m_);

                    // wait for new job or stop signal
                    cv_tasks_.wait(lk, [this] {
                        return stopped_ || !jobs_.empty();
                    });

                    // check if there are any jobs left in the queue
                    if (jobs_.empty())
                        continue;

                    // take job from the queue and signal that worker is busy
                    std::function<void()> job = std::move(jobs_.front());
                    jobs_.pop();
                    num_busy_++;
                    lk.unlock();

                    // execute job
                    job();

                    // signal that job is done
                    lk.lock();
                    num_busy_--;
                    cv_busy_.notify_one();
                }
            }
            );
        }
    }

    ~ThreadPool()
    {
        join();
    }

    // assignment operators
    ThreadPool& operator=(const ThreadPool&) = delete;
    ThreadPool& operator=(ThreadPool&& other) = default;

    //! pushes new jobs to the thread pool.
    //! @param f a function taking an arbitrary number of arguments.
    //! @param args a comma-seperated list of the other arguments that shall
    //!   be passed to `f`.
    //! @return an `std::shared_future`, where the user can get the result and
    //!   rethrow the catched exceptions.
    template<class F, class... Args>
    auto push(F&& f, Args&&... args) -> std::future<decltype(f(args...))>
    {
        // create packaged task on the heap to avoid stack overflows.
        auto job = std::make_shared<std::packaged_task<decltype(f(args...))()>>(
            [&f, args...] { return f(args...); }
        );

        // if there are no workers, just do the job in the main thread
        if (pool_.size() == 0) {
            (*job)();
            return job->get_future();
        }

        // add job to the queue
        {
            std::unique_lock<std::mutex> lk(m_);
            if (stopped_)
                throw std::runtime_error("cannot push to stopped thread pool");
            jobs_.emplace([job] () { (*job)(); });
        }

        // signal a waiting worker that there's a new job
        cv_tasks_.notify_one();

        // return future result of the job
        return job->get_future();
    }

    //! waits for all jobs to finish, but does not join the threads.
    void wait()
    {
        std::unique_lock<std::mutex> lk(m_);
        cv_busy_.wait(lk, [&] {return (num_busy_ == 0) && jobs_.empty();});
    }

    //! waits for all jobs to finish and joins all threads.
    void join()
    {
        // signal all threads to stop
        {
            std::unique_lock<std::mutex> lk(m_);
            stopped_ = true;
        }
        cv_tasks_.notify_all();

        // join threads if not done already
        if (pool_.size() > 0) {
            if (pool_[0].joinable()) {
                for (auto &worker : pool_)
                    worker.join();
            }
        }
    }

    //! maps a function on a list of items, possibly running tasks in parallel.
    //! @param f function to be mapped.
    //! @param items an objects containing the items on which `f` shall be mapped;
    //!     must allow for `auto` loops (i.e., `std::begin(I)`/`std::end(I)` must be
    //!     defined).
    template<class F, class I>
    void map(F&& f, I &&items)
    {
        for (const auto &item : items)
            push(f, item);
    }

private:
    std::vector<std::thread> pool_;           // worker threads in the pool
    std::queue<std::function<void()>> jobs_;  // the task queue

    // variables for synchronization between workers
    std::mutex m_;
    std::condition_variable cv_tasks_;
    std::condition_variable cv_busy_;
    std::atomic_uint num_busy_;
    bool stopped_;
};

}

}
