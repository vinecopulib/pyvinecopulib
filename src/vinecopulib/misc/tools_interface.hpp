// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter
//
// This file is part of the vinecopulib library and licensed under the terms of
// the MIT license. For a copy, see the LICENSE file in the root directory of
// vinecopulib or https://vinecopulib.github.io/vinecopulib/.

#pragma once

// interface specfifc #defines can be set here
// (R package does: #define INTERFACED_FROM_R)

// interface specific headers
#ifdef INTERFACED_FROM_R
#include <RcppThread.h>
#define cout Rcout
namespace std {
static RcppThread::RPrinter Rcout = RcppThread::RPrinter();
}
#else
#include <iostream>
#endif

// parallel backend
#include <vinecopulib/misc/tools_batch.hpp>
#ifdef INTERFACED_FROM_R
namespace vinecopulib {
namespace tools_thread {
typedef RcppThread::ThreadPool ThreadPool;
}
}
#else
#include <vinecopulib/misc/tools_thread.hpp>
#endif

namespace vinecopulib {

namespace tools_interface {

inline void
check_user_interrupt(bool do_check = true)
{
  if (do_check) {
#ifdef INTERFACED_FROM_R
    RcppThread::checkUserInterrupt();
#endif
  }
}
}
}
