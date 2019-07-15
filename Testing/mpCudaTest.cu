/*=============================================================================

  PYVINECOPULIB: A python interface to vinecopulib.

  Copyright (c) University College London (UCL). All rights reserved.

  This software is distributed WITHOUT ANY WARRANTY; without even
  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
  PURPOSE.

  See LICENSE.txt in the top level directory for details.

=============================================================================*/

#include "catch.hpp"
#include "pvCatchMain.h"
#include <iostream>
#include <thrust/version.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/copy.h>
#include <algorithm>
#include <cstdlib>
#include <time.h>

TEST_CASE( "Check Thrust Exists", "[CUDA]" ) {
  int major = THRUST_MAJOR_VERSION;
  int minor = THRUST_MINOR_VERSION;
  std::cout << "Thrust v" << major << "." << minor << std::endl;
  REQUIRE(major > 0);
  REQUIRE(minor > 0);
}

TEST_CASE( "Copy to device, sort", "[CUDA]" ) {
  thrust::host_vector<int> h(2);
  h[0] = 3;
  h[1] = 1;
  thrust::device_vector<int> d = h;
  thrust::sort(d.begin(), d.end());
  REQUIRE(d[0] == 1);
  REQUIRE(d[1] == 3);
}

TEST_CASE( "Sort 32M numbers on GPU", "[CUDA]" ) {

  // Create test array.
  thrust::host_vector<int> h(32 << 20);
  std::generate(h.begin(), h.end(), rand);

  // Start clock.
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // Copy to card, compute, copy back.
  thrust::device_vector<int> d = h;
  thrust::sort(d.begin(), d.end());
  thrust::copy(d.begin(), d.end(), h.begin());

  // Stop clock.
  std::chrono::steady_clock::time_point end= std::chrono::steady_clock::now();
  std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() <<std::endl;

  REQUIRE(h.size() == d.size());
}
