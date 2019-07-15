#/*============================================================================
#
#  PYVINECOPULIB: A python interface to vinecopulib.
#
#  Copyright (c) University College London (UCL). All rights reserved.
#
#  This software is distributed WITHOUT ANY WARRANTY; without even
#  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
#  PURPOSE.
#
#  See LICENSE.txt in the top level directory for details.
#
#============================================================================*/

set(DATA_DIR ${CMAKE_SOURCE_DIR}/Testing/Data)
set(TMP_DIR ${CMAKE_BINARY_DIR}/Testing/Temporary)

set(TEST_CASES
  pvBasicTest
  pvCommandLineArgsTest
)

foreach(_test_case ${TEST_CASES})
  add_executable(${_test_case} ${_test_case}.cpp pvCatchMain.cpp)
  target_link_libraries(${_test_case} ${ALL_LIBRARIES})
endforeach()

if (PYVINECOPULIB_USE_CUDA)
  set(CUDA_TEST_CASES
    mpCudaTest
  )
  foreach(_test_case ${CUDA_TEST_CASES})
    cuda_add_executable(${_test_case} ${_test_case}.cu pvCatchMain.cpp)
    target_link_libraries(${_test_case} ${ALL_LIBRARIES})
  endforeach()
endif()

# Adding tests is done outside the above loop,
# as each test will need different arguments.

##################################################################################################
# Dont forget its:  add_test(<test name (globally unique) > <exe name> <argument1> <argument2>
##################################################################################################

add_test(NoArgs ${EXECUTABLE_OUTPUT_PATH}/pvBasicTest)
add_test(1File ${EXECUTABLE_OUTPUT_PATH}/pvCommandLineArgsTest ${DATA_DIR}/input.txt)
if (PYVINECOPULIB_USE_CUDA)
  add_test(Cuda ${EXECUTABLE_OUTPUT_PATH}/mpCudaTest)
endif()
