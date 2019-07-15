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

OPTION(PYVINECOPULIB_USE_CPPCHECK "Enable the use of CppCheck for static analysis." OFF)
IF(PYVINECOPULIB_USE_CPPCHECK)

  # Set the required CppCheck version
  SET(CPPCHECK_REQ_MAJOR 1)
  SET(CPPCHECK_REQ_MINOR 68)

  FIND_PROGRAM(CPPCHECK_EXECUTABLE
    NAMES cppcheck 
    PATHS
    /usr/local/bin
  )

  IF(CPPCHECK_EXECUTABLE)
  
    EXECUTE_PROCESS(
      COMMAND ${CPPCHECK_EXECUTABLE} --version
      OUTPUT_VARIABLE CPPCHECK_VERSION_TEXT
    ) 

    string(STRIP ${CPPCHECK_VERSION_TEXT} CPPCHECK_VERSION_TEXT)
    string(LENGTH ${CPPCHECK_VERSION_TEXT} CPPCHECK_VERSION_LENGTH)
    math(EXPR CPPCHECK_VERSION_FINAL_LENGTH "${CPPCHECK_VERSION_LENGTH}-9")
    string(SUBSTRING ${CPPCHECK_VERSION_TEXT} 9 ${CPPCHECK_VERSION_FINAL_LENGTH} CPPCHECK_VERSION)
    
    # now parse the parts of the user given version string into variables
    STRING(REGEX REPLACE "^([0-9]+)\\.[0-9]" "\\1" CPPCHECK_MAJOR_VERSION "${CPPCHECK_VERSION}")
    STRING(REGEX REPLACE "^[0-9]+\\.([0-9])" "\\1" CPPCHECK_MINOR_VERSION "${CPPCHECK_VERSION}")

    MATH(EXPR CPPCHECK_REQ_VERSION "${CPPCHECK_REQ_MAJOR}*10000 + ${CPPCHECK_REQ_MINOR}*100")
    MATH(EXPR CPPCHECK_LONG_VERSION "${CPPCHECK_MAJOR_VERSION}*10000 + ${CPPCHECK_MINOR_VERSION}*100")
  
    # Set the minimum require version for batchmake
    IF(CPPCHECK_LONG_VERSION LESS CPPCHECK_REQ_VERSION)
      MESSAGE(FATAL_ERROR "This project requires a newer version of CppCheck. Please upgrade the CppCheck executable.")
    ELSE(CPPCHECK_LONG_VERSION LESS CPPCHECK_REQ_VERSION)
      SET(CPPCHECK_FOUND 1)
    ENDIF(CPPCHECK_LONG_VERSION LESS CPPCHECK_REQ_VERSION)

    IF(CPPCHECK_FOUND)
      #
      #  Define file names
      #
      SET(CPPCHECK_PYVINECOPULIB_FILES_LIST
        ${PROJECT_BINARY_DIR}/Utilities/CppCheck/PYVINECOPULIBFiles.txt)

      #
      # Configure the files
      #
      CONFIGURE_FILE(
        ${PROJECT_SOURCE_DIR}/Utilities/CppCheck/PYVINECOPULIBFiles.txt.in
        ${CPPCHECK_PYVINECOPULIB_FILES_LIST})

      SET(CPPCHECK_ARGUMENTS_CODE
        --file-list=${CPPCHECK_PYVINECOPULIB_FILES_LIST}
        --enable=warning,style,performance,portability,missingInclude
        --error-exitcode=2
        --quiet
      )

      ADD_CUSTOM_COMMAND(
        OUTPUT  ${PYVINECOPULIB_BINARY_DIR}/CppCheckCodeReport.txt
        COMMAND ${CPPCHECK_EXECUTABLE}
        ARGS    ${CPPCHECK_ARGUMENTS_CODE}
        COMMENT "Coding Style Checker"
      )

      ADD_CUSTOM_TARGET(StaticallyAnalyzeCode DEPENDS ${PYVINECOPULIB_BINARY_DIR}/CppCheckCodeReport.txt)

      ADD_TEST(CppCheck ${CPPCHECK_EXECUTABLE} ${CPPCHECK_ARGUMENTS_CODE})

    ENDIF(CPPCHECK_FOUND)
  ENDIF(CPPCHECK_EXECUTABLE)
ENDIF(PYVINECOPULIB_USE_CPPCHECK)
