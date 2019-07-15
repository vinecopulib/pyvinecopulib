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

#########################################################################################
#
# Usage: mpMacroGetCommitHashOfCurrentFile(commit_hash_var)
#
# Retrieves the hash of the commit of the last modification of the CMake list file
# from which the macro is called.
# The macro stores the result in the 'commit_hash_var' variable.
#
#########################################################################################


macro(mpMacroGetCommitHashOfCurrentFile commit_hash_var)

  execute_process(COMMAND ${GIT_EXECUTABLE} log -n 1 --pretty=format:%h -- ${CMAKE_CURRENT_LIST_FILE}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    ERROR_VARIABLE GIT_error
    OUTPUT_VARIABLE ${commit_hash_var}
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(NOT ${GIT_error} EQUAL 0)
    message(SEND_ERROR "Command \"${GIT_EXECUTBALE} log -n 1 --pretty=format:\"%h\" -- ${CMAKE_CURRENT_LIST_FILE} in directory ${CMAKE_SOURCE_DIR} failed with output:\n${GIT_error}")
  endif()

endmacro()


#########################################################################################
#
# Usage: mpMacroDefineExternalProjectVariables(project_name version_number)
#
# Defines variables that are needed to set up an external project.
# The proj_DEPENDENCIES variable is set to an empty list. If the project depends
# on other external projects, it needs to be updated after the call of this macro.
#
#########################################################################################

macro(mpMacroDefineExternalProjectVariables project version location)

  if (EP_DIRECTORY_PER_VERSION)

    mpMacroGetCommitHashOfCurrentFile(config_version)
    string(SUBSTRING ${config_version} 0 5 config_version)

    set(version_subdir "/${config_version}")

  else()

    set(version_subdir "")

  endif()

  set(proj ${project})
  set(proj_VERSION ${version})
  set(proj_LOCATION ${location})
  set(proj_DIR ${EP_BASE}/${proj}${version_subdir})
  set(proj_CONFIG ${proj_DIR}/cmake)
  set(proj_SOURCE ${proj_DIR}/src)
  set(proj_BUILD ${proj_DIR}/build)
  set(proj_INSTALL ${proj_DIR}/install)
  set(proj_DEPENDENCIES "")
  set(${project}_DEPENDS ${project})

  if (${location} MATCHES "^.*(\\.tar\\.gz|\\.tar\\.bz2)$")
    mpMacroGetChecksum(proj_CHECKSUM ${proj_LOCATION})
  endif()

  set(${project}_VERSION ${version})
  set(${project}_LOCATION ${location})

endmacro()


#########################################################################################
#
# Usage: mpMacroGetChecksum(RESULT_VAR FILE_URI)
#
# Downloads the md5 checksum file for the file and stores the checksum
# in RESULT_VAR. It expects that the checksum file has the same name as
# the original file plus the '.md5' extension.
#
#########################################################################################

macro(mpMacroGetChecksum RESULT_VAR FILE_URI)

  # We expect that the checksum has the name of the original file plus
  # the '.md5' extension.
  set(MD5_FILE_URI "${FILE_URI}.md5")

  # Cuts the host name and directory and keeps the file name only:
  string(REGEX REPLACE ".*/" "" MD5_FILE ${MD5_FILE_URI})

  # Downloads the md5 file:
  file(DOWNLOAD "${MD5_FILE_URI}" "${proj_CONFIG}/src/${MD5_FILE}")

  # Reads the first 32B to the output variable. (MD5 checksums are 128b.)
  file(STRINGS "${proj_CONFIG}/src/${MD5_FILE}" checksum LIMIT_INPUT 32)

  set(${RESULT_VAR} ${checksum})
endmacro()
