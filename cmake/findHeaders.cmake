set(vinecopulib_includes ${CMAKE_SOURCE_DIR}/include)

if(VINECOPULIB_SHARED_LIB)
    set(vinecopulib_generated_sources ${CMAKE_BINARY_DIR}/generated/src)
    set(vinecopulib_generated_includes ${CMAKE_BINARY_DIR}/generated/include)

    string(CONCAT license "// Copyright © 2016-2019 Thomas Nagler and Thibault Vatter \n"
            "// \n"
            "// This file is part of the vinecopulib library and licensed under the terms of \n"
            "// the MIT license. For a copy, see the LICENSE file in the root directory of \n"
            "// vinecopulib or https://vinecopulib.github.io/vinecopulib/. \n")

    file(GLOB_RECURSE vinecopulib_ipp ${vinecopulib_includes}/*.ipp)
    foreach (file ${vinecopulib_ipp})

        # Get directory, name and path for header/source files
        get_filename_component(name_without_extension ${file} NAME_WE)
        get_filename_component(directory ${file} DIRECTORY)
        string(REGEX REPLACE "/implementation" "" directory ${directory})
        set(header_file ${directory}/${name_without_extension}.hpp)
        string(REGEX REPLACE "${vinecopulib_includes}/" ""
                header_file ${header_file})
        string(REGEX REPLACE ${vinecopulib_includes}/vinecopulib
                ${vinecopulib_generated_sources} source_folder ${directory})
        set(source_file "${source_folder}/${name_without_extension}.cpp")

        # Scrap file content, remove inline & add include
        file(READ ${file} file_content)
        string(REGEX REPLACE "inline " "" file_content "${file_content}")
        string(SUBSTRING "${file_content}" 282 -1 file_content)
        string(CONCAT file_content "${license}"
                "\n#include <${header_file}>" "${file_content}")

        # If source does not exists or has changed, generate new source file
        if(EXISTS "${source_file}")
            file(READ ${source_file} old_content)
            string(COMPARE NOTEQUAL "${file_content}" "${old_content}" has_changed)
            if(has_changed)
                file(WRITE ${source_file} "${file_content}")
            endif()
        else()
            file(WRITE ${source_file} "${file_content}")
        endif()

    endforeach ()

    file(GLOB_RECURSE vinecopulib_hpp ${vinecopulib_includes}/vinecopulib/*.hpp)
    foreach (file ${vinecopulib_hpp})
        # Get directory, name and path for header/source files
        get_filename_component(name_without_extension ${file} NAME_WE)
        get_filename_component(directory ${file} DIRECTORY)
        set(ipp_file ${directory}/implementation/${name_without_extension}.ipp)
        string(REGEX REPLACE "${vinecopulib_includes}/" "" ipp_file ${ipp_file})
        string(REGEX REPLACE ${vinecopulib_includes} ${vinecopulib_generated_includes}
                header_folder ${directory})
        set(header_file "${header_folder}/${name_without_extension}.hpp")

        # File scrap content and remove ipp include
        file(READ ${file} file_content)
        string(REGEX REPLACE "#include <${ipp_file}>" ""
                file_content "${file_content}")

        # If header does not exists or has changed, generate new header file
        if(EXISTS "${header_file}")
            file(READ ${header_file} old_content)
            string(COMPARE NOTEQUAL "${file_content}" "${old_content}" has_changed)
            if(has_changed)
                file(WRITE ${header_file} "${file_content}")
            endif()
        else()
            file(WRITE ${header_file} "${file_content}")
        endif()
    endforeach ()

    file(GLOB_RECURSE vinecopulib_main_hpp ${vinecopulib_includes}/vinecopulib.hpp)
    file(GLOB_RECURSE vinecopulib_version_hpp ${vinecopulib_includes}/version.hpp)
    file(COPY ${vinecopulib_main_hpp} DESTINATION ${vinecopulib_generated_includes})
    file(COPY ${vinecopulib_version_hpp} DESTINATION ${vinecopulib_generated_includes})

    set(vinecopulib_includes ${vinecopulib_generated_includes})
    file(GLOB_RECURSE vinecopulib_sources ${vinecopulib_generated_sources}/*.cpp)
else()
    file(GLOB_RECURSE vinecopulib_bicop_ipp
            ${vinecopulib_includes}/vinecopulib/bicop/implementation/*.ipp)
    file(GLOB_RECURSE vinecopulib_vinecop_ipp
            ${vinecopulib_includes}/vinecopulib/vinecop/implementation/*.ipp)
    file(GLOB_RECURSE vinecopulib_misc_ipp
            ${vinecopulib_includes}/vinecopulib/misc/implementation/*.ipp)
endif()

file(GLOB_RECURSE bicop_hpp ${vinecopulib_includes}/vinecopulib/bicop/*.hpp)
file(GLOB_RECURSE vinecop_hpp ${vinecopulib_includes}/vinecopulib/vinecop/*.hpp)
file(GLOB_RECURSE misc_hpp ${vinecopulib_includes}/vinecopulib/misc/*.hpp)
file(GLOB_RECURSE main_hpp ${vinecopulib_includes}/vinecopulib.hpp)
file(GLOB_RECURSE version_hpp ${vinecopulib_includes}/version.hpp)
