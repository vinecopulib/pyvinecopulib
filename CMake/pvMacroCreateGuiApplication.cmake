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

macro(pvMacroCreateGuiApplication APP_NAME ADDITIONAL_SEARCH_PATHS)

  # Installs the main application.

  set(plugin_dest_dir bin)
  set(APPS "\${CMAKE_INSTALL_PREFIX}/bin/${APP_NAME}")
  if(APPLE)
    set(plugin_dest_dir ${APP_NAME}.app/Contents/MacOS)

    set(APPS "\${CMAKE_INSTALL_PREFIX}/${APP_NAME}.app")
  endif()
  if(WIN32)
    set(APPS "\${CMAKE_INSTALL_PREFIX}/bin/${APP_NAME}.exe")
  endif()
  install(TARGETS ${APP_NAME}
          BUNDLE DESTINATION . COMPONENT Runtime
          RUNTIME DESTINATION bin COMPONENT Runtime
         )

  # Writes out qt.conf which sets up the paths.

  set(qtconf_dest_dir bin)
  if(APPLE)
    set(qtconf_dest_dir ${APP_NAME}.app/Contents/Resources)
  endif()
  set(_qt_conf_plugin_install_prefix "Prefix=.")
  if(APPLE)
    set(_qt_conf_plugin_install_prefix "Prefix=./MacOS")
  endif()
  set(_qt_conf_lib_prefix)
  if(NOT APPLE)
    set(_qt_conf_lib_prefix "Libraries=.")
  endif()
  install(CODE "
          file(WRITE \"\${CMAKE_INSTALL_PREFIX}/${qtconf_dest_dir}/qt.conf\" \"
  [Paths]
  ${_qt_conf_plugin_install_prefix}
  ${_qt_conf_lib_prefix}
  \")" COMPONENT Runtime)

  # Setup the Icon for Macs.

  if(APPLE)
    set_target_properties(${APP_NAME} PROPERTIES MACOSX_BUNDLE_NAME "${APP_NAME}")
    set(icon_name "icon.icns")
    set(icon_full_path "${CMAKE_CURRENT_SOURCE_DIR}/icons/${icon_name}")
    if(EXISTS "${icon_full_path}")
      set_target_properties(${APP_NAME} PROPERTIES MACOSX_BUNDLE_ICON_FILE "${icon_name}")
      file(COPY ${icon_full_path} DESTINATION "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${APP_NAME}.app/Contents/Resources/")
      install(FILES ${icon_full_path} DESTINATION "${APP_NAME}.app/Contents/Resources/")
    endif()
  endif()

  # Installs ALL Qt plugins.

  if(Qt5_DIR)
    get_property(_qmake_location TARGET ${Qt5Core_QMAKE_EXECUTABLE}
                 PROPERTY IMPORT_LOCATION)
    get_filename_component(_qmake_path "${_qmake_location}" DIRECTORY)

    set(_plugin_source_dir ${_qmake_path}/../plugins)
    set(_plugin_dest_dir bin/)
    if(APPLE)
      set(_plugin_dest_dir ${APP_NAME}.app/Contents/MacOS/)
    endif()
    file(GLOB children RELATIVE ${_plugin_source_dir} ${_plugin_source_dir}/*)
    foreach( c ${children})
      if( IS_DIRECTORY ${_plugin_source_dir}/${c} )
        install(DIRECTORY ${_plugin_source_dir}/${c} DESTINATION ${_plugin_dest_dir})
      else()
        install(FILES ${_plugin_source_dir}/${c} DESTINATION ${_plugin_dest_dir})
      endif()
    endforeach()
  endif()

  # Calls fixup_bundle.

  install(CODE "
          file(GLOB_RECURSE QTPLUGINS
          \"\${CMAKE_INSTALL_PREFIX}/${plugin_dest_dir}/platforms/*${CMAKE_SHARED_LIBRARY_SUFFIX}\")
          file(GLOB_RECURSE QMLPLUGINS
          \"\${CMAKE_INSTALL_PREFIX}/${plugin_dest_dir}/qml/*${CMAKE_SHARED_LIBRARY_SUFFIX}\")
          list(APPEND QTPLUGINS \"\${QMLPLUGINS}\")
          include(BundleUtilities)
          fixup_bundle(\"${APPS}\" \"\${QTPLUGINS}\" \"${ADDITIONAL_SEARCH_PATHS}\")
          " COMPONENT Runtime)

endmacro()
