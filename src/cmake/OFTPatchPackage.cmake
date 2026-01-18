if(APPLE)
    install(
      CODE
      "
      function(gp_item_default_embedded_path_override item path)
        set(path @executable_path PARENT_SCOPE)
      endfunction()
      include(BundleUtilities)
      set(libs \"\")
      foreach(lib ${OFT_SHARED_LIBS})
        set(libs \${libs} \"\${CMAKE_INSTALL_PREFIX}/bin/lib\${lib}.dylib\")
      endforeach()
      fixup_bundle(\"\${CMAKE_INSTALL_PREFIX}/bin/oft_mesh_check\" \"\${libs}\" \"\")
      "
      COMPONENT app
    )
    # Fixup runtime dependencies for libraries
    install(
      CODE
      "
      function(patch_library_rpath bin_lib)
        file(GET_RUNTIME_DEPENDENCIES UNRESOLVED_DEPENDENCIES_VAR deps_var LIBRARIES \${bin_lib})
        find_program(CMAKE_INSTALL_NAME_TOOL NAMES install_name_tool HINTS ${_CMAKE_TOOLCHAIN_LOCATION})
        set(changes \"\")
        foreach(pr \${deps_var})
          if(\${pr} MATCHES \"@executable_path*\")
            string(REPLACE \"@executable_path\" \"@loader_path\" pr_new \"\${pr}\")
            set(changes \${changes} \"-change\" \"\${pr}\" \"\${pr_new}\")
          endif()
        endforeach()
        list(LENGTH changes nchanges)
        if( nchanges GREATER_EQUAL 1 )
          set(cmd \${CMAKE_INSTALL_NAME_TOOL} \${changes} \"\${bin_lib}\")
          execute_process(COMMAND \${cmd} RESULT_VARIABLE install_name_tool_result ERROR_QUIET )
        endif()
        # Update code signature
        if( \"${CMAKE_HOST_SYSTEM_PROCESSOR}\" STREQUAL \"arm64\" )
          execute_process(COMMAND \"codesign\" \"--force\" \"-o\" \"linker-signed\" \"-s\" \"-\" \"\${bin_lib}\" )
        endif()
      endfunction()
      # Patch @loader_path for all libraries
      file(GLOB libs \${CMAKE_INSTALL_PREFIX}/bin/*.dylib)
      foreach(lib \${libs})
        patch_library_rpath(\${lib})
      endforeach()
      "
      COMPONENT app
    )
else()
  install(
    CODE
    "
    message(STATUS \"Patching executables\")
    execute_process(
      COMMAND ${Python_EXECUTABLE} ${CMAKE_SOURCE_DIR}/cmake/patch_package.py \${CMAKE_BINARY_DIR}/bin
      WORKING_DIRECTORY \${CMAKE_INSTALL_PREFIX}/bin
    )
    "
    COMPONENT app
  )
endif()
