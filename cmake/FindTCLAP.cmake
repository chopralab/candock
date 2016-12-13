find_path(TCLAP_INCLUDE_DIR tclap/CmdLine.h
          DOC "Include for TCLAP"
         )

set(TCLAP_INCLUDE_DIRS ${TCLAP_INCLUDE_DIR} )

mark_as_advanced(TCLAP_INCLUDE_DIR)

find_package_handle_standard_args(TCLAP DEFAULT_MSG TCLAP_INCLUDE_DIR)
