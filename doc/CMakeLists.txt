configure_file( ${PROJECT_SOURCE_DIR}/doc/README.md
                ${PROJECT_BINARY_DIR}/doc/README.md
                @ONLY
)

add_custom_command( OUTPUT candock_defaults.cfg
                    COMMAND cd_check_options
                    ARGS --help > candock_defaults.cfg
)

add_custom_target(  documentation ALL
                    COMMAND perl ${PROJECT_SOURCE_DIR}/doc/boost_po_to_md_table.pl candock_defaults.cfg >> README.md
                    DEPENDS candock_defaults.cfg README.md
)