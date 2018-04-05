string (REPLACE ";" ":" CANDOCK_ADDITONAL_LD_PATH "${CANDOCK_ADDITONAL_LD_PATH}")

if (APPLE)
  set(ADDITIONAL_LOOKUP "DYLD_LIBRARY_PATH")
else()
  set(ADDITIONAL_LOOKUP "LD_LIBRARY_PATH")
endif()

configure_file( ${PROJECT_SOURCE_DIR}/scripts/load_variables.sh.in
                ${PROJECT_BINARY_DIR}/load_variables.sh
                @ONLY
)

configure_file( ${PROJECT_SOURCE_DIR}/scripts/module.sh.in
                ${PROJECT_BINARY_DIR}/module.sh
                @ONLY
)

configure_file( ${PROJECT_SOURCE_DIR}/scripts/dock_many_proteins.sh.in
                ${PROJECT_BINARY_DIR}/dock_many_proteins.sh
                @ONLY
)

install( PROGRAMS ${PROJECT_BINARY_DIR}/load_variables.sh
                  ${PROJECT_BINARY_DIR}/module.sh
                  ${PROJECT_BINARY_DIR}/dock_many_proteins.sh

                  ${PROJECT_SOURCE_DIR}/scripts/extract_candock_to_scores_and_models.pl
                  ${PROJECT_SOURCE_DIR}/scripts/submit_candock_module.sh
                  ${PROJECT_SOURCE_DIR}/scripts/dock_multiple.sh
                  ${PROJECT_SOURCE_DIR}/scripts/submit_many_proteins.sh
                  ${PROJECT_SOURCE_DIR}/scripts/extract_result.sh
                  
         DESTINATION
                  ${CMAKE_INSTALL_PREFIX}/${CANDOCK_VERSION}/modules
)
