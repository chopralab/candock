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
                  ${PROJECT_SOURCE_DIR}/scripts/splitfragments.pl
                  ${PROJECT_SOURCE_DIR}/scripts/submit_candock_module.sh
                  ${PROJECT_SOURCE_DIR}/scripts/dock_multiple.sh
                  ${PROJECT_SOURCE_DIR}/scripts/submit_many_proteins.sh
                  
                  ${PROJECT_SOURCE_DIR}/scripts/find_centroids.sh
                  ${PROJECT_SOURCE_DIR}/scripts/prep_fragments.sh
                  ${PROJECT_SOURCE_DIR}/scripts/make_fragments.sh
                  ${PROJECT_SOURCE_DIR}/scripts/dock_fragments.sh
                  ${PROJECT_SOURCE_DIR}/scripts/link_fragments.sh
                  ${PROJECT_SOURCE_DIR}/scripts/design_ligands.sh
                  ${PROJECT_SOURCE_DIR}/scripts/extract_result.sh
         DESTINATION
                  ${CMAKE_INSTALL_PREFIX}/modules/v${CANDOCK_MAJOR_VERSION}.${CANDOCK_MINOR_VERSION}.${CANDOCK_TWEAK_VERSION}
)
