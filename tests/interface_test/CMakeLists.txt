add_executable(test_interface test_interface.c)

target_link_libraries(test_interface
     candock
    ${Boost_LIBRARIES}
    ${OPENMM_LIBRARY}
    ${GSL_LIBRARIES}
)