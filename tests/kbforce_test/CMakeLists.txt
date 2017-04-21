add_executable(test_kbforce TestReferenceExampleForce.cpp)

target_link_libraries(test_kbforce
candock ${Boost_LIBRARIES}
        ${OPENMM_LIBRARY}
        ${GSL_LIBRARIES}
        )
