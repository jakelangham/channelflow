# Unit tests of the settling module (run by ctest and the regression suite's unit case)
set(settling_TESTS settling_packingTest settling_scalarsolverTest)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/settle_eq3.nc DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)

foreach (program ${settling_TESTS})
    install_channelflow_application(${program} OFF)
    target_link_libraries(${program}_app PUBLIC settling)
endforeach (program)

add_serial_test(settling_packing settling_packingTest)
add_serial_test(settling_scalarsolver settling_scalarsolverTest)

if (USE_MPI)
    add_mpi_test(mpi_settling_packing settling_packingTest)
endif ()
