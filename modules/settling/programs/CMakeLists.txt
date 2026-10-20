# Programs of the settling module. They replace the core programs of the same name (simulateflow,
# findsoln, continuesoln, findeigenvals), which are then not built, and go to the same place in the
# build tree (programs/) so scripts keep working.
set(
    settling_APPS
    simulateflow
    findsoln
    continuesoln
    findeigenvals
)

foreach (program ${settling_APPS})
    install_channelflow_application(settling_${program} bin)
    target_link_libraries(settling_${program}_app PUBLIC settling)
    set_target_properties(settling_${program}_app PROPERTIES
        OUTPUT_NAME ${program}
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/programs)
endforeach (program)
