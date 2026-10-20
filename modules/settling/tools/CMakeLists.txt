# Tools of the settling module. Like the module's programs, they replace the core tools of the same
# name, which are then not built, and go to the same place in the build tree (tools/).
set(
    settling_TOOLS
    changegrid
)

foreach (program ${settling_TOOLS})
    install_channelflow_application(settling_${program} bin)
    target_link_libraries(settling_${program}_app PUBLIC settling)
    set_target_properties(settling_${program}_app PROPERTIES
        OUTPUT_NAME ${program}
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/tools)
endforeach (program)
