# Programs of the settling module. While the module is built alongside the original solver they
# carry a settling_ prefix; at the switch-over they replace the programs of the same name.
set(
    settling_APPS
    settling_simulateflow
    settling_findsoln
    settling_continuesoln
    settling_findeigenvals
)

foreach (program ${settling_APPS})
    install_channelflow_application(${program} bin)
    target_link_libraries(${program}_app PUBLIC settling)
endforeach (program)
