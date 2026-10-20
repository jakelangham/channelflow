#  
# Include and configure the channelflow programs
#
# This file is a part of channelflow version 2.0 https://channelflow.ch.
# License is GNU GPL version 2 or later: ./LICENCE
#
# With the settling module (the default) its programs take these names; the ones here are the
# velocity-only solver.
set(program_APPS
#  edgetracking
  )

if (NOT WITH_SETTLING)
    set(program_APPS ${program_APPS} simulateflow)
    if (WITH_NSOLVER)
        set(program_APPS ${program_APPS}
             continuesoln
             findsoln
             findeigenvals)
    endif ()
endif ()


foreach (program ${program_APPS})
    install_channelflow_application(${program} bin)
endforeach ()
