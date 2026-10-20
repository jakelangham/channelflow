#  
# Include and configure the channelflow tools
#
# This file is a part of channelflow version 2.0, https://channelflow.ch .
# License is GNU GPL version 2 or later: ./LICENSE
#
set(tool_APPS
  addfields
#  benchmark
  correlation
#  diffop
#  extrapolatefields
  fieldcompare
  fieldconvert
  fieldprops
  findsymmetries
  L2op
#  optphaseshift
  perturbfield
# pressure
  randomfield
  settlingprofile
#  symmetrize
#  symmetryop
  wall_shear
  )

# With the settling module (the default) its version of changegrid is built instead
if (NOT WITH_SETTLING)
    set(tool_APPS ${tool_APPS} changegrid)
endif ()

foreach (program ${tool_APPS})
    install_channelflow_application(${program} bin)
endforeach ()
