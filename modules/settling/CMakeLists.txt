#
# This file is a part of channelflow version 2.0 https://channelflow.ch.
# License is GNU GPL version 2 or later: ./LICENCE
#
# The settling-suspension module: Navier-Stokes coupled to settling scalars (Langham & Hogg 2026).
#
set(
    settling_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/boundarycond.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/densehelmholtz.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/realdensematrix.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/scalarsolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/settlingnse.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/settlingdns.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/settlingdsi.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/settlingstats.cpp
)

set(
    settling_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/boundarycond.h
    ${CMAKE_CURRENT_SOURCE_DIR}/densehelmholtz.h
    ${CMAKE_CURRENT_SOURCE_DIR}/realdensematrix.h
    ${CMAKE_CURRENT_SOURCE_DIR}/scalarsolver.h
    ${CMAKE_CURRENT_SOURCE_DIR}/settlingnse.h
    ${CMAKE_CURRENT_SOURCE_DIR}/settlingdns.h
    ${CMAKE_CURRENT_SOURCE_DIR}/settlingdsi.h
    ${CMAKE_CURRENT_SOURCE_DIR}/settlingstats.h
)

install_channelflow_library(settling)
target_link_fftw(settling PRIVATE)
target_link_libraries(settling PUBLIC chflow)
if (WITH_NSOLVER)
    target_link_libraries(settling PUBLIC nsolver)
endif ()

install(FILES ${settling_HEADERS} DESTINATION include/modules/settling)
