
find_package(HDF5 REQUIRED)


LIST(APPEND sources "main-eiger.cpp")

#add_definitions(-DCHEETAH_SIT_DATA="${ANA_SIT_DATA}")

include_directories("include")
include_directories(${CHEETAH_INCLUDES})
include_directories(${HDF5_INCLUDE_DIR})

#add_library(cheetah-sacla-1 SHARED ${sources})

add_executable(cheetah-eiger ${sources})

add_dependencies(cheetah-eiger cheetah)

target_link_libraries(cheetah-eiger ${CHEETAH_LIBRARY})

#set_target_properties(
# cheetah-sacla
# PROPERTIES SOVERSION 1
# VERSION 1
#)

install(TARGETS cheetah-eiger
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib${LIB_SUFFIX})
