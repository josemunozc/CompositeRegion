CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.5.0 REQUIRED
  HINTS ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR} /opt/dealii/
  )
DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(composite_region)

INCLUDE_DIRECTORIES(/home/zerpiko/git/libraries)

ADD_LIBRARY(mylib /home/zerpiko/git/libraries/DataTools.cpp /home/zerpiko/git/libraries/Names.cpp /home/zerpiko/git/libraries/Material.cpp /home/zerpiko/git/libraries/PorousMaterial.cpp)

#ADD_LIBRARY(mylib /home/zerpiko/git/libraries/Material.cpp /home/zerpiko/git/libraries/PorousMaterial.cpp)

ADD_EXECUTABLE(mycode composite_region.cc)
DEAL_II_SETUP_TARGET(mycode)

TARGET_LINK_LIBRARIES(mycode mylib)

ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

ADD_CUSTOM_TARGET(run
    COMMAND ./mycode input.prm
  )
