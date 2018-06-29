PROJECT(mba)

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)


# Include directories

INCLUDE_DIRECTORIES(
  ${mba_SOURCE_DIR}/include
  )

# Make the mba library

FILE(GLOB_RECURSE mba_SRCS src/*.C src/*.cpp include/*.h)
ADD_LIBRARY(mba ${mba_SRCS})

set_target_properties(mba PROPERTIES
  CXX_STANDARD 11)

# Set various compiler flags

IF(NOT WIN32)
  # Set warning level to '-Wall' on Linux
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
ELSE(NOT WIN32)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP8
  -D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS")
ENDIF(NOT WIN32)

# Get rid of MESSAGEs in Release mode
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNVERBOSE")


# Apps, examples, tests, ...?

FILE(GLOB mba_APPS app/*.C app/*.cpp)
FOREACH(app ${mba_APPS})
  GET_FILENAME_COMPONENT(appname ${app} NAME_WE)
  ADD_EXECUTABLE(${appname} ${app})
  TARGET_LINK_LIBRARIES(${appname} mba)  
ENDFOREACH(app)

# Copy data
ADD_CUSTOM_COMMAND(
  TARGET mba
  POST_BUILD
  COMMAND ${CMAKE_COMMAND}
  ARGS -E copy_directory ${mba_SOURCE_DIR}/Data 
  ${mba_BINARY_DIR}/Data
  )


# 'install' target

IF(WIN32)
  # Windows
  SET(CMAKE_INSTALL_PREFIX CACHE INTERNAL "")
  SET(mba_INSTALL_PREFIX "$ENV{PROGRAMFILES}/MBA"
    CACHE PATH "Path to install MBA")
  INSTALL(TARGETS mba DESTINATION ${mba_INSTALL_PREFIX}/lib)
  INSTALL(DIRECTORY include
    DESTINATION ${mba_INSTALL_PREFIX}/include
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ELSE(WIN32)
  # Linux
  # lib
  INSTALL(TARGETS mba DESTINATION lib/mba)
  # include
  INSTALL(DIRECTORY include
    DESTINATION include/mba
    FILES_MATCHING PATTERN "*.h"
    PATTERN ".svn" EXCLUDE
    )
ENDIF(WIN32)


# CPack stuff
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "MBA-1.3.0")
SET(CPACK_SOURCE_IGNORE_FILES "/CVS/;/\\.svn/;\\.swp$;\\.#;/#;.*~")
INCLUDE(CPack)
