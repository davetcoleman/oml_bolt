cmake_minimum_required(VERSION 2.8.3)
project(ompl_bolt)

# C++ 11
set(CMAKE_CXX_FLAGS "-std=c++11 -Wall ${CMAKE_CXX_FLAGS}")

find_package(catkin REQUIRED COMPONENTS
)

find_package(OMPL REQUIRED)
find_package(Boost REQUIRED)

catkin_package(
  CATKIN_DEPENDS
    ompl
  INCLUDE_DIRS
    src
  LIBRARIES
    ${PROJECT_NAME}
)

###########
## Build ##
###########

include_directories(SYSTEM
  ${Boost_INCLUDE_DIR}
  ${OMPL_INCLUDE_DIRS}
  #${EIGEN_INCLUDE_DIRS}
)

include_directories(
  src
  ${catkin_INCLUDE_DIRS}
  ${OMPL_INCLUDE_DIRS}
)

# Secondary Library
add_library(${PROJECT_NAME}_aux
  src/ompl/tools/bolt/src/BoltStorage.cpp
  src/ompl/tools/bolt/src/EdgeCache.cpp
  src/ompl/tools/bolt/src/Discretizer.cpp
)
# Specify libraries to link a library or executable target against
target_link_libraries(${PROJECT_NAME}_aux
  ${catkin_LIBRARIES}
  ${OMPL_LIBRARIES}
  ${Boost_LIBRARIES}
)

# Library
add_library(${PROJECT_NAME}
  src/ompl/tools/bolt/src/Bolt.cpp
  src/ompl/tools/bolt/src/BoltRetrieveRepair.cpp
  src/ompl/tools/bolt/src/DenseDB.cpp
  src/ompl/tools/bolt/src/SparseDB.cpp
)
# Specify libraries to link a library or executable target against
target_link_libraries(${PROJECT_NAME}
  ${PROJECT_NAME}_aux
  ${catkin_LIBRARIES}
  ${OMPL_LIBRARIES}
  ${Boost_LIBRARIES}
)
