cmake_minimum_required(VERSION 3.15)
project(graph_reduce)

set(CMAKE_CXX_STANDARD 14)

set(CPLEX_INCLUDE "/opt/ibm/ILOG/CPLEX_Studio129/cplex/include")
set(CONCERT_INCLUDE "/opt/ibm/ILOG/CPLEX_Studio129/concert/include")
set(CPLEX_LIB "/opt/ibm/ILOG/CPLEX_Studio129/cplex/lib/x86-64_linux/static_pic")
set(CONCERT_LIB "/opt/ibm/ILOG/CPLEX_Studio129/concert/lib/x86-64_linux/static_pic")

include_directories(
        ${CPLEX_INCLUDE}
        ${CONCERT_INCLUDE}
)

link_directories(
        ${CPLEX_LIB}
        ${CONCERT_LIB}
)

# 生成debug版本
SET(CMAKE_BUILD_TYPE "Release")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -DIL_STD -lilocplex -lconcert -lcplex -lm -lpthread -ldl")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -DIL_STD -lilocplex -lconcert -lcplex -lm -lpthread -ldl")

add_executable(gr gr.h gr.cpp)

target_link_libraries(gr ilocplex)
target_link_libraries(gr concert)
target_link_libraries(gr cplex)
# 其他系统库，如果链接失败也一起链接：
TARGET_LINK_LIBRARIES(gr dl)
TARGET_LINK_LIBRARIES(gr pthread)
TARGET_LINK_LIBRARIES(gr m)