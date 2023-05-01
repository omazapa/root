############################################################################
# CMakeLists.txt file for building ROOT math/ipopt package
############################################################################



ROOT_STANDARD_LIBRARY_PACKAGE(Ipopt
                              HEADERS 
                                Math/IpoptMinimizer.h
                              LINKDEF 
                                Math/LinkDef.h
                              SOURCES
                                src/IpoptMinimizer.cxx
                              LIBRARIES ${IPOPT_LIBRARIES}
                              DEPENDENCIES Core MathCore RIO MathMore)
#target_compile_definitions(Ipopt PUBLIC ${IPOPT_DEFINITIONS})
target_include_directories(Ipopt PUBLIC ${IPOPT_INCLUDE_DIRS})
ROOT_ADD_TEST_SUBDIRECTORY(test)
